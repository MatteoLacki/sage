use serde::{Deserialize, Serialize};

use crate::mass::Tolerance;

/// A piecewise-linear function sampled on an equally spaced grid.
///
/// Outside `[grid_start, grid_start + grid_step * (values.len() - 1)]`,
/// `eval` clamps to the nearest edge value (flat extrapolation) rather than
/// extending the boundary segment's slope.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct LinearSpline {
    pub grid_start: f32,
    pub grid_step: f32,
    pub values: Vec<f32>,
}

impl LinearSpline {
    /// `values` must be non-empty and `grid_step` must be positive.
    pub fn validate(&self) -> Result<(), String> {
        if self.values.is_empty() {
            return Err("LinearSpline.values must not be empty".into());
        }
        if !(self.grid_step > 0.0) {
            return Err(format!(
                "LinearSpline.grid_step must be positive, got {}",
                self.grid_step
            ));
        }
        Ok(())
    }

    /// Linear-interpolate the spline at `x`. Panics if `values` is empty —
    /// call [`LinearSpline::validate`] once at config-load time to rule
    /// that out up front.
    pub fn eval(&self, x: f32) -> f32 {
        assert!(!self.values.is_empty(), "LinearSpline has no grid points");
        if self.values.len() == 1 {
            return self.values[0];
        }

        let last = self.values.len() - 1;
        let pos = (x - self.grid_start) / self.grid_step;
        if pos <= 0.0 {
            return self.values[0];
        }
        if pos >= last as f32 {
            return self.values[last];
        }

        let i = pos.floor() as usize;
        let t = pos - i as f32;
        self.values[i] * (1.0 - t) + self.values[i + 1] * t
    }
}

/// Fragment ppm tolerance as a function of fragment mass, given as two
/// independent linear splines for the lower/left and upper/right edges of
/// the ppm window (asymmetric error distributions are expected, same as
/// `Tolerance::Ppm(lo, hi)`).
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct FragmentTolSpline {
    pub ppm_lo: LinearSpline,
    pub ppm_hi: LinearSpline,
}

impl FragmentTolSpline {
    pub fn validate(&self) -> Result<(), String> {
        self.ppm_lo.validate()?;
        self.ppm_hi.validate()?;
        Ok(())
    }

    pub fn tolerance_at(&self, mass: f32) -> Tolerance {
        Tolerance::Ppm(self.ppm_lo.eval(mass), self.ppm_hi.eval(mass))
    }
}

#[cfg(test)]
mod test {
    use super::*;

    fn spline() -> LinearSpline {
        // grid points at x = 0, 10, 20, 30 -> values 0, 10, 20, 60
        LinearSpline {
            grid_start: 0.0,
            grid_step: 10.0,
            values: vec![0.0, 10.0, 20.0, 60.0],
        }
    }

    #[test]
    fn eval_at_exact_grid_points() {
        let s = spline();
        assert_eq!(s.eval(0.0), 0.0);
        assert_eq!(s.eval(10.0), 10.0);
        assert_eq!(s.eval(20.0), 20.0);
        assert_eq!(s.eval(30.0), 60.0);
    }

    #[test]
    fn eval_interpolates_midpoints() {
        let s = spline();
        assert_eq!(s.eval(5.0), 5.0);
        assert_eq!(s.eval(25.0), 40.0); // halfway between 20 and 60
        assert!((s.eval(21.0) - 24.0).abs() < 1e-4); // 20 + 0.1*(60-20)
    }

    #[test]
    fn eval_clamps_outside_grid_range() {
        let s = spline();
        assert_eq!(s.eval(-100.0), 0.0);
        assert_eq!(s.eval(1000.0), 60.0);
    }

    #[test]
    fn eval_single_value_is_constant() {
        let s = LinearSpline {
            grid_start: 0.0,
            grid_step: 10.0,
            values: vec![7.5],
        };
        assert_eq!(s.eval(-50.0), 7.5);
        assert_eq!(s.eval(0.0), 7.5);
        assert_eq!(s.eval(1234.0), 7.5);
    }

    #[test]
    fn validate_rejects_empty_values() {
        let s = LinearSpline {
            grid_start: 0.0,
            grid_step: 10.0,
            values: vec![],
        };
        assert!(s.validate().is_err());
    }

    #[test]
    fn validate_rejects_non_positive_step() {
        let s = LinearSpline {
            grid_start: 0.0,
            grid_step: 0.0,
            values: vec![1.0, 2.0],
        };
        assert!(s.validate().is_err());

        let s = LinearSpline {
            grid_start: 0.0,
            grid_step: -5.0,
            values: vec![1.0, 2.0],
        };
        assert!(s.validate().is_err());
    }

    #[test]
    fn fragment_tol_spline_tolerance_at() {
        let fts = FragmentTolSpline {
            ppm_lo: LinearSpline {
                grid_start: 0.0,
                grid_step: 10.0,
                values: vec![-5.0, -10.0],
            },
            ppm_hi: LinearSpline {
                grid_start: 0.0,
                grid_step: 10.0,
                values: vec![5.0, 20.0],
            },
        };
        match fts.tolerance_at(5.0) {
            Tolerance::Ppm(lo, hi) => {
                assert_eq!(lo, -7.5);
                assert_eq!(hi, 12.5);
            }
            other => panic!("expected Ppm, got {other:?}"),
        }
    }
}
