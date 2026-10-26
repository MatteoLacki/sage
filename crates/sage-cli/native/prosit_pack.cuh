// One block per spectrum. Dense outputs use intensity-rank order; legacy FP16
// outputs retain model order for byte compatibility with existing libraries.
#include <cuda_fp16.h>
#include <cstdint>

__device__ float intensity_value(uint16_t bits) { return __half2float(__ushort_as_half(bits)); }
__device__ float intensity_value(float value) { return value; }

template<class T, bool Ranked = false>
__global__ void pack_prosit(const T* raw, const int32_t* tokens,
                           const float* charges, int k, int width, uint8_t* packed) {
    const int row = blockIdx.x;
    const int id = threadIdx.x;
    uint8_t* result = packed + row * (1 + k * (width + 1));
    for (int i = id; i < 1 + k * (width + 1); i += blockDim.x) result[i] = 0;
    __shared__ float values[174];
    __shared__ bool valid[174], selected[174];
    __shared__ int length, charge, bad;
    __shared__ float maximum;
    if (id == 0) {
        length = 0;
        while (length < 30 && tokens[row * 30 + length] != 0) ++length;
        charge = 0;
        for (int z = 0; z < 6; ++z)
            if (charges[row * 6 + z] == 1.0f) charge = z + 1;
        bad = length < 1 || charge == 0;
        maximum = 0.0f;
        for (int i = 0; i < 174; ++i) {
            values[i] = intensity_value(raw[row * 174 + i]);
            if (!isfinite(values[i])) bad = 1;
            valid[i] = i / 6 + 1 < length && i % 3 + 1 <= charge && values[i] >= 0.0f;
            if (valid[i]) maximum = fmaxf(maximum, values[i]);
            selected[i] = false;
        }
    }
    __syncthreads();
    if (bad) {
        if (id == 0) result[0] = 255;
        return;
    }
    if (width == 2 && maximum > 65504.0f) {
        if (id == 0) result[0] = 255;
        return;
    }
    int rank = 0;
    if (id < 174 && valid[id]) {
        const float magnitude = values[id];
        for (int j = 0; j < 174; ++j) {
            const float other = values[j];
            if (valid[j] && (other > magnitude || (other == magnitude && j < id))) ++rank;
        }
        selected[id] = rank < k;
    }
    __syncthreads();
    if (id == 0) {
        int count = 0;
        for (int j = 0; j < 174; ++j) count += selected[j];
        result[0] = static_cast<uint8_t>(count);
    }
    if (id < 174 && selected[id]) {
        int offset = rank;
        if constexpr (!Ranked) {
            offset = 0;
            for (int j = 0; j < id; ++j) offset += selected[j];
        }
        if (width == 1) {
            const float value = values[id];
            // Explicit correctly rounded operations match the FP32 CPU codec.
            const float scaled = maximum == 0.0f ? 0.0f :
                __fmul_rn(255.0f, __fsqrt_rn(__fdiv_rn(value, maximum)));
            result[1 + offset] = static_cast<uint8_t>(__float2uint_rn(scaled));
        } else {
            const uint16_t bits = __half_as_ushort(__float2half_rn(values[id]));
            result[1 + 2 * offset] = bits & 255;
            result[2 + 2 * offset] = bits >> 8;
        }
        const int ordinal = id / 6 + 1;
        const int z = id % 3;
        const int annotation = id % 6 < 3 ? 87 + (length - 1 - ordinal) * 3 + z
                                         : (ordinal - 1) * 3 + z;
        result[1 + k * width + offset] = static_cast<uint8_t>(annotation);
    }
}
