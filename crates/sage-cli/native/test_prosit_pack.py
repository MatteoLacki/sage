"""Independent NumPy reference for the actual CUDA kernel.

Build test_prosit_pack.cu with nvcc and pass its executable as argv[1].
No inference server is needed; one CUDA device is required.
"""
import json
import struct
import subprocess
import sys

import numpy as np


def reference(raw, tokens, charges, k, width, ranked=False):
    dense = raw.dtype == np.float32
    result = np.zeros((len(raw), 1 + k * (width + 1)), dtype=np.uint8)
    model_ids = np.arange(174)
    for row, bits in enumerate(raw):
        length = np.count_nonzero(tokens[row])
        precursor = np.flatnonzero(charges[row])
        intensities = bits if dense else bits.view(np.float16).astype(np.float32)
        if not np.all(np.isfinite(intensities)) or length == 0 or len(precursor) == 0:
            result[row, 0] = 255
            continue
        charge = precursor[0] + 1
        valid = (model_ids // 6 + 1 < length) & (model_ids % 3 < charge) & (intensities >= 0)
        ids = model_ids[valid]
        if width == 2 and len(ids) and intensities[ids].max() > 65504:
            result[row, 0] = 255
            continue
        order = np.lexsort((ids, -intensities[ids]))
        ids = ids[order[:k]]
        if not (dense or ranked):
            ids = np.sort(ids)
        count = len(ids)
        values = intensities[ids]
        result[row, 0] = count
        if width == 2:
            result[row, 1:1 + count * 2] = values.astype('<f2').view(np.uint8)
        elif count and values.max() > 0:
            result[row, 1:1 + count] = np.rint(np.sqrt(values / values.max()) * np.float32(255)).astype(np.uint8)
        annotations = np.where(ids % 6 < 3, 87 + (length - 2 - ids // 6) * 3 + ids % 3,
                               ids // 6 * 3 + ids % 3)
        result[row, 1 + k * width:1 + k * width + count] = annotations
    return result


def compare(executable, raw, tokens, charges, k, width, label, ranked=False):
    format_width = width | (256 if raw.dtype == np.float32 else 0) | (512 if ranked else 0)
    payload = struct.pack('<III', len(raw), k, format_width) + raw.tobytes() + tokens.tobytes() + charges.tobytes()
    actual = subprocess.run([executable], input=payload, stdout=subprocess.PIPE, check=True).stdout
    expected = reference(raw, tokens, charges, k, width, ranked=ranked).tobytes()
    if actual != expected:
        a, b = np.frombuffer(actual, np.uint8), np.frombuffer(expected, np.uint8)
        index = np.flatnonzero(a != b)[0]
        raise AssertionError((label, k, width, 'byte', int(index), int(a[index]), int(b[index])))
    print(json.dumps(dict(fixture=label, rows=len(raw), k=k, width=width, exact=True)), flush=True)


if __name__ == '__main__':
    executable = sys.argv[1]
    rng = np.random.default_rng(20261022)
    raw = rng.integers(0, 65536, size=(2048, 174), dtype=np.uint16)
    raw[raw & 0x7c00 == 0x7c00] = 0xbc00
    tokens = np.zeros((len(raw), 30), dtype='<i4')
    charges = np.zeros((len(raw), 6), dtype='<f4')
    for row in range(len(raw)):
        tokens[row, :row % 30 + 1] = 1
        charges[row, row % 6] = 1
    # All ties, signed zeros, empty spectra, above-one values, and nonfinite
    # values in slots that length/charge filtering would otherwise remove.
    raw[:16] = 0xbc00
    raw[0] = 0x3c00
    raw[1] = 0
    raw[2] = 0x8000
    raw[3, [0, 3, 6]] = [0x4800, 0x4000, 0]
    raw[4, [0, 3]] = [0x8000, 0]
    raw[5:9, 173] = [0x7c00, 0xfc00, 0x7e00, 0xfe00]
    tokens[9] = 0
    charges[10] = 0
    raw[11] = 0x0001  # Smallest positive FP16 subnormal.
    for k in [1, 2, 12, 173, 174]:
        for width in [1, 2]:
            compare(executable, raw, tokens, charges, k, width, 'edges-and-random')
            compare(executable, raw, tokens, charges, k, width, 'dense-fp16-ranked', ranked=True)
    # Every FP16 bit pattern, with 8.0 as a competing maximum. Exercises
    # rejection, signs, subnormals, extrema, and quantization rounding.
    raw = np.full((65536, 174), 0xbc00, dtype='<u2')
    raw[:, 0] = np.arange(65536, dtype=np.uint16)
    raw[:, 3] = 0x4800
    tokens = np.zeros((len(raw), 30), dtype='<i4')
    tokens[:, :3] = 1
    charges = np.zeros((len(raw), 6), dtype='<f4')
    charges[:, 0] = 1
    for width in [1, 2]:
        compare(executable, raw, tokens, charges, 12, width, 'all-fp16-patterns')

    raw = rng.normal(0, 4, (2048, 174)).astype('<f4')
    tokens = np.zeros((len(raw), 30), dtype='<i4')
    charges = np.zeros((len(raw), 6), dtype='<f4')
    for row in range(len(raw)):
        tokens[row, :row % 30 + 1] = 1
        charges[row, row % 6] = 1
    raw[0] = 1
    raw[1] = 0
    raw[2] = -0.0
    raw[3] = -1
    raw[3, [0, 3, 6]] = [8, 2, 0]
    raw[4] = -1e-12  # Must be rejected before conversion to FP16 minus zero.
    raw[4, 3] = 0
    raw[5:9, 173] = [np.inf, -np.inf, np.nan, -np.nan]
    tokens[9] = 0
    charges[10] = 0
    raw[11] = np.nextafter(np.float32(0), np.float32(1))
    raw[12] = np.finfo(np.float32).max
    for k in [1, 2, 12, 173, 174]:
        for width in [1, 2]:
            compare(executable, raw, tokens, charges, k, width, 'dense-fp32-ranked')
