#!/usr/bin/env bash
set -euo pipefail
source_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
destination="${1:?usage: build_prosit_bridge.sh /path/to/libprosit_trt.so}"
"${NVCC:-${CUDA_HOME:-/usr/local/cuda}/bin/nvcc}" -x cu -std=c++17 -O3 -shared \
    -Xcompiler=-fPIC,-Wno-deprecated-declarations --fmad=false \
    -I"${TENSORRT_INCLUDE_DIR:-/usr/include/x86_64-linux-gnu}" \
    -I"${CUDA_HOME:-/usr/local/cuda}/include" \
    "$source_dir/prosit_trt.cpp" \
    -L"${TENSORRT_LIBRARY_DIR:-/usr/lib/x86_64-linux-gnu}" \
    -L"${CUDA_HOME:-/usr/local/cuda}/lib64" \
    -lnvinfer -lcudart -Xlinker=-rpath,'$ORIGIN' -o "$destination"
