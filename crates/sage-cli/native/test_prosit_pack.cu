// Standalone GPU fixture runner. No TensorRT/model or production test hooks.
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <vector>
#include "prosit_pack.cuh"

void checked(cudaError_t status) {
    if (status != cudaSuccess) throw std::runtime_error(cudaGetErrorString(status));
}
template<class T> std::vector<T> read_values(size_t n) {
    std::vector<T> values(n);
    if (!std::cin.read(reinterpret_cast<char*>(values.data()), n * sizeof(T)))
        throw std::runtime_error("incomplete test fixture");
    return values;
}
template<class T> T* upload(const std::vector<T>& values) {
    T* device = nullptr;
    checked(cudaMalloc(&device, values.size() * sizeof(T)));
    checked(cudaMemcpy(device, values.data(), values.size() * sizeof(T), cudaMemcpyHostToDevice));
    return device;
}
int main() {
    try {
        auto header = read_values<uint32_t>(3);
        const bool dense = header[2] & 256;
        const bool ranked_half = header[2] & 512;
        int n = header[0], k = header[1], width = header[2] & 255;
        if (n < 1 || n > 100000 || k < 1 || k > 174 || (width != 1 && width != 2))
            throw std::runtime_error("invalid fixture dimensions");
        void* d_raw = dense ? static_cast<void*>(upload(read_values<float>(n * 174)))
                            : static_cast<void*>(upload(read_values<uint16_t>(n * 174)));
        auto tokens = read_values<int32_t>(n * 30);
        auto charges = read_values<float>(n * 6);
        auto d_tokens = upload(tokens);
        auto d_charges = upload(charges);
        std::vector<uint8_t> result(n * (1 + k * (width + 1)));
        uint8_t* d_result = nullptr;
        checked(cudaMalloc(&d_result, result.size()));
        if (dense) pack_prosit<float, true><<<n, 256>>>(static_cast<float*>(d_raw), d_tokens, d_charges, k, width, d_result);
        else if (ranked_half) pack_prosit<uint16_t, true><<<n, 256>>>(static_cast<uint16_t*>(d_raw), d_tokens, d_charges, k, width, d_result);
        else pack_prosit<<<n, 256>>>(static_cast<uint16_t*>(d_raw), d_tokens, d_charges, k, width, d_result);
        checked(cudaGetLastError());
        checked(cudaDeviceSynchronize());
        checked(cudaMemcpy(result.data(), d_result, result.size(), cudaMemcpyDeviceToHost));
        checked(cudaFree(d_raw));
        checked(cudaFree(d_tokens));
        checked(cudaFree(d_charges));
        checked(cudaFree(d_result));
        std::cout.write(reinterpret_cast<const char*>(result.data()), result.size());
        return std::cout ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
