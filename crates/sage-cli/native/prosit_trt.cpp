// Optional TensorRT 8.6 bridge. No Triton or Python runtime is used.
#include <NvInfer.h>
#include <cuda_runtime_api.h>
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <vector>
#include "prosit_pack.cuh"

namespace {
void check(cudaError_t error) {
    if (error != cudaSuccess) throw std::runtime_error(cudaGetErrorString(error));
}
void set_error(char* buffer, size_t capacity, const char* message) noexcept {
    if (!buffer || capacity == 0) return;
    size_t n = std::min(capacity - 1, std::strlen(message));
    std::memcpy(buffer, message, n);
    buffer[n] = 0;
}
struct Logger : nvinfer1::ILogger {
    std::mutex mutex;
    std::string last_error;
    void log(Severity severity, const char* message) noexcept override {
        if (severity <= Severity::kERROR) {
            std::lock_guard<std::mutex> lock(mutex);
            last_error = message;
        }
    }
};
struct Model {
    Logger logger;
    std::unique_ptr<nvinfer1::IRuntime> runtime;
    std::unique_ptr<nvinfer1::ICudaEngine> engine;
    std::unique_ptr<nvinfer1::IExecutionContext> context;
    cudaStream_t stream = nullptr;
    std::vector<void*> buffers;
    int peptide = -1, charge = -1, energy = -1, output = -1;
    int capacity = 0;
    uint8_t* packed = nullptr;
    bool dense_output = false;
    bool dense_float = false;

    ~Model() {
        if (stream) cudaStreamSynchronize(stream);
        for (void* p : buffers) if (p) cudaFree(p);
        if (packed) cudaFree(packed);
        if (stream) cudaStreamDestroy(stream);
    }

    int input(const char* name, nvinfer1::DataType type, int width) {
        int index = engine->getBindingIndex(name);
        if (index < 0 || !engine->bindingIsInput(index) || engine->getBindingDataType(index) != type ||
            engine->getBindingFormat(index) != nvinfer1::TensorFormat::kLINEAR)
            throw std::runtime_error(std::string("unexpected model input: ") + name);
        auto dims = engine->getBindingDimensions(index);
        auto max = engine->getProfileDimensions(index, 0, nvinfer1::OptProfileSelector::kMAX);
        auto min = engine->getProfileDimensions(index, 0, nvinfer1::OptProfileSelector::kMIN);
        if (dims.nbDims != 2 || dims.d[1] != width || min.d[0] != 1 || max.d[0] < capacity)
            throw std::runtime_error(std::string("incompatible input shape/profile: ") + name);
        return index;
    }

    void load(const char* path, int batch_size) {
        if (batch_size < 1 || batch_size > 1000) throw std::runtime_error("batch size must be 1..1000");
        capacity = batch_size;
        std::ifstream file(path, std::ios::binary | std::ios::ate);
        if (!file) throw std::runtime_error("cannot open TensorRT engine");
        auto size = file.tellg();
        if (size <= 0) throw std::runtime_error("empty TensorRT engine");
        std::vector<char> plan(static_cast<size_t>(size));
        file.seekg(0);
        if (!file.read(plan.data(), size)) throw std::runtime_error("cannot read TensorRT engine");
        runtime.reset(nvinfer1::createInferRuntime(logger));
        if (!runtime) throw std::runtime_error("cannot create TensorRT runtime");
        engine.reset(runtime->deserializeCudaEngine(plan.data(), plan.size()));
        if (!engine) throw std::runtime_error("cannot deserialize engine: " + logger.last_error);
        if (engine->hasImplicitBatchDimension() || engine->getNbBindings() != 4)
            throw std::runtime_error("expected a four-binding explicit-batch Prosit engine");
        peptide = input("peptides_in", nvinfer1::DataType::kINT32, 30);
        charge = input("precursor_charge_in", nvinfer1::DataType::kFLOAT, 6);
        energy = input("collision_energy_in", nvinfer1::DataType::kFLOAT, 1);
        output = engine->getBindingIndex("dense_intensities");
        dense_output = output >= 0;
        if (!dense_output) output = engine->getBindingIndex("intensities");
        dense_float = dense_output && engine->getBindingDataType(output) == nvinfer1::DataType::kFLOAT;
        const auto dtype = dense_float ? nvinfer1::DataType::kFLOAT : nvinfer1::DataType::kHALF;
        if (output < 0 || engine->bindingIsInput(output) || engine->getBindingDataType(output) != dtype ||
            engine->getBindingFormat(output) != nvinfer1::TensorFormat::kLINEAR)
            throw std::runtime_error("expected FP16/FP32 dense_intensities or legacy FP16 intensities output");
        context.reset(engine->createExecutionContext());
        if (!context) throw std::runtime_error("cannot create execution context: " + logger.last_error);
        check(cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking));
        buffers.resize(4, nullptr);
        check(cudaMalloc(&buffers[peptide], capacity * 30 * sizeof(int32_t)));
        check(cudaMalloc(&buffers[charge], capacity * 6 * sizeof(float)));
        check(cudaMalloc(&buffers[energy], capacity * sizeof(float)));
        check(cudaMalloc(&buffers[output], capacity * 174 * (dense_float ? sizeof(float) : sizeof(uint16_t))));
        check(cudaMalloc(&packed, capacity * (1 + 174 * 3)));
    }

    void infer(int n, const int32_t* tokens, const float* charges, const float* energies,
               int k, int encoding, uint8_t* result) {
        if (n < 1 || n > capacity) throw std::runtime_error("invalid inference batch size");
        if (k < 1 || k > 174 || (encoding != 0 && encoding != 1))
            throw std::runtime_error("invalid postprocessing configuration");
        const int width = encoding == 0 ? 2 : 1;
        if (!context->setBindingDimensions(peptide, nvinfer1::Dims2(n, 30)) ||
            !context->setBindingDimensions(charge, nvinfer1::Dims2(n, 6)) ||
            !context->setBindingDimensions(energy, nvinfer1::Dims2(n, 1)))
            throw std::runtime_error("cannot set inference dimensions");
        auto shape = context->getBindingDimensions(output);
        const bool expected_shape = dense_output
            ? shape.nbDims == 3 && shape.d[0] == n && shape.d[1] == 29 && shape.d[2] == 6
            : shape.nbDims == 2 && shape.d[0] == n && shape.d[1] == 174;
        if (!expected_shape)
            throw std::runtime_error("unexpected output shape");
        check(cudaMemcpyAsync(buffers[peptide], tokens, n * 30 * sizeof(int32_t), cudaMemcpyHostToDevice, stream));
        check(cudaMemcpyAsync(buffers[charge], charges, n * 6 * sizeof(float), cudaMemcpyHostToDevice, stream));
        check(cudaMemcpyAsync(buffers[energy], energies, n * sizeof(float), cudaMemcpyHostToDevice, stream));
        if (!context->enqueueV2(buffers.data(), stream, nullptr))
            throw std::runtime_error("TensorRT enqueue failed: " + logger.last_error);
        if (dense_float) {
            pack_prosit<float, true><<<n, 256, 0, stream>>>(static_cast<const float*>(buffers[output]),
                static_cast<const int32_t*>(buffers[peptide]), static_cast<const float*>(buffers[charge]),
                k, width, packed);
        } else if (dense_output) {
            pack_prosit<uint16_t, true><<<n, 256, 0, stream>>>(static_cast<const uint16_t*>(buffers[output]),
                static_cast<const int32_t*>(buffers[peptide]), static_cast<const float*>(buffers[charge]),
                k, width, packed);
        } else {
            pack_prosit<<<n, 256, 0, stream>>>(static_cast<const uint16_t*>(buffers[output]),
                static_cast<const int32_t*>(buffers[peptide]), static_cast<const float*>(buffers[charge]),
                k, width, packed);
        }
        check(cudaGetLastError());
        check(cudaMemcpyAsync(result, packed, n * (1 + k * (width + 1)), cudaMemcpyDeviceToHost, stream));
        check(cudaStreamSynchronize(stream));
    }
};
}

extern "C" uint32_t sage_prosit_abi_version() { return 3; }
extern "C" uint32_t sage_prosit_output_kind(void* handle) noexcept {
    if (!handle) return 0;
    const auto* model = static_cast<Model*>(handle);
    return model->dense_output ? (model->dense_float ? 2 : 1) : 0;
}
extern "C" void* sage_prosit_create(const char* path, int capacity, char* error, size_t error_size) noexcept {
    try {
        auto model = std::make_unique<Model>();
        model->load(path, capacity);
        return model.release();
    } catch (const std::exception& e) { set_error(error, error_size, e.what()); }
    catch (...) { set_error(error, error_size, "unknown native initialization error"); }
    return nullptr;
}
extern "C" int sage_prosit_infer(void* handle, int n, const int32_t* tokens, const float* charges,
                                  const float* energies, int k, int encoding, uint8_t* output,
                                  char* error, size_t error_size) noexcept {
    try {
        if (!handle || !tokens || !charges || !energies || !output) throw std::runtime_error("null inference argument");
        static_cast<Model*>(handle)->infer(n, tokens, charges, energies, k, encoding, output);
        return 0;
    } catch (const std::exception& e) {
        // Rust may recycle or drop the host buffers as soon as this call returns,
        // including on errors after an asynchronous copy has been submitted.
        if (handle) cudaStreamSynchronize(static_cast<Model*>(handle)->stream);
        set_error(error, error_size, e.what());
    } catch (...) {
        if (handle) cudaStreamSynchronize(static_cast<Model*>(handle)->stream);
        set_error(error, error_size, "unknown native inference error");
    }
    return -1;
}
extern "C" void sage_prosit_destroy(void* handle) noexcept { delete static_cast<Model*>(handle); }
