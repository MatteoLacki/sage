"""Build with TensorRT 8.6 on the target GPU, using the old batch profile."""
import argparse
import json
from pathlib import Path

import tensorrt as trt

parser = argparse.ArgumentParser()
parser.add_argument('onnx')
parser.add_argument('destination')
parser.add_argument('--reference-engine', required=True)
parser.add_argument('--output-dtype', choices=['fp16', 'fp32'], default='fp16')
args = parser.parse_args()
destination = Path(args.destination)
if destination.exists():
    raise FileExistsError(destination)
logger = trt.Logger(trt.Logger.INFO)
trt.init_libnvinfer_plugins(logger, '')
runtime = trt.Runtime(logger)
reference = runtime.deserialize_cuda_engine(Path(args.reference_engine).read_bytes())
assert reference is not None
shapes = {reference.get_binding_name(i): reference.get_profile_shape(0, i)
          for i in range(reference.num_bindings) if reference.binding_is_input(i)}
del reference
builder = trt.Builder(logger)
network = builder.create_network(1 << int(trt.NetworkDefinitionCreationFlag.EXPLICIT_BATCH))
reader = trt.OnnxParser(network, logger)
if not reader.parse(Path(args.onnx).read_bytes()):
    raise RuntimeError('\n'.join(str(reader.get_error(i)) for i in range(reader.num_errors)))
assert network.num_outputs == 1
network.get_output(0).dtype = trt.float16 if args.output_dtype == 'fp16' else trt.float32
config = builder.create_builder_config()
config.set_flag(trt.BuilderFlag.FP16)
config.set_memory_pool_limit(trt.MemoryPoolType.WORKSPACE, 1024 ** 3)
config.profiling_verbosity = trt.ProfilingVerbosity.DETAILED
profile = builder.create_optimization_profile()
for name, dimensions in shapes.items():
    profile.set_shape(name, *dimensions)
assert bool(profile), 'invalid optimization profile'
config.add_optimization_profile(profile)
serialized = builder.build_serialized_network(network, config)
assert serialized is not None, 'TensorRT build failed'
destination.write_bytes(bytes(serialized))
engine = runtime.deserialize_cuda_engine(serialized)
inspection = engine.create_engine_inspector().get_engine_information(trt.LayerInformationFormat.JSON)
destination.with_suffix('.layers.json').write_text(inspection)
summary = dict(tensorrt=trt.__version__, profile={k:[list(d) for d in v] for k,v in shapes.items()},
               output_dtype=args.output_dtype,
               fp16_tactics_enabled=True, workspace_limit_bytes=1024**3,
               bindings=[dict(name=engine.get_binding_name(i), shape=list(engine.get_binding_shape(i)),
                              dtype=str(engine.get_binding_dtype(i))) for i in range(engine.num_bindings)])
destination.with_suffix('.build.json').write_text(json.dumps(summary, indent=2) + '\n')
print(json.dumps(summary, indent=2))
