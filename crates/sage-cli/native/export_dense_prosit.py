"""Expose the existing trained Dense output; remove only the inference tail."""
import argparse
import hashlib
import json
from pathlib import Path

import onnx
from onnx import TensorProto, helper

parser = argparse.ArgumentParser()
parser.add_argument('source')
parser.add_argument('destination')
args = parser.parse_args()
model = onnx.load(args.source)
nodes = list(model.graph.node)
assert [n.op_type for n in nodes[-3:]] == ['LeakyRelu', 'Reshape', 'Cast']
activation, flatten, cast = nodes[-3:]
assert flatten.input[0] == activation.output[0]
assert cast.input[0] == flatten.output[0]
assert cast.output[0] == model.graph.output[0].name
dense = activation.input[0]
assert nodes[-4].op_type == 'Reshape' and nodes[-4].output[0] == dense
del model.graph.node[:]
model.graph.node.extend(nodes[:-3])
model.graph.node.append(helper.make_node('Identity', [dense], ['dense_intensities']))
del model.graph.output[:]
model.graph.output.append(helper.make_tensor_value_info('dense_intensities', TensorProto.FLOAT, ['batch', 29, 6]))
onnx.checker.check_model(model)
destination = Path(args.destination)
if destination.exists():
    raise FileExistsError(destination)
onnx.save(model, destination)
provenance = dict(source=str(Path(args.source).resolve()),
                  source_sha256=hashlib.sha256(Path(args.source).read_bytes()).hexdigest(),
                  output_sha256=hashlib.sha256(destination.read_bytes()).hexdigest(),
                  removed_nodes=[dict(name=n.name, op=n.op_type) for n in nodes[-3:]],
                  output='dense_intensities', dtype='float32', shape=['batch', 29, 6],
                  trained_initializers_unchanged=True)
destination.with_suffix('.provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')
print(json.dumps(provenance, indent=2))
