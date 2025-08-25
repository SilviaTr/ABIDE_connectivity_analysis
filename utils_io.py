# utils_io.py
import numpy as np, pandas as pd, json
from pathlib import Path

def read_1d(path):
    return pd.read_csv(path, sep=r"\s+", comment="#", header=None)

def save_npz(path, **arrays):
    np.savez_compressed(path, **arrays)

def load_npz(path):
    return np.load(path, allow_pickle=True)

def save_memmap_cube(path, shape, dtype="float32"):
    path = Path(path)
    return np.memmap(path, dtype=dtype, mode="w+", shape=shape)

def open_memmap_cube(path, shape, dtype="float32"):
    return np.memmap(path, dtype=dtype, mode="r", shape=shape)

def dump_json(path, obj):
    Path(path).write_text(json.dumps(obj))

def load_json(path):
    return json.loads(Path(path).read_text())
