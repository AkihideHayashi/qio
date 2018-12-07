import numpy as np
from . import poscar


def read_chgcar(f):
    system, unit, cell, elem_num, coordinates, binds = poscar.read_poscar(f)
    next(f)  # blank line
    grid_shape = np.array([int(w) for w in next(f).split()])
    grid = []
    for _ in range(int(np.ceil(np.prod(grid_shape) / 5))):
        grid.extend([float(w) for w in next(f).split()])
    grid = np.array(grid).reshape(grid_shape)
    augment = [read_augment(f) for _ in range(coordinates.shape[0])]
    return system, unit, cell, elem_num, coordinates, binds, grid, augment


def read_augment(f):
    n = int(next(f).split()[-1])
    aug = []
    for _ in range(int(np.ceil(n / 5))):
        aug.extend([float(w) for w in next(f).split()])
    return np.array(aug)

