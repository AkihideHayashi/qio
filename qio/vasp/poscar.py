import numpy as np
from collections import OrderedDict


def read_poscar(f):
    def read_switch(f):
        line = next(f).strip()
        selective = True if line[0] in ('S', 's') else False
        if selective:
            line = next(f).strip()
        direct = False if line[0] in ('C', 'c', 'K', 'k') else True
        return selective, direct

    def read_tf(tf):
        if tf in ('T', 't'):
            return True
        elif tf in ('F', 'f'):
            return False
        else:
            raise RuntimeError("{} encountered in read_tf".format(tf))

    def read_coordinates(f, n, selective, direct, unit, cell):
        lines = [next(f).split() for _ in range(n)]
        coordinates = np.array([[float(w) for w in line[:3]] for line in lines])
        if selective:
            bind = [[read_tf(w) for w in line[3:6]] for line in lines]
        else:
            bind = [[True for _ in range(3)] for _ in lines]
        if direct:
            return np.dot(coordinates, cell), bind
        else:
            return coordinates * unit, bind

    system = next(f).strip()
    unit = float(next(f).strip())
    cell = np.array([[
        float(w) for w in next(f).split()] for _ in range(3)]) * unit
    elems = next(f).split()
    nums = np.array([int(w) for w in next(f).split()])
    selective, direct = read_switch(f)
    coordinates, binds = read_coordinates(
            f, sum(nums), selective, direct, unit, cell)
    return system, unit, cell, OrderedDict(zip(elems, nums)), coordinates, binds


def write_poscar(f, system, unit, cell, elem_num, coordinates, binds, direct):
    def write_tf(tf):
        if tf:
            return "T"
        else:
            return "F"

    f.write("{}\n".format(system))
    f.write("{:24.20}\n".format(unit))
    for l in cell / unit:
        f.write("{:< 024.20}  {:< 024.20}  {:< 024.20}\n".format(*l))
    f.write("  ".join(elem_num.keys()) + "\n")
    f.write("  ".join(map(str, elem_num.values())) + "\n")
    f.write("Selective Dynamics\n")
    if direct:
        f.write("Direct\n")
        for c, b in zip(coordinates / unit, binds):
            f.write("{:< 024.20} {:< 024.20} {:< 024.20} {} {} {}\n".format(
                *nl.solve(cell.T, c), *map(write_tf, b)))
    else:
        f.write("Cartesian\n")
        for c, b in zip(coordinates / unit, binds):
            f.write("{:< 024.20} {:< 024.20} {:< 024.20} {} {} {}\n".format(
                *c, *map(write_tf, b)))
