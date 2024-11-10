"""
Input:
- <protein.nmd>:     A file with normal modes to create a trajectory from
- <protein.nmd.pdb>: The output filename to write
"""

import sys
import os
import itertools
from pathlib import Path

import MDAnalysis as mda
from MDAnalysis.coordinates.memory import MemoryReader as MDAMemoryReader
import numpy as np

FRAME_COUNT     = 100
MAGNITUDE_SCALE = 0.1
MODE_COUNT      = 1

if len(sys.argv) < 3:
    print(f"USAGE: python {os.path.basename(__file__)} <protein.nmd> <protein.nmd.pdb>")
    sys.exit(1)

def main():
    nmd_file    = sys.argv[1]
    output_file = sys.argv[2]

    coordinates = None
    atomnames   = None
    resnames    = None
    resids      = None
    modes       = []

    # NMD structure:
    #
    # atomnames CA CA CA ...
    # resnames SER ARG LEU ...
    # resids 0 1 2 3 4 5 6 7 8 11 12 ... <- Note: possible to have gaps
    # coordinates 54.260 50.940 73.060 ...
    # mode 1 29.61 -0.008 -0.005 0.012 ...
    # mode 2 18.28 -0.009 0.004 -0.008 ...
    # mode 3 17.50 -0.027 -0.025 0.023 ...
    # mode <N> <magnitude> <x1> <y1> <z1> <x2> <y2> <z2> ...

    with open(nmd_file) as f:
        for line in f:
            line = line.strip()
            section, _, line = line.partition(' ')

            if section == 'atomnames':
                atomnames = np.array(line.split(' '))

            elif section == 'resnames':
                resnames = np.array(line.split(' '))

            elif section == 'resids':
                resids = np.array(line.split(' '))

            elif section == 'coordinates':
                coordinates = (float(c) for c in line.split(' '))
                coordinates = np.array(list(batched(coordinates, n=3, strict=True)))

            elif section == 'mode':
                mode,      _, line = line.partition(' ')
                magnitude, _, line = line.partition(' ')
                vector_coordinates = line.split(' ')

                magnitude          = float(magnitude) * MAGNITUDE_SCALE
                # magnitude          = 1.0
                vector_coordinates = (magnitude * float(vc) for vc in vector_coordinates)
                vectors            = np.array(list(batched(vector_coordinates, n=3, strict=True)))

                modes.append((mode, vectors))

    for mode, vectors in modes:
        assert(vectors.shape == coordinates.shape)

    trajectory = []

    for _ in range(0, FRAME_COUNT // 2):
        for mode, vectors in modes[:MODE_COUNT]:
            coordinates = coordinates + vectors
        trajectory.append(coordinates)

    for _ in range(0, FRAME_COUNT // 2):
        for mode, vectors in modes[:MODE_COUNT]:
            coordinates = coordinates - vectors
        trajectory.append(coordinates)

    n_atoms = coordinates.shape[0]

    u = mda.Universe.empty(
        n_atoms=n_atoms,
        n_residues=n_atoms,
        n_frames=FRAME_COUNT,
        atom_resindex=np.arange(n_atoms),
        trajectory=True,
    )

    assert(len(atomnames) == len(resids))
    assert(len(resids) == len(resnames))

    u.add_TopologyAttr('names', atomnames)
    u.add_TopologyAttr('resids', resids)
    u.add_TopologyAttr('resnames', resnames)
    # u.add_TopologyAttr('tempfactors', tempfactors)

    u.load_new(np.array(trajectory), format=MDAMemoryReader)

    u.atoms.write(output_file, frames='all')


# See <https://docs.python.org/3/library/itertools.html#itertools.batched>
def batched(iterable, n, *, strict=False):
    # batched('ABCDEFG', 3) → ABC DEF G
    if n < 1:
        raise ValueError('n must be at least one')

    iterator = iter(iterable)

    while batch := tuple(itertools.islice(iterator, n)):
        if strict and len(batch) != n:
            raise ValueError('batched(): incomplete batch')
        yield batch


if __name__ == "__main__":
    main()
