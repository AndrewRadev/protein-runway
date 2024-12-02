from pathlib import Path

class TrajectoryWriter:
    def __init__(self, mda_universe):
        self.mda_universe = mda_universe

    def write_static_file(self, path: Path|str, selection='all'):
        atoms = self.mda_universe.select_atoms(selection)
        atoms.write(str(path))

    def write_trajectory_file(self, path: Path|str, selection='all'):
        atoms = self.mda_universe.select_atoms(selection)
        atoms.write(str(path), frames='all')
