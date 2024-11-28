class TrajectoryWriter:
    def __init__(self, mda_universe):
        self.mda_universe = mda_universe

    def write_static_file(self, path, selection='all'):
        atoms = self.mda_universe.select_atoms(selection)
        atoms.write(path)

    def write_trajectory_file(self, path, selection='all'):
        atoms = self.mda_universe.select_atoms(selection)
        atoms.write(path, frames='all')
