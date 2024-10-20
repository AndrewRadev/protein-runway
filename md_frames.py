import MDAnalysis as mda


def write_frames(universe, file_prefix, frame_set):
    atoms = universe.select_atoms('protein')

    # Get the max number of digits to format with leading zeroes:
    digit_count = len(str(len(universe.trajectory)))

    for timestep in universe.trajectory:
        frame = timestep.frame

        if frame in frame_set:
            with mda.Writer(file_prefix + f'_{frame:0{digit_count}}.pdb') as w:
                w.write(atoms)


def main():
    u = mda.Universe('data/crambin_md.pdb', 'data/crambin_md.xtc')
    # u = mda.Universe('data/spike_protein_md.pdb', 'data/spike_protein_md.xtc')

    frame_count  = len(u.trajectory)
    last_index   = frame_count - 1
    middle_index = frame_count // 2

    write_frames(u, 'output/crambin_frame', {0, middle_index, last_index})
    # write_frames(u, 'output/spike_protein_frame', {0, middle_index, last_index})


if __name__ == "__main__":
    main()
