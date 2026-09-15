import argparse
import os
import shutil
import subprocess
import json
import itertools
import numpy as np
import fit
import parse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
            description="set up aggregate simulation",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # required arguments
    parser.add_argument('-pf', '--protein_file', type=str,
            default='protein.json',
            help=r'File to load protein data from in JSON format')
    parser.add_argument('-sf', '--simulation_file', type=str,
            default='simulation.json',
            help=r'File to load protein data from in JSON format')
    # optional arguments
    parser.add_argument('-p', '--protein_name', type=str, default=None,
            help=r'Name of protein within protein file, if there are multiple')
    parser.add_argument('-o', '--outdir', type=str, default='out',
            help="Output directory (default: 'out')")
    parser.add_argument('-c', '--connected', type=bool, default=False,
            help="Toggle connectedness (if False, prevent excitation hopping)")
    parser.add_argument('-n', '--n_procs', type=int, default=0,
            help="Number of MPI processes to use (default is whatever os.cpu_count() returns")

    args = parser.parse_args()

    with open(args.protein_file, "r") as f:
        protein_json = json.load(f)

    with open(args.simulation_file, "r") as f:
        simulation_json = json.load(f)

    if args.protein_name is not None:
        protein_name = args.protein_name
        if args.protein_name in protein_json:
            protein = protein_json[args.protein_name]
        else:
            raise KeyError("Invalid protein choice. "
            f"Add it to {args.protein_file} or choose a different file.")
    else:
        # if there's only one protein in this file use that and warn
        if len(protein_json) == 1:
            k = list(protein_json.keys())[0]
            protein = protein_json[k]
            print(f"Warning: no protein name given; using protein name {k}.")
            print(f"Protein parameters: {protein}")
        else:
            # try using the filename as a key
            k = os.path.splitext(os.path.basename(args.protein_file))[0]
            if k in protein_json:
                protein = protein_json[k]
                print(f"Warning: no protein name given; using protein data {k}.")
                print(f"Protein parameters: {protein}")
            else:
                raise KeyError(f"Couldn't find a valid set of protein data. Check options -pf and -p.")
        protein_name = k

    outdir = parse.generate_output_dirs(protein_json, simulation_json,
             protein_name, args.connected, args.outdir)

    print("Setting up the input files for the fortran...")
    protein_file, simulation_file = parse.write_fortran_inputs(protein_json,
            simulation_json, protein_name, outdir)

    print("Running make on the fortran...")
    subprocess.run(['make', 'all'], check=True)

    if args.n_procs == 0:
        n_procs = os.cpu_count()
        print(f"Number of MPI processes not given. Using n = {n_procs}")
    else:
        n_procs = args.n_procs

    # oh go on then let's do a little experiment shall we
    print("Running...")
    subprocess.run(['mpirun',
        '-np', f"{n_procs}",
        './stop', protein_file, simulation_file], check=True)

    tau_init = [protein["intra"][0][0], np.min(protein["ann"]), 500e-12]
    for i in range(simulation_json["n_repeats"]):
        hist_file = os.path.join(outdir, f"{args.protein}_run_{i + 1:1d}.csv")
        for j in range(len(tau_init)):
            stuff = fit.do_fit(hist_file, tau_init[:j + 1], "simulation.json", None)
