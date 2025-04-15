import argparse
import os
import shutil
import subprocess
import json
import itertools
import numpy as np
import fit

def write_fortran_inputs(pj, sj, name, outdir):
    p = pj[name] # get the JSON data for our specific protein
    print(p)
    pf = os.path.join(outdir, "protein_params")
    with open(pf, 'w') as f:
        f.write(f"{name}\n")
        f.write(f"{p['n_p']}\n")
        f.write(f"{p['n_s']}\n")
        f.write(" ".join(p["pigment_names"]) + "\n")
        f.write(" ".join(p["state_names"]) + "\n")
        f.write(" ".join([f"{i:4d}" for i in p["which_pigment"]]) + "\n")
        f.write(" ".join([f"{i:.4e}" for i in p["abundance"]]) + "\n")
        f.write(" ".join([str(i) for j in p["dist"] for i in j]) + "\n")
        f.write(" ".join([f"{i:4d}" for i in p["n_tot"]]) + "\n")
        f.write(" ".join([f"{i:4d}" for i in p["n_thermal"]]) + "\n")
        f.write(" ".join([f"{i:.4e}" for i in p["hop"]]) + "\n")
        f.write(" ".join([f"{i:.4e}" for j in p["intra"] for i in j]) + "\n")
        f.write(" ".join([f"{i:.4e}" for j in p["ann"] for i in j]) + "\n")
        f.write(" ".join([f"{i:4d}" for j in p["ann_remainder"] for i in j]) + "\n")
        f.write(" ".join([f"{i:.4e}" for i in p["xsec"]]) + "\n")
        f.write(" ".join([str(e) for e in p["emissive"]]) + "\n")

    sf = os.path.join(outdir, "simulation_params")
    with open(sf, 'w') as f:
        f.write(f"{sj['fwhm']}\n")
        f.write(f"{sj['fluence']}\n")
        f.write(f"{sj['n_sites']}\n")
        f.write(f"{sj['lattice']}\n")
        f.write(f"{sj['rep_rate']}\n")
        f.write(f"{sj['tmax']}\n")
        f.write(f"{sj['dt1']}\n")
        f.write(f"{sj['dt2']}\n")
        f.write(f"{sj['binwidth']}\n")
        f.write(f"{sj['n_counts']}\n")
        f.write(f"{sj['n_repeats']}\n")
        f.write(f"{outdir}\n")
    return pf, sf

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
            description="set up aggregate simulation",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # required arguments
    parser.add_argument('-p', '--protein', type=str, required=True,
            help=r'Choice of protein (check protein.json)')
    # optional arguments
    parser.add_argument('-o', '--outdir', type=str, default='out',
            help="Output directory (default: 'out')")
    parser.add_argument('-n', '--n_procs', type=int, default=0,
            help="Number of MPI processes to use (default is whatever os.cpu_count() returns")

    args = parser.parse_args()

    with open("protein.json", "r") as f:
        protein_json = json.load(f)

    with open("simulation.json", "r") as f:
        simulation_json = json.load(f)

    if args.protein in protein_json:
        protein = protein_json[args.protein]
    else:
        raise KeyError("Invalid protein choice. Add it to protein.json!")

    # fortran will not know about OS directory separators, so
    # just add an extra one at the end of the path for it
    fstr = np.format_float_scientific(simulation_json['fluence'])
    outdir = os.path.join(args.outdir,
            f"{args.protein}_{fstr}", "")
    os.makedirs(outdir, exist_ok=True)

    # setup the files
    print("Setting up the input files for the fortran...")
    protein_file, simulation_file = write_fortran_inputs(protein_json,
            simulation_json, args.protein, outdir)


    # check the fortran is compiled and up to date
    print("Running make on the fortran...")
    subprocess.run(['make', 'all'], check=True)

    if args.n_procs == 0:
        n_procs = os.cpu_count()
        print(f"Number of MPI processes not given. Using {n_procs}")
    else:
        n_procs = args.n_procs

    # oh go on then let's do a little experiment shall we
    print("Running...")
    subprocess.run(['mpirun',
        '-np', f"{n_procs}",
        './stop', protein_file, simulation_file], check=True)

    for i in range(simulation_json["n_repeats"]):
        hist_file = os.path.join(outdir, f"{args.protein}_run_{i + 1:1d}.csv")
        stuff = fit.do_fit(hist_file, [4.0e-9], "simulation.json", None)
        stuff = fit.do_fit(hist_file, [0.5e-9, 4.0e-9], "simulation.json", None)
