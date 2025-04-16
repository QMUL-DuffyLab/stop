STOP: Simulated TCSPC On Proteins
=================================

do you like to shoot lasers at proteins? have you ever thought that maybe the whole process of getting hold of a protein and a laser in real life to do that was a bit too much effort? well, do i have the software for you! it simulates the process of firing an expensive laser at your weird little protein, and you can do it from the comfort of your own desk or couch or wherever you like, you wastrel.

on a more serious note, this is some fortran code to simulate TCSPC experiments on arbitrary proteins (with some constraints, which i'll go into). the idea is that you set up your protein using JSON or python, there's some python code to write all the relevant parameters to fortran-friendly files and then an MPI fortran kernel to do the heavy lifting. once the simulated experiments are done, reconvolution fits are performed to get amplitude-weighted lifetimes and so on.

Installation
============

You'll need a relatively modern version of Python (I think >3.8 should do) with numpy, scipy, matplotlib and sympy, as well as a modern Fortran compiler with OpenMPI. I'm using GCC 14.2.1 on Linux and it works fine but all the Fortran is strictly F2008 so any modern Fortran compiler should work.

Setup
=====

The Python code uses two JSON files, one containing protein parameters and one containing simulation parameters; examples are included in `protein.json` and `simulation.json`.

Protein parameters
------------------

The JSON file is structured as a list of named proteins each with a set of parameters:
```JSON
{
    "protein_1" : {
      {protein 1 parameters}
    },
    "protein_2" : {
      {protein 2 parameters}
    },
    ...
}
```

For a given protein we have:
```JSON
"protein_name" : {
    "n_p": n_p, # the number of pigments in the protein
    "n_s": n_s, # the number of total pigment states
    "pigment_names": [pigment_names], # names of the pigments
    "state_names": [state_names], # names of the states
    "which_pigment": [which_pigment], # which pigment is each state on
    "abundance": [abundances], # are all the states present on all proteins?
    "dist": [[dist]], # which states are distinguishable?
    "n_tot": [n_tot], # total number of each pigment in the protein
    "n_thermal": [n_thermal], # thermally accessible number of each pigment
    "hop": [hop], # 1-d array of hopping rates between neighbouring proteins for each state
    "intra": [[intra]], # 2-d array of intra-protein rates
    "ann": [[ann]], # 2-d array of annihilation rates
    "ann_remainder": [[ann_remainder]], # which state is left behind after annihilation
    "xsec": [xsec], # cross-sections of each state
    "emissive": [emissive] # are the decays of each state emissive
}
```

A few things to note:
- all of the numbers here are given as times in SI units; the Fortran inverts them to get rates.
- the reason for the names is to generate columns for the histogram that the Fortran outputs, to make it more obvious which column corresponds to which thing.
- for `intra`, the convention I've used is that off-diagonal elements denote transition times between different states in the protein, and the diagonal elements indicate decay times.
- `n_tot` and `n_thermal` should be of length `n_p` since they refer to the actual pigments; for example, there are 42 chlorophylls in an LHCII trimer, of which 24 (the Chl *a*s) are thermally accessible. We use these to check stimulated emission and intra-protein transfer rates.
- `abundance` doesn't do anything yet - I haven't added the logic for it. will update when that's done.
- `ann_remainder` is set up the way it is because it is not always trivial which state will be left after an annihilation event; this removes any uncertainty by making you figure it out and write it down, rather than me writing some weird bit of code to just have a guess for you.

Check the included file for examples, as mentioned above.

Simulation parameters
---------------------

This is a simpler little file. I think that most if not all of these should be self-explanatory, but just in case:

- `fwhm` is the FWHM of your pulse (again in SI units)
- `fluence` is in units of photons per pulse per square centimetre ($ \gamma \text{ pulse }^{-1} cm^{-2} $)
- `n_sites` is what it sounds like - it shouldn't make much difference really, but it's there
- `lattice` can be either "hex", "square" or "line". unless you have a good reason to think your protein collection is specifically a line, you can probably just leave this.
- `rep_rate` is the rep rate in hertz
- `tmax` is the maximum time to bin for. Internally we use one time step `dt1` up to `tmax`, then another (longer, preferably) time step `dt2` up until the next pulse starts, with the necessary interval calculated using `rep_rate` given above.
- `binwidth` is the binwidth of the histogram
- `n_counts` is how many counts you want. really shouldn't have to explain this one
- `n_repeats` is how many repeats you want to do

Note: The fitting code will output various CSV and text files containing the fitted arrays, details of the fits and errors and so on, but it doesn't do anything sophisticated to compare the fit between repeats.

Running
=======

Call it from a command line with between one and three arguments. They're documented in `main.py`; you can do `python main.py -h` to see details. The compulsory first argument is the name of the protein you want to simulate, and should be one of the names in `protein.json`. Optional ones are an output path, which if not given will default to `./out` (i.e. a new folder named "out" in your current directory). Third is the number of cores you want to run it on - if you don't give this one it'll do `os.cpu_count()` and just use that, which may not be optimal. The Python will create the output path, make Fortran-friendly parameter files and copy them to the output directory, run `make all` on the Fortran if necessary, run it for you, and then perform reconvolution fits when the Fortran returns. By default it will use the decay times of the states as a starting point and try 1- to n-exponential fits based on those.

TODO: make `tau_init` a parameter.


FAQS
====

Q: where is the test suite, the continuous integration, all that kind of stuff?  
A: what are you, a coward?
