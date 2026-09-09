STOP: Simulated TCSPC On Proteins
=================================

do you like to shoot lasers at proteins? have you ever thought that maybe the whole process of getting hold of a protein and a laser in real life to do that was a bit too much effort? well, do i have the software for you! it simulates the process of firing an expensive laser over and over again at your poor little protein, and you can do it from the comfort of your own desk or couch or wherever you like.

the idea is that you figure out the details of your protein, what it's made up of and how you think energy transfer works. Then you encode that into a protein parameter file either by writing some JSON or using a GUI, and then specify some simulation (experiment) parameters like the laser fluence and rep rate and so on. There's some python code to write all those relevant parameters to fortran-friendly files, and then an MPI fortran kernel to do the heavy lifting. once the simulated experiments are done, reconvolution fits are performed to get amplitude-weighted lifetimes and so on.

Installation
============

You'll need a modern version of Python (I'm using 3.14) with numpy, scipy, matplotlib, pandas, sympy, as well as pyqt6 if you want to use the GUI and a modern Fortran compiler with OpenMPI to compile the fortran kernel. I'm using GCC 14.2.1 on Linux and it works fine, but all the Fortran is strictly F2008 so any modern Fortran compiler should work. It also all works on WSL; I have not and will not figure out how to get it to work natively on windows so do not ask me.

The easiest way to ensure you have everything you need is to first install miniforge (see [here](https://github.com/conda-forge/miniforge?tab=readme-ov-file#install)), clone this repo and then do
```
mamba create -f environment.yml
mamba activate stop
```
from the root directory of the cloned repo.
Note that e.g. ipython and jupyter are not included in that YAML file because they're not strictly necessary so you'll have to do `mamba install ipython` or similar if you want those.

Setup
=====

The Python code uses two JSON files, one containing protein parameters and one containing simulation parameters; examples are included in `protein.json` and `simulation.json`. You can set these up manually if you'd like; skip to the [protein](#protein-parameters) and [simulation](#simulation-parameters) explain how the JSON is structured. Alternatively, run `python gui.py` from the root directory of the repo. This will take you through the process of creating/importing/modifying the relevant input files and then finally will run the code for you with a box to show the terminal output.

Running
=======

As mentioned above, if you do `python gui.py` and set up the input files, eventually you'll reach a screen with some final options and a big button that says "Run simulation". Click this button to Run simulation. The terminal output will appear in the box below the button.

By default, the Python wrapper will create a folder `out` inside the code folder, and put the output files in there with a subfolder hierarchy based on the name of the protein, the fluence and rep rate, and so on. Once you've run a simulation and it's finished, click around in that output folder, hopefully it will be fairly self-explanatory.

Alternatively you can run it manually in a terminal by calling `python main.py ` with various options. They're documented in `main.py`; you can do `python main.py -h` to see details. The compulsory first argument is `-pf FILE`, where `FILE` should be a filename containing JSON data for the protein you want to simulate.
Optional ones are:

- `-p PROTEIN_NAME`, where `PROTEIN_NAME` should be a string. The code will search for the key `PROTEIN_NAME` in the protein file given. This is so that you can (for example) make one big protein.json file with lots of details of different proteins in it and then pick one. If not given, the code will search for keys in the protein file: if there's only one, or if there's one that matches the filename, it will use that (and print a warning).
- `-o OUTPUT_PATH`, where `OUTPUT_PATH` is a string. This changes the root directory for the output files to be placed in. If not given, will default to `./out` (i.e. a new folder named "out" in your current directory).
- `-c CONNECTED`, where `CONNECTED` should be True or False. if True, hopping rates will be left as-is in order to simulate a connected aggregate of proteins; if False, this will zero out all hopping rates to simulate an ensemble of unconnected proteins (e.g. if you're experimenting on membrane proteins in detergent).
- `-n NUM_CORES`, where `N_CORES` should be a number. This is the number of cores that will be passed to the OpenMPI call in the fortran - more cores will make it run faster so long as you have the cores available on your machine. Setting this to use more cores than you actually have will make the performance worse, though. The GUI will show what `os.cpu_count()` returns, which you should probably take as a maximum recommended value, and actually the best value is quite probably smaller than that. On a modern laptop try 2 or 4 and see how it goes.

The Python wrapper will take all these options, make Fortran-friendly parameter files and copy them to the output directory, run `make all` on the Fortran if necessary, run it for you, and then perform reconvolution fits when the Fortran returns. By default it will use the decay times of the protein states as a starting point and try 1- to n-exponential fits based on those.

Note that the fitting script is not that sophisticated really; it cuts off the trace at the peak and fits 1- to n-exponentials to the tail, then fixes the fitted time constants in place and does a reconvolution fit with the IRF for the amplitudes.
It does not do anything more complicated than that; the output traces are just CSVs and therefore hopefully it should be possible to load them into another fitting program fairly easily.

Parameter JSON setup
====================

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
- `abundance` controls what fraction of sites have the given state on them. the fortran will take these fractions and randomise which specific sites have which states on them, separately for each core and each repeat
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
- `debug` should be true or false. if true, the fortran will output some extra stuff about move statistics that you probably don't need

Note: The fitting code will output various CSV and text files containing the fitted arrays, details of the fits and errors and so on, but it doesn't do anything sophisticated to compare the fit between repeats.

FAQS
====

Q: where is the test suite, the continuous integration, all that kind of stuff?  
A: what are you, a coward?
