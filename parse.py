# callum - 14/09/26
# -*- coding: utf-8 -*-

import os
import json
import numpy as np
'''
parsers for the protein and simulation JSON files, as well as
functions to generate the output directory hierarchies and then
generate and place the fortran-friendly input files within those
directories. moved here from main.py because it was all very
haphazard before.

for protein_spec and simulation_spec the logic is to define the name
and then a set of specifications like the type, lower bound, upper bound etc.
the reason for doing it this way is that i realised i was writing the same
sets of error messages over and over for different parameters and it was
pointless. also i don't do any checking of the pigment and state names so
don't just paste in random strings you found on the internet or anything, i
will not be held responsible for that
'''

protein_spec = {
    "n_p":           {'ptype': int,   'lb': 1},
    "n_s":           {'ptype': int,   'lb': 1},
    "pigment_names": {'ptype': str,   'size': 'n_p'},
    "state_names":   {'ptype': str,   'size': 'n_s'},
    "which_pigment": {'ptype': int,   'lb': 1, 'size': 'n_s'},
    "xsec":          {'ptype': float, 'lb': 0.0, 'size': 'n_s'},
    "emissive":      {'ptype': bool,  'size': 'n_s'},
    "abundance":     {'ptype': float, 'lb': 0.0, 'ub': 1.0, 'size': 'n_s'},
    "dist":          {'ptype': bool,  'size': ['n_s', 'n_s']},
    "n_tot":         {'ptype': int,   'size': 'n_p'},
    "n_thermal":     {'ptype': int,   'size': 'n_p'},
    "hop":           {'ptype': float, 'lb': 0.0, 'size': 'n_s'},
    "intra":         {'ptype': float, 'lb': 0.0, 'size': ['n_s', 'n_s']},
    "ann":           {'ptype': float, 'lb': 0.0, 'size': ['n_s', 'n_s'],
        'symmetry': True},
    "ann_remainder": {'ptype': int, 'lb': 0, 'size': ['n_s', 'n_s'],
        'symmetry': True},
}

simulation_spec = {
  "fwhm":      {'ptype': float, 'lb': 0.0},
  "fluence":   {'ptype': float, 'lb': 0.0},
  "n_sites":   {'ptype': int,   'lb': 1},
  "lattice":   {'ptype': str,
      'choices': ["hex", "honeycomb", "square", "line"]},
  "rep_rate":  {'ptype': float, 'lb': 0.0},
  "burn_reps": {'ptype': int,   'lb': 0},
  "tmax":      {'ptype': float, 'lb': 0.0},
  "dt1":       {'ptype': float, 'lb': 0.0},
  "dt2":       {'ptype': float, 'lb': 0.0},
  "binwidth":  {'ptype': float, 'lb': 0.0},
  "n_counts":  {'ptype': int,   'lb': 1},
  "n_repeats": {'ptype': int,   'lb': 1},
  "debug":     {'ptype': bool},
}

def check(data, key, ptype=None, lb=None, ub=None,
        size=None, symmetry=False, choices=None):
    '''
    check that the requirements for a given parameter which are specified
    above in protein_spec and simulation_spec are met. return False and
    a set of descriptive error messages if they are not.
    note that i don't do any checking on ptype; as long as it's only
    used with the gui code i've written this shouldn't matter but if you
    start messing around with the JSON and putting weird stuff in there
    it'll probably break, that's on you. also note that the bounds
    lb and ub are exclusive, i.e. if you want > 0, set lb = 1
    '''
    valid = True
    msgs = []
    # make the data iterable even if it's a single value, so that
    # the list comp for type checking will always work
    arr = np.array(data[key]).astype(ptype)
    if any([isinstance(ptype, type(item)) for item in arr.flatten()]):
        valid = False
        msgs.append(f"Non-{ptype} value given in {key}: "
        f"{data[key]} has type {type(arr.flatten()[0])}")
    if lb is not None:
        if np.any(arr < lb):
            valid = False
            msgs.append(f"{key} = {data[key]} has values < {lb}.")
    if ub is not None:
        if np.any(arr > ub):
            valid = False
            msgs.append(f"{key} = {data[key]} has values > {ub}.")
    if size is not None:
        st = size if type(size) is list else [size]
        shape = tuple([data[dim] for dim in st])
        if arr.shape != shape:
            valid = False
            msgs.append(f"{key} = {data[key]} has invalid shape {shape}, {arr.shape}.")
    if choices is not None: # only for lattice atm so we can just check one
        if data[key] not in choices:
            valid = False
            msgs.append(f"{key} = {data[key]} has items not in {choices}")
    if symmetry:
        if not np.allclose(arr, arr.T):
            valid = False
            msgs.append(f"{key} = {data[key]} should be symmetric.")
    return valid, msgs

def parse(data, spec, keys=None, necessary_keys=None):
    '''
    check an dictionary `data` against a specification `spec`, as
    given above in protein_spec and simulation_spec, and make sure that:
    - it does not contain any keys that shouldn't be there
    - all the keys that are there have corresponding values that
      are valid wrt the rest of the code.
    if a list of keys is given, check only those keys (this allows us
    to use the same function for each wizard page); otherwise check them all.
    '''
    validated = True
    msgs = []
    allowed_keys = spec.keys()
    if keys is None:
        # if no keys are given, assume we're checking everything;
        # in which case, all keys should be present
        keys = data.keys()
        if allowed_keys != keys:
            validated = False
            # check which keys must be changed in data for them to match
            ktr = list(keys - allowed_keys)
            kta = list(allowed_keys - keys)
            msg = ["Data's keys don't match allowed keys."]
            if len(ktr) > 0:
                msg.append(f"Remove keys {ktr} from data.")
            if len(kta) > 0:
                msg.append(f"Add keys {kta} to data.")
            msgs.append((' ').join(msg))
            return validated, msgs
    '''
    keys should now be either the list passed to the function, or
    the set of all allowed keys (if not, the function will have returned).
    however, if a list was passed, we should do a quick sanity check: 
    any necessary_keys, if given, should be in the data regardless.
    '''
    if necessary_keys is not None:
        for key in necessary_keys:
            if key not in data.keys():
                validated = False
                msgs.append(f"Key '{key}' is missing from data.")
                return validated, msgs

    # the list of keys passed should be a subset of the allowed keys
    if not set(keys) <= allowed_keys:
        validated = False
        msgs.append("List of keys passed to parse.parse() is "
                "not a subset of the allowed keys. Invalid keys: "
                f"{keys - allowed_keys}.")
        return validated, msgs

    # now more detailed conditions for each key
    for key in keys:
        reqs = spec[key]
        valid, curr_msgs = check(data, key, **reqs)
        if not valid:
            validated = False
        # if validated is True and curr_msgs is [], this is a no-op
        msgs.extend(curr_msgs)
    return validated, msgs

def parse_protein(data, keys=None):
    '''
    convenience wrapper for parsing protein data specifically.
    the reason for doing it this way is that the protein requires the
    keys 'n_p' and 'n_s' because other parameters depend on those,
    whereas none of the simulation parameters strictly depend on any
    of the others. this way i can write one parse function and just
    call it with different args rather than duplicating code.
    '''
    return parse(data, protein_spec, keys, ['n_p', 'n_s'])

def parse_simulation(data, keys=None):
    '''
    convenience wrapper for parsing simulation data specifically.
    see above.
    '''
    return parse(data, simulation_spec, keys)

