import glob
import os
import shutil
import random

import numpy as np
import atpy

from mysg.ranges import read_ranges, write_ranges, select_required_ranges
from mysg.util import create_dir, random_id
from mysg.parameters import write_parfile
from mysg.odict import odict

VALID = []
VALID.append(['-', 's'])
VALID.append(['-', 'p', 'a'])
VALID.append(['-', 'p', 'u'])
VALID.append(['-', 'b'])
VALID.append(['s', 'h'])
VALID.append(['-', 'm'])
VALID.append(['-', 'i'])
VALID.append(['-', 'p'])
VALID.append(['-'])
VALID.append(['-'])
VALID.append(['-'])
VALID.append(['-'])

# Set path with dust files
DATA_PATH = os.path.join(os.path.dirname(__file__), 'data/dust/')


def export_dust_files():

    # Copy all dust files there
    for filename in glob.glob(os.path.join(DATA_PATH, '*.hdf5')):
        shutil.copy(filename, os.path.basename(filename))


def _check_set_name(set_name):
    for i, letter in enumerate(set_name):
        if letter not in VALID[i]:
            raise Exception("Letter %i cannot be %s" % (i + 1, letter))


def sample_set_models(directory, set_name, number_function):

    # Ensure reproducibility
    random.seed(abs(hash(set_name)))
    np.random.seed(seed=abs(hash(set_name)))

    # Read in the ranges
    ranges = read_ranges("%s/%s/ranges.conf" % (directory, set_name))

    # Find total number of free parameters
    n_free = 0
    for name in ranges:
        par = ranges[name]
        if par['sampling'] in ['linear', 'log10']:
            n_free += 1

    number = number_function(n_free)

    # Sample all values that need to be sampled
    values = atpy.Table()
    for name in ranges:
        par = ranges[name]
        if par['sampling'] == 'linear':
            values.add_column(name, np.random.uniform(par['lower'],
                                                      par['upper'],
                                                      number))
        elif par['sampling'] == 'log10':
            values.add_column(name,
                              10. ** np.random.uniform(np.log10(par['lower']),
                                                       np.log10(par['upper']),
                                                       number))
        elif par['sampling'] in ['fixed', 'str']:
            values.add_column(name, np.repeat(par['value'], number))
        elif par['sampling'] == 'linked':
            values.add_column(name, values[ranges[name]['parameter']])
        else:
            raise Exception("Unknown sampling: %s" % par['sampling'])

    # Write out parameter files
    create_dir("%s/%s/par" % (directory, set_name))
    create_dir("%s/%s/log" % (directory, set_name))
    create_dir("%s/%s/input" % (directory, set_name))
    create_dir("%s/%s/output" % (directory, set_name))

    for i in range(len(values)):
        model_name = random_id()
        write_parfile("%s/%s/par/%s.par" % (directory, set_name, model_name), odict(zip(values.keys(), values[i])))

    # Write out table
    values.write("%s/%s/parameters.hdf5" % (directory, set_name), verbose=False)


def make_set_dir(directory, set_name):
    '''
    Given the name of a model set and a file containing ranges, set up the
    directory with a subset of the ranges file
    '''

    # Check set name
    _check_set_name(set_name)

    # Select ranges that are actually needed
    ranges = select_required_ranges(set_name)

    # Create directory
    create_dir('%s/%s' % (directory, set_name))

    # Write out ranges file
    write_ranges('%s/%s/ranges.conf' % (directory, set_name), ranges)
