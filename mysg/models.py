import random

import numpy as np
import atpy

from mysg.odict import odict
from mysg.ranges import read_ranges, write_ranges
from mysg.util import create_dir, random_id
from mysg.parameters import required_parameters


def sample_set_models(set_name, number, seed=123456789):

    # Ensure reproducibility
    random.seed(seed)
    np.random.seed(seed=seed)

    # Read in the ranges
    ranges = read_ranges("models/%s/ranges.conf" % set_name)

    # First pass, sample all values that need to be sampled
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
        elif par['sampling'] == 'fixed':
            values.add_column(name, np.repeat(par['value'], number))
        elif par['sampling'] == 'linked':
            pass
        else:
            raise Exception("Unknown sampling: %s" % par['sampling'])

    # Second pass, deal with linked parameters
    for name in ranges:
        if ranges[name]['sampling'] == 'linked':
            values.add_column(name, values[ranges[name]['parameter']])

    # Write out parameter files
    create_dir("models/%s/par" % set_name)

    for i in range(len(values)):
        model_name = random_id()
        f = open("models/%s/par/%s.par" % (set_name, model_name), 'wb')
        for name in values.columns:
            f.write("%s = %9.3e\n" % (name, values[name][i]))
        f.close()

    # Write out table
    values.write("models/%s/parameters.hdf5" % set_name)


def make_set_dir(set_name, ranges_file):
    '''
    Given the name of a model set and a file containing ranges, set up the
    directory with a subset of the ranges file
    '''

    # Find out which parameters are required
    required = required_parameters(set_name)

    # Read in ranges
    ranges = read_ranges(ranges_file)

    # Define new subset of ranges
    ranges_new = odict()
    for parameter in required:
        if ranges[parameter]['sampling'] == 'linked' \
           and ranges[parameter]['parameter'] not in required:
            ranges_new[parameter] = ranges[ranges[parameter]['parameter']]
        else:
            ranges_new[parameter] = ranges[parameter]

    # Create directory
    create_dir('models/%s' % set_name)

    # Write out ranges file
    write_ranges('models/%s/ranges.conf' % set_name, ranges_new)
