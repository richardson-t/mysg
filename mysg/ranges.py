from mysg.odict import odict


def write_ranges(filename, ranges):
    "Write a file containing information about ranges for each parameter"
    f = open(filename, 'wb')
    for name in ranges:
        par = ranges[name]
        if par['sampling'] in ['log10', 'linear']:
            f.write("%-20s %-10s %9.3e  %9.3e\n"
                    % (name, par['sampling'], par['lower'], par['upper']))
        elif par['sampling'] in ['linked']:
            f.write("%-20s %-10s %s\n"
                    % (name, par['sampling'], par['parameter']))
        else:
            raise Exception("Unknown sampling: %s" % par['sampling'])
    f.close()


def read_ranges(filename):
    "Read a file containing information about ranges for each parameter"
    ranges = odict()
    f = open(filename, 'rb')
    for line in f.readlines():
        name, sampling, details = line.strip().split(None, 2)
        ranges[name] = odict()
        ranges[name]['sampling'] = sampling
        if sampling in ['log10', 'linear']:
            lower, upper = details.split()
            ranges[name]['lower'] = float(lower)
            ranges[name]['upper'] = float(upper)
        elif sampling == 'linked':
            ranges[name]['parameter'] = details
        else:
            raise Exception("Unknown sampling: %s" % sampling)
    f.close()
    return ranges

# In future, can define function to find intersection of two ranges or
# compare the parameter space volume of two set of ranges.
