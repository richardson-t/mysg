from mysg.odict import odict


def _fixed(value):
    return {'sampling': 'fixed', 'value': value}


def write_ranges(filename, ranges):
    "Write a file containing information about ranges for each parameter"
    f = open(filename, 'wb')
    for name in ranges:
        par = ranges[name]
        if par['sampling'] in ['log10', 'linear']:
            f.write("%-25s %-10s %10.3e  %10.3e\n"
                    % (name, par['sampling'], par['lower'], par['upper']))
        elif par['sampling'] in ['fixed']:
            f.write("%-25s %-10s %10.3e\n"
                    % (name, par['sampling'], par['value']))
        elif par['sampling'] in ['str']:
            f.write("%-25s %-10s  %s\n"
                    % (name, par['sampling'], par['value']))
        elif par['sampling'] in ['linked']:
            f.write("%-25s %-10s  %s\n"
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
        elif sampling in ['fixed']:
            ranges[name]['value'] = float(details)
        elif sampling in ['str']:
            ranges[name]['value'] = details
        elif sampling in ['linked']:
            ranges[name]['parameter'] = details
        else:
            raise Exception("Unknown sampling: %s" % sampling)
    f.close()
    return ranges


def select_required_ranges(set_name, ranges):

    # Add parameters
    required = []
    if set_name[0] == 's':
        required += ['star.radius', 'star.temperature']
    if set_name[1] == 'd':
        required += ['disk.rmax', 'disk.mass', 'disk.beta', 'disk.h100']
    if set_name[2] == 'u':
        required += ['envelope.rc', 'envelope.rho0']
    if set_name[2] == 'p':
        required += ['envelope.rmin', 'envelope.power', 'envelope.mass']
    if set_name[3] == 'c':
        required += ['cavity.theta_0', 'cavity.rho_0', 'cavity.rho_exp', 'cavity.power']
    if set_name[4] == 'h':
        if set_name[1] == 'd':
            required += ['disk.rmin']
        if set_name[1] == 'e':
            required += ['envelope.rmin']
    if set_name[6] == 'a':
        required += ['disk.rtrunc', 'disk.lacc', 'star.fspot']
    if set_name[7] == 'g':
        required += ['disk.eta']

    for parameter in required:
        if parameter not in ranges:
            raise Exception("Parameter %s missing from ranges" % parameter)

    ranges_new = odict()
    for parameter in required:
        if ranges[parameter]['sampling'] == 'linked' \
           and ranges[parameter]['parameter'] not in required:
            ranges_new[parameter] = ranges[ranges[parameter]['parameter']]
        else:
            ranges_new[parameter] = ranges[parameter]

    # Sets parameters that depend on flags
    if set_name[5] == 'm':
        ranges_new['ambient.density'] = _fixed(1.e-22)
        ranges_new['ambient.temperature'] = _fixed(10.)
    if set_name[8] == 'g':
        ranges_new['disk.dust.amax'] = _fixed(1000.)
    if set_name[9] == 'p':
        ranges_new['disk.dust.pah_frac'] = _fixed(0.18)
        ranges_new['envelope.dust.pah_frac'] = _fixed(0.18)
        ranges_new['cavity.dust.pah_frac'] = _fixed(0.18)
        ranges_new['ambient.dust.pah_frac'] = _fixed(0.18)

    return ranges_new

# In future, can define function to find intersection of two ranges or
# compare the parameter space volume of two set of ranges.
