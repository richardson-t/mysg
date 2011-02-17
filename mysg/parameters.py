import numpy as np

from mysg.odict import odict


def write_parfile(filename, parameters):
    f = open(filename, 'wb')
    for name in parameters:
        if type(parameters[name]) not in [str, np.string_]:
            f.write("%s = %9.3e\n" % (name, parameters[name]))
        else:
            f.write("%s = %s\n" % (name, parameters[name]))
    f.close()


def read_parfile(filename, nested=False):
    parameters = odict()
    f = open(filename, 'rb')
    for line in f.readlines():
        key, value = line.strip().split('=')
        try:
            key, value = key.strip(), float(value)
        except ValueError:
            key, value = key.strip(), value.strip()
        if nested:
            p = parameters
            components = key.split('.')
            for component in components[:-1]:
                if component not in p:
                    p[component] = odict()
                p = p[component]
            p[components[-1]] = value
        else:
            parameters[key] = value
    f.close()
    return parameters
