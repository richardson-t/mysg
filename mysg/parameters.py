def required_parameters(model):
    required = []
    if model[0] == 's':
        required += ['star.radius', 'star.temperature']
    if model[1] == 'd':
        required += ['disk.rmin', 'disk.rmax', 'disk.mass', 'disk.beta',
                     'disk.h100']
    if model[2] == 'e':
        required += ['envelope.rmin', 'envelope.rc', 'envelope.rho0']
    if model[3] == 'c':
        required += ['cav.theta0', 'cav.rho0']
    if model[4] == 'h':
        if model[1] == 'd':
            required += ['disk.rmin']
        if model[1] == 'e':
            required += ['envelope.rmin']
    if model[5] == 'a':
        required += ['disk.lacc']
    required = list(set(required))
    required.sort()
    return required
