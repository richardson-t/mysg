import numpy as np

from hyperionrt.model import AnalyticalYSOModel
from hyperionrt.util.constants import msun, rsun, au, pi, sigma, lsun
from hyperionrt.util.convenience import OptThinRadius

from mysg.parameters import read_parfile
from mysg.atmosphere import interp_atmos

# Set dust sublimation temperature
tsub = 1600.


def setup_model(parfile):

    # Read in model parameters
    par = read_parfile(parfile, nested=True)

    # Set default dust if needed
    for component in ['disk', 'envelope', 'cavity', 'ambient']:
        if component in par and 'dust' not in par[component]:
            par[component] = 'kmh.hdf5'

    # Set up model
    m = AnalyticalYSOModel()

    if not 'star' in par:
        raise Exception("Cannot compute a model without a central source")

    # Set radius and luminosity
    m.star.radius = par['star']['radius'] * rsun
    m.star.luminosity = 4. * pi * par['star']['radius'] * rsun * sigma \
                        * par['star']['temperature'] ** 4.

    # Interpolate and set spectrum
    nu, fnu = interp_atmos(par['star']['temperature'])
    m.star.spectrum = (nu, fnu)

    if 'disk' in par:

        # Add the flared disk component
        disk = m.add_flared_disk()

        # Basic parameters
        disk.mass = par['disk']['mass'] * msun
        disk.rmax = par['disk']['rmaxd'] * au
        disk.alpha = par['disk']['beta'] + 1
        disk.beta = par['disk']['beta']
        disk.h0 = par['disk']['h100']
        disk.r0 = 100. * au

        # Set dust
        disk.dust = par['disk']['dust']

        # Inner radius
        if 'rmin' in par['disk']:
            disk.rmin = par['disk']['rmin'] * OptThinRadius(tsub)
        else:
            disk.rmin = OptThinRadius(tsub)

        # Settling
        if 'eta' in par['disk']:
            raise Exception("Dust settling implemented")

        # Accretion luminosity
        if 'lacc' in par['disk']:
            m.add_magnetospheric_accretion(par['disk']['lacc'] * lsun,
                                           par['rtrunc'], par['fspot'], disk)

    if 'envelope' in par:

        if 'rc' in par['envelope']:  # Ulrich envelope

            envelope = m.add_ulrich_envelope()
            envelope.rho0 = par['envelope']['rho0']
            envelope.rc = par['envelope']['rc'] * au

        elif 'power' in par['envelope']:  # Power-law envelope

            envelope = m.add_powerlaw_envelope()
            envelope.power = par['envelope']['power']
            envelope.mass = par['envelope']['mass'] * msun

        # Set dust
        envelope.dust = par['envelope']['dust']

        # Inner radius
        if 'rmin' in par['envelope']:
            envelope.rmin = par['envelope']['rmin'] * OptThinRadius(tsub)
        else:
            envelope.rmin = OptThinRadius(tsub)

    if 'cavity' in par:

        if not 'envelope' in par:
            raise Exception("Can't have a bipolar cavity without an envelope")

        # Add the bipolar cavity component
        cavity = envelope.add_bipolar_cavity()

        # Basic parameters
        cavity.power = par['cavity']['power']
        cavity.r_0 = 10000 * au
        cavity.theta_0 = par['cavity']['theta_0']
        cavity.rho_0 = par['cavity']['rho_0']
        cavity.rho_exp = par['cavity']['rho_exp']

        # Set dust
        cavity.dust = par['cavity']['dust']

    if 'ambient' in par:

        # This is where we set the envelope outer radius
        pass

    # Set up run-time parameters
    m.set_raytracing(True)
    m.set_mrw(True, gamma=2.)

    # Set up SEDs
    image = m.add_peeled_images()
    image.set_wavelength_range(250, 0.001, 5000.)
    image.set_image_size(1, 1)
    image.set_image_limits(-np.inf, np.inf, -np.inf, np.inf)
    image.set_aperture_range(1, np.inf, np.inf)  # needs changing
    image.set_output_bytes(4)
