import numpy as np

from hyperionrt.model import AnalyticalYSOModel
from hyperionrt.util.constants import msun, rsun, au, pi, sigma, c
from hyperionrt.dust import SphericalDust
from mysg.parameters import read_parfile
from mysg.atmosphere import interp_atmos


def setup_model(parfile):

    # Set dust sublimation temperature
    tsub = 1600.

    # Read in model parameters
    par = read_parfile(parfile, nested=True)

    # Set default dust if needed
    for component in ['disk', 'envelope', 'cavity', 'ambient']:
        if component in par and 'dust' not in par[component]:
            par[component] = 'kmh.hdf5'

    # Set up model
    m = AnalyticalYSOModel()

    if 'star' in par:

        # Set radius and luminosity
        m.star.radius = par['star']['radius'] * rsun
        m.star.luminosity = 4. * pi * par['star']['radius'] * rsun * sigma \
                            * par['star']['temperature'] ** 4.

        # Interpolate and set spectrum
        nu, fnu = interp_atmos(par['star']['temperature'])
        m.star.spectrum = (nu, fnu)

    else:

        raise Exception("Cannot compute a model without a central source")

    if 'disk' in par:

        # Add the flared disk component
        disk = m.add_flared_disk()

        # Accretion luminosity
        if 'lacc' in par['disk']:

            # Tell the model that we are including accretion
            m.accretion = True

            # Find the total luminosity of accretion shock on star
            lshock = par['disk']['lacc'] / 2.

            # Hot spot parameters
            fspot = 0.05
            fluxratio = 0.5 * lshock / m.star.luminosity / fspot
            tshock = par['star']['temperature'] * (1 + fluxratio) ** 0.25  # K

            # Set the hot spot source
            m.star.sources['uv'].luminosity = lshock / 2. \
                                            + m.star.luminosity * fspot
            m.star.sources['uv'].temperature = tshock

            # X-rays
            wav = np.logspace(-3., -2., 100)[::-1]
            nu = c * 1.e4 / wav
            fnu = np.repeat(1., nu.shape)

            # Set the X-ray source
            m.star.sources['xray'].luminosity = lshock / 2.
            m.star.sources['xray'].spectrum = (nu, fnu)

            # Reduce the total luminosity from the original source
            m.star.luminosity *= 1 - fspot

            # Set luminosity from viscous dissipation in disk
            disk.lacc = par['disk']['lacc'] / 2.  # incorrect

        # Basic parameters
        disk.mass = par['disk']['mass'] * msun
        disk.rmax = par['disk']['rmaxd'] * au
        disk.alpha = par['disk']['beta'] + 1
        disk.beta = par['disk']['beta']
        disk.h0 = par['disk']['h100']
        disk.r0 = 100. * au

        # Set dust
        disk.dust = par['disk']['dust']

        # Read in dust file
        d = SphericalDust(disk.dust)

        # Find the effective temperature, and spectrum of the central source
        # including accretion emission
        teff = m.star.effective_temperature()
        nu, fnu = m.star.total_spectrum()

        # Find the dust sublimation radius
        disk.rmin = m.star.radius * (1. - (1. - 2. * (tsub / teff) ** 4. \
                    * d.kappa_planck_temperature(tsub) \
                    / d.kappa_planck_spectrum(nu, fnu)) ** 2.) ** -0.5

        # Inner radius
        if 'rmin' in par['disk']:
            disk.rmin *= par['disk']['rmin']

        # Settling
        if 'eta' in par['disk']:
            raise Exception("Dust settling implemented")

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

        # Read in dust file
        d = SphericalDust(envelope.dust)

        # Find the effective temperature, and spectrum of the central source
        # including accretion emission
        teff = m.star.effective_temperature()
        nu, fnu = m.star.total_spectrum()

        # Find the dust sublimation radius
        envelope.rmin = m.star.radius * (1. - (1. - 2. * (tsub / teff) ** 4. \
                        * d.kappa_planck_temperature(tsub) \
                        / d.kappa_planck_spectrum(nu, fnu)) ** 2.) ** -0.5

        # Inner radius
        if 'rmin' in par['envelope']:
            envelope.rmin *= par['envelope']['rmin']

    if 'cavity' in par:

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
