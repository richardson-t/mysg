import numpy as np

from hyperionrt.model import AnalyticalYSOModel
from hyperionrt.util.constants import msun, rsun, au, pi, sigma, c
from hyperionrt.dust import SphericalDust
from mysg.parameters import read_parfile
from mysg.atmosphere import interp_atmos


def set_up_model(parfile):

    # Set dust sublimation temperature
    tsub = 1600.

    # Read in model parameters
    par = read_parfile(parfile, nested=True)

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

        # Accretion
        if 'lacc' in par['disk']:

            # Total luminosity of accretion shock on star
            lshock = par['disk']['lacc'] / 2.

            # Hot spot parameters
            fspot = 0.05
            fluxratio = 0.5 * lshock / m.star.luminosity / fspot
            tshock = par['star']['temperature'] * (1 + fluxratio) ** 0.25  # K
            m.star.luminosity *= 1 - fspot
            m.star.sources['uv'].luminosity = lshock / 2. \
                                            + m.star.luminosity * fspot
            m.star.sources['uv'].temperature = tshock

            # X-rays
            wav = np.logspace(-3., -2., 100)[::-1]
            nu = c * 1.e4 / wav
            fnu = np.repeat(1., nu.shape)
            m.star.sources['xray'].luminosity = lshock / 2.
            m.star.sources['xray'].spectrum = (nu, fnu)

        # Add the flared disk component
        disk = m.add_flared_disk()

        # Basic parameters
        disk.mass = par['disk']['mass'] * msun
        disk.rmax = par['disk']['rmaxd'] * au
        disk.alpha = par['disk']['beta'] + 1
        disk.beta = par['disk']['beta']
        disk.h0 = par['disk']['h100']
        disk.r0 = 100. * au

        # Read in dust file
        d = SphericalDust(par['disk']['dust'])

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

        # Read in dust file
        d = SphericalDust(par['envelope']['dust'])

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

        pass

    if 'ambient' in par:

        pass
