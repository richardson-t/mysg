from __future__ import absolute_import

import numpy as np

from hyperion.model import AnalyticalYSOModel
from hyperion.util.constants import msun, rsun, au, pi, sigma, lsun
from hyperion.util.convenience import OptThinRadius
from hyperion.dust import SphericalDust

from mysg.parameters import read_parfile
from mysg.atmosphere import interp_atmos

# Set dust sublimation temperature
tsub = 1600.


def setup_model(parfile, output):

    # Read in model parameters
    par = read_parfile(parfile, nested=True)

    # Set default dust if needed
    for component in ['disk', 'envelope', 'cavity', 'ambient']:
        if component in par and 'dust' not in par[component]:
            par[component]['dust'] = 'kmh.hdf5'

    # Set up model
    m = AnalyticalYSOModel(output)

    if not 'star' in par:
        raise Exception("Cannot compute a model without a central source")

    # Set radius and luminosity
    m.star.radius = par['star']['radius'] * rsun
    m.star.luminosity = 4. * pi * (par['star']['radius'] * rsun) ** 2. \
                        * sigma * par['star']['temperature'] ** 4.

    # Interpolate and set spectrum
    nu, fnu = interp_atmos(par['star']['temperature'])
    m.star.spectrum = (nu, fnu)

    if 'disk' in par:

        # Add the flared disk component
        disk = m.add_flared_disk()

        # Basic parameters
        disk.mass = par['disk']['mass'] * msun
        disk.rmax = par['disk']['rmax'] * au
        disk.p = par['disk']['p']
        disk.beta = par['disk']['beta']
        disk.h_0 = par['disk']['h100'] * au
        disk.r_0 = 100. * au

        # Set dust
        disk.dust = SphericalDust(par['disk']['dust'])

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
            m.setup_magnetospheric_accretion(par['disk']['lacc'] * lsun,
                                             par['disk']['rtrunc'],
                                             par['star']['fspot'], disk)

    if 'envelope' in par:

        if 'rc' in par['envelope']:  # Ulrich envelope

            envelope = m.add_ulrich_envelope()
            envelope.rho_0 = par['envelope']['rho_0']
            envelope.rc = par['envelope']['rc'] * au

        elif 'power' in par['envelope']:  # Power-law envelope

            envelope = m.add_power_law_envelope()
            envelope.power = par['envelope']['power']
            envelope.rho_0 = par['envelope']['rho_0']
            envelope.r_0 = par['envelope']['r_0']

        # Set dust
        envelope.dust = SphericalDust(par['envelope']['dust'])

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
        cavity.dust = SphericalDust(par['cavity']['dust'])

    if 'ambient' in par:

        # Add the ambient medium contribution
        ambient = m.add_ambient_medium()

        # Set the density, temperature, and dust properties
        ambient.rho = par['ambient']['density']
        ambient.temperature = par['ambient']['temperature']
        ambient.dust = SphericalDust(par['ambient']['dust'])

        # If there is an envelope, set the outer radius to where the
        # optically thin temperature would transition to the ambient medium
        # temperature
        if 'envelope' in par:
            envelope.rmax = OptThinRadius(ambient.temperature)

        # The inner radius for the ambient medium should be the largest of
        # the inner radii for the disk and envelope
        if 'disk' in par and 'rmin' in par['disk']:
            ambient.rmin = par['disk']['rmin'] * OptThinRadius(tsub)
        if 'envelope' in par and 'rmin' in par['envelope']:
            if 'disk' in par and 'rmin' in par['disk']:
                ambient.rmin = max(par['disk']['rmin'], \
                                   par['envelope']['rmin']) \
                               * OptThinRadius(tsub)
            else:
                ambient.rmin = par['envelope']['rmin'] * OptThinRadius(tsub)
        elif 'disk' in par and 'rmin' in par['disk']:
            ambient.rmin = par['disk']['rmin'] * OptThinRadius(tsub)
        else:
            ambient.rmin = OptThinRadius(tsub)

        # The ambient medium needs to go out to sqrt(3.) times the envelope
        # radius to make sure the slab is full
        ambient.rmax = OptThinRadius(ambient.temperature) * np.sqrt(3.)

        # Make sure that the temperature in the model is always at least
        # the ambient temperature
        m.set_minimum_temperature(ambient.temperature)

    else:

        if 'envelope' in par:
            envelope.rmax = OptThinRadius(2.725)

    # Use raytracing to improve s/n of thermal/source emission
    m.set_raytracing(True)

    # Use the modified random walk
    m.set_mrw(True, gamma=2.)

    # Use the partial diffusion approximation
    m.set_pda(True)

    # Improve s/n of scattering by forcing the first interaction
    m.set_forced_first_scattering(True)

    # Use slow dust sublimation
    m.set_dust_sublimation('slow')

    # Set physical array output to 32-bit
    m.set_output_bytes(4)

    # Set up grid
    m.set_spherical_polar_grid_auto(399, 199, 1)

    # Find the range of radii spanned by the grid
    rmin, rmax = m.radial_range()

    # Set up SEDs
    image = m.add_peeled_images(sed=True, image=False)
    image.set_wavelength_range(250, 1.e-2, 1.e+4)
    image.set_aperture_range(10, rmin, rmax)
    image.set_output_bytes(4)
    image.set_track_origin(True)
    image.set_uncertainties(True)
    image.set_viewing_angles(np.linspace(0., 90., 10), np.repeat(45., 10))
    image.set_depth(-np.inf, np.inf)

    # Set number of photons
    m.set_n_photons(temperature=1000000, imaging=1000000,
                    raytracing_sources=1000000, raytracing_dust=1000000,
                    stats=10000)

    # Request 32-bit output
    m.set_output_bytes(4)

    # Only request certain arrays to be output
    m.conf.output.output_temperature = 'none'
    m.conf.output.output_density = 'none'
    m.conf.output.output_specific_energy_abs = 'last'
    m.conf.output.output_n_photons = 'none'
    m.conf.output.output_density_diff = 'last'

    # Set number of temperature iterations and convergence criterion
    m.set_n_temperature_iterations(10)
    m.set_convergence(True, percentile=99.0, absolute=2.0, relative=1.1)

    # Write out file
    m.write(copy_dust=False, absolute_paths=False,
            physics_dtype=np.float32, wall_dtype=float)
