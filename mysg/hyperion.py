from __future__ import absolute_import

import numpy as np

from hyperion.model import AnalyticalYSOModel
from hyperion.util.constants import msun, rsun, au, pi, sigma
from hyperion.util.convenience import OptThinRadius
from hyperion.dust import SphericalDust

from mysg.parameters import read_parfile
from mysg.atmosphere import interp_atmos

# Set dust sublimation temperature
tsub = 1600.


def setup_model(parfile, output):

    # Read in model parameters
    par = read_parfile(parfile, nested=True)

    # Find dimensionality of problem:
    if 'disk' in par:
        ndim = 2
    elif 'cavity' in par:
        ndim = 2
    elif 'envelope' in par and 'rc' in par['envelope']:
        ndim = 2
    else:
        ndim = 1

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
            raise Exception("Dust settling not implemented")

        # Accretion luminosity
        if 'lacc' in par['disk']:
            raise Exception("Accretion luminosity not implemented")
            # m.setup_magnetospheric_accretion(par['disk']['lacc'] * lsun,
            #                                  par['disk']['rtrunc'],
            #                                  par['star']['fspot'], disk)

    if 'envelope' in par:

        if 'rc' in par['envelope']:  # Ulrich envelope

            envelope = m.add_ulrich_envelope()
            envelope.rho_0 = par['envelope']['rho_0']
            envelope.rc = par['envelope']['rc'] * au

        elif 'power' in par['envelope']:  # Power-law envelope

            envelope = m.add_power_law_envelope()
            envelope.power = par['envelope']['power']
            envelope.rho_0 = par['envelope']['rho_0']
            envelope.r_0 = 100. * au

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
        cavity.rho_exp = 0.

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

            # Find radius where the optically thin temperature drops to the
            # ambient temperature. We can do this only if we've already set
            # up all the sources of emission beforehand (which we have)
            rmax_temp = OptThinRadius(ambient.temperature).evaluate(m.star, envelope.dust)

            # Find radius where the envelope density drops to the ambient density
            rmax_dens = envelope.outermost_radius(ambient.rho)

            # Pick the largest
            if rmax_temp < rmax_dens:
                print "Setting envelope outer radius to that where rho(r) = rho_amb"
                envelope.rmax = rmax_dens
            else:
                print "Setting envelope outer radius to that where T_thin(r) = T_amb"
                envelope.rmax = OptThinRadius(ambient.temperature)

            ambient.rmax = envelope.rmax

        else:

            ambient.rmax = OptThinRadius(ambient.temperature)

        # The inner radius for the ambient medium should be the largest of
        # the inner radii for the disk and envelope
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

        # The ambient medium needs to go out to sqrt(2.) times the envelope
        # radius to make sure the slab is full (don't need to do sqrt(3)
        # because we only need a cylinder along line of sight)
        ambient.rmax *= np.sqrt(2.)

        # Make sure that the temperature in the model is always at least
        # the ambient temperature
        m.set_minimum_temperature(ambient.temperature)

    else:

        # Make sure that the temperature in the model is always at least
        # the CMB temperature
        m.set_minimum_temperature(2.725)

        if 'envelope' in par:
            raise Exception("Can't have an envelope without an ambient medium")

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

    # Set up grid.
    if ndim == 1:
        m.set_spherical_polar_grid_auto(400, 1, 1)
    else:
        m.set_spherical_polar_grid_auto(400, 300, 1)

    # Find the range of radii spanned by the grid
    rmin, rmax = m.radial_range()

    # Set up SEDs
    image = m.add_peeled_images(sed=True, image=False)
    image.set_wavelength_range(250, 1.e-2, 1.e+4)

    if 'ambient' in par:
        image.set_aperture_range(10, rmin, rmax / np.sqrt(2.))
    else:
        image.set_aperture_range(10, rmin, rmax)

    image.set_output_bytes(8)
    image.set_track_origin(True)
    image.set_uncertainties(True)

    if ndim == 1:
        image.set_viewing_angles([45.], [45.])
    else:
        image.set_viewing_angles(np.linspace(0., 90., 10), np.repeat(45., 10))

    if 'ambient' in par:  # take a slab to avoid spherical geometrical effects
        w = ambient.rmax / np.sqrt(2.)
        image.set_depth(-w, w)
    else:  # don't need to take a slab, as no ambient material or envelope
        image.set_depth(-np.inf, np.inf)

    # Set number of photons
    if ndim == 1:
        m.set_n_photons(temperature=1000, imaging=1000000,
                        raytracing_sources=1000000, raytracing_dust=1000000,
                        stats=100)
    else:
        m.set_n_photons(temperature=1000000, imaging=1000000,
                        raytracing_sources=1000000, raytracing_dust=1000000,
                        stats=10000)

    print 'Dimensionality: ', ndim

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
