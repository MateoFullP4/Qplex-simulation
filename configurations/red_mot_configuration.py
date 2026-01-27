"""
This file sets up the configuration for the red MOT (Magneto-Optical Trap) for Ytterbium atoms used in QPlex.

Conventions and setup:
- Atom: Ytterbium, using the "main" transition.
- Lasers: Multiple Gaussian laser beams with specified wavelength, waist, power, direction, and polarization.
  * l399_ZS_1: Zeeman slower beam
  * l399_MM_d1–d4: MOT beams in 3D configuration
- Polarization conventions: CircularLeft, CircularRight, Horizontal
- Magnetic field: Imported from 'magnetic_configuration' as 'mag_field_asymmetric'

Atom-light interactions are added with specific detunings (in units of gamma) for each beam.

This configuration can be imported in other scripts as:

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from configurations.red_mot_configuration import configuration_MM

All beams are positioned relative to the origin (0,0,0). Directions are given as (x, y, z) unit vectors.
"""

import numpy as np
from atomsmltr.atoms import Ytterbium
from atomsmltr.environment.lasers import GaussianLaserBeam
from atomsmltr.environment.zones import Limits, Box
from atomsmltr.simulation import Configuration
from atomsmltr.environment.lasers.polarization import CircularLeft, CircularRight, Horizontal
import os
import sys


# --- Import Magnetic Configuration ---
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from configurations.magnetic_configuration import mag_field_asymmetric

# --- Atom Setup ---
atom = Ytterbium()
main = atom.trans["main"]
gamma = main.Gamma

# --- Laser Setup ---
l399_ZS_1 = GaussianLaserBeam(
    wavelength = 399e-9,
    waist = 6.7e-3, 
    power = 225e-3,
    waist_position = (0, 0, 0),
    direction = (0, 0, -1), 
    polarization = Horizontal(),
    tag = "laser_ZS_1"
)

l399_MM_d1 = GaussianLaserBeam(
    wavelength = 399e-9,
    waist = 25e-3,
    power = 100e-3,
    waist_position = (0, 0, 0), 
    direction = (1, 1, 1),
    polarization = CircularRight(),
    tag = "laser_MM_d1"
)

l399_MM_d2 = GaussianLaserBeam(
    wavelength = 399e-9,
    waist = 25e-3,
    power = 100e-3,
    waist_position = (0, 0, 0), 
    direction = (-1, 1, -1),
    polarization = CircularRight(),
    tag = "laser_MM_d2"
)

l399_MM_d3 = GaussianLaserBeam(
    wavelength = 399e-9,
    waist = 25e-3,
    power = 100e-3,
    waist_position = (0, 0, 0), 
    direction = (1, -1, -1),
    polarization = CircularLeft(),
    tag = "laser_MM_d3"
)

l399_MM_d4 = GaussianLaserBeam(
    wavelength = 399e-9,
    waist = 25e-3,
    power = 100e-3,
    waist_position = (0, 0, 0), 
    direction = (-1, -1, 1),
    polarization = CircularLeft(),
    tag = "laser_MM_d4"
)

# --- Configuration Setup ---
configuration_MM = Configuration()
configuration_MM.atom = atom
configuration_MM += l399_MM_d1, l399_MM_d2, l399_MM_d3, l399_MM_d4, l399_ZS_1, mag_field_asymmetric

# --- Add Atom-light Coupling ---
configuration_MM.add_atomlight_coupling("laser_MM_d1", "main", -1.0*gamma)
configuration_MM.add_atomlight_coupling("laser_MM_d2", "main", -1.0*gamma)
configuration_MM.add_atomlight_coupling("laser_MM_d3", "main", -1.4*gamma)
configuration_MM.add_atomlight_coupling("laser_MM_d4", "main", -1.4*gamma)
configuration_MM.add_atomlight_coupling("laser_ZS_1", "main", -16.5*gamma)


