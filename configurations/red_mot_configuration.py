"""
This file sets up the configuration for the blue MOT (Magneto-Optical Trap) for Ytterbium atoms used in QPlex.

Conventions and setup:
- Atom: Ytterbium, using the "main" transition.
- Lasers: Multiple Gaussian laser beams with specified wavelength, waist, power, direction, and polarization.
  * l399_NB_1-4: Narrow Band beams in 2D MOT configuration 
  * l399_BB_1-4: Broad Band beams in 2D MOT configuration
- Polarization conventions: CircularLeft, CircularRight, Horizontal

Atom-light interactions are added with specific detunings (in units of gamma) for each beam.
The atom_ligth couplings used here are lists of tuples to simulate the effects of a laser comb.

This configuration can be imported in other scripts as:

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from configurations.red_mot_configuration import configuration_MM

All beams are positioned relative to the origin (0,0,0). Directions are given as (x, y, z) unit vectors.

"""

from atomsmltr.atoms import Strontium
from atomsmltr.environment.lasers import GaussianLaserBeam
from atomsmltr.environment.lasers.polarization import CircularLeft, CircularRight, Horizontal
from atomsmltr.simulation import Configuration
import numpy as np


# --- Constants ---
TOTAL_RADIUS = 12e-3
NARROW_RADIUS = 4e-3
BROAD_RADIUS = TOTAL_RADIUS - NARROW_RADIUS
TOTAL_POWER = 6e-3
NARROW_POWER = 3e-3
BROAD_POWER = TOTAL_POWER - NARROW_POWER
PI = np.pi
FREQUENCY_GAP = 5
NUMBER_OF_DETUNINGS = 200//FREQUENCY_GAP
WEIGHT = 1 / (NUMBER_OF_DETUNINGS + 1)


# --- Atom Setup ---
atom = Strontium()
intercombination = atom.trans["intercombination"]
gamma = intercombination.Gamma


# --- Detunings ---
DETUNINGS_BB = [(-FREQUENCY_GAP*i*gamma*PI, WEIGHT) for i in range(1, NUMBER_OF_DETUNINGS)]
DETUNINGS_BB.append((-2*gamma*PI, WEIGHT))
DETUNINGS_NB = -2*gamma*PI


# --- Laser Setup ---
l399_NB_1 = GaussianLaserBeam(
    wavelength = 689e-9,
    waist = NARROW_RADIUS, 
    power = NARROW_POWER/4,
    waist_position = (0, 0, 0),
    direction = (0, 0, -1), 
    polarization = Horizontal(),
    tag = "laser_NB_1"
)


l399_NB_2 = GaussianLaserBeam(
    wavelength = 689e-9,
    waist = NARROW_RADIUS, 
    power = NARROW_POWER/4,
    waist_position = (0, 0, 0),
    direction = (0, 0, 1), 
    polarization = Horizontal(),
    tag = "laser_NB_2"
)


l399_NB_3 = GaussianLaserBeam(
    wavelength = 689e-9,
    waist = NARROW_RADIUS, 
    power = NARROW_POWER/4,
    waist_position = (0, 0, 0),
    direction = (1, 0, 0), 
    polarization = Horizontal(),
    tag = "laser_NB_3"
)


l399_NB_4 = GaussianLaserBeam(
    wavelength = 689e-9,
    waist = NARROW_RADIUS, 
    power = NARROW_POWER/4,
    waist_position = (0, 0, 0),
    direction = (-1, 0, 0), 
    polarization = Horizontal(),
    tag = "laser_NB_4"
)


l399_BB_1 = GaussianLaserBeam(
    wavelength = 689e-9,
    waist = BROAD_RADIUS, 
    power = BROAD_POWER/4,
    waist_position = (0, 0, 0),
    direction = (0, 0, -1), 
    polarization = Horizontal(),
    tag = "laser_BB_1"
)


l399_BB_2 = GaussianLaserBeam(
    wavelength = 689e-9,
    waist = BROAD_RADIUS, 
    power = BROAD_POWER/4,
    waist_position = (0, 0, 0),
    direction = (0, 0, 1), 
    polarization = Horizontal(),
    tag = "laser_BB_2"
)


l399_BB_3 = GaussianLaserBeam(
    wavelength = 689e-9,
    waist = BROAD_RADIUS, 
    power = BROAD_POWER/4,
    waist_position = (0, 0, 0),
    direction = (1, 0, 0), 
    polarization = Horizontal(),
    tag = "laser_BB_3"
)


l399_BB_4 = GaussianLaserBeam(
    wavelength = 689e-9,
    waist = BROAD_RADIUS, 
    power = BROAD_POWER/4,
    waist_position = (0, 0, 0),
    direction = (-1, 0, 0), 
    polarization = Horizontal(),
    tag = "laser_BB_4"
)


# --- Configuration Setup ---
configuration_NB_BB = Configuration()
configuration_NB_BB.atom = atom
configuration_NB_BB += l399_BB_1, l399_BB_2, l399_BB_3, l399_BB_4, l399_NB_1, l399_NB_2, l399_NB_3, l399_NB_4


# --- Add Atom-light Coupling ---
configuration_NB_BB.add_atomlight_coupling("laser_NB_1", "main", DETUNINGS_NB)
configuration_NB_BB.add_atomlight_coupling("laser_NB_2", "main", DETUNINGS_NB)
configuration_NB_BB.add_atomlight_coupling("laser_NB_3", "main", DETUNINGS_NB)
configuration_NB_BB.add_atomlight_coupling("laser_NB_4", "main", DETUNINGS_NB)
configuration_NB_BB.add_atomlight_coupling("laser_BB_1", "main", DETUNINGS_BB)
configuration_NB_BB.add_atomlight_coupling("laser_BB_2", "main", DETUNINGS_BB)
configuration_NB_BB.add_atomlight_coupling("laser_BB_3", "main", DETUNINGS_BB)
configuration_NB_BB.add_atomlight_coupling("laser_BB_4", "main", DETUNINGS_BB)











