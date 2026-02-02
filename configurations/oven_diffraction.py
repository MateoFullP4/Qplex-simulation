"""
This file is designed to create the initial distribution for atoms right outside of the oven. 

In this file, the convention are such as :
- Axial : z-axis
- Radial : v_x = v_r*cos(θ), v_y = v_r*sin(θ)

The initial conditions follow the distributions found in these papers : 

- Axial Distribution : Alexandre Dareau. Manipulation cohérente d’un condensat de Bose-Einstein d’ytterbium sur la transition ”d’horloge” : de la spectroscopie au magnétisme artificiel. Physique [physics]. Ecole normale supérieure - ENS PARIS, 2015. Français. ffNNT : 2015ENSU0018ff. fftel-01194429v3
- Radial Distribution : P T Greenland et al 1985 J. Phys. D: Appl. Phys. 18 1223, equation (13)
- Initial Spread : Gaussian Variables along x-axis and y-axis inside the hole of diameter D centered in (0, 0, 0).

To use this configuration in an other folder, you can use this syntax : 

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from configurations.oven_diffraction import u0_diffraction

Note that you can not input V_MIN = 0 as it is physically impossible, therefore leads to division by 0 in the code.
"""

import numpy as np
import random as rd
from scipy import constants as csts
from atomsmltr.atoms import Ytterbium, Rubidium, Strontium
import matplotlib.pyplot as plt


# --- Constants Definition ---
N_ATOMS = 1000
OVEN_TEMPERATURE = 823
MASS = Strontium().mass
KB = csts.Boltzmann
SIGMA = np.sqrt((KB*OVEN_TEMPERATURE)/MASS)
PI = np.pi
V_MIN_AXIAL = 5         # Minimum axial speed reached by the atoms right outside the oven, used for computation purposes
V_MAX_AXIAL = 1500       # Maximum axial speed reached by the atoms right outside the oven
V_MIN_ORTHO = 0.00001     # Minimum radial speed reached by the atoms right outside the oven
V_MAX_ORTHO = 50     # Maximum radial speed reached by the atoms right outside the oven
D = 0.0004                 # Diameter of the slit (oven)
L = 0.01                   # Length of the slit (oven)
SIGMA_X = 15e-3         # Standard deviation of the atoms along x-axis at z=0
SIGMA_Y = 15e-3         # Standard deviation of the atoms along y-axis at z=0


# --- Distribution Functions ---
def orthogonal_distribution(v):
    """
    Transverse (orthogonal) velocity distribution after the oven aperture.
    IMPORTANT:
        This function represents a 1D distribution in v_r only.
        It does NOT include the 2π v_r Jacobian factor required for a 2D radial probability density.

    Args:
        - v (float or ndarray): Transverse (radial) velocity v_r.

    Returns:
        - pdf (float or ndarray): Transverse velocity probability density f(v_r), normalized.
    """
    s = D/L
    a1 = np.exp(((-v**2)/(2*SIGMA**2))*((1+s**2)/s**2))
    a2 = np.sqrt(2*np.pi)*(np.sqrt(1+s**2) - 1)
    a3 = (SIGMA*s**3)/v

    return (a1/a2)*a3


def axial_distribution(v):
    """
    Axial velocity distribution for atoms exiting the oven.
    Physically, this describes the flux-weighted velocity distribution of atoms that successfully pass through the aperture.

    Args:
        - v (float or ndarray): Axial velocity v_z.

    Returns:
        - pdf (float or ndarray): Axial velocity probability density f(v_z), normalized.
    """
    a1 = np.sqrt(MASS/(2*PI*KB*OVEN_TEMPERATURE))
    a2 = np.exp((-MASS*v**2)/(2*KB*OVEN_TEMPERATURE))
    a3 = 1-np.exp(-(MASS*(v**2)*(D**2))/(2*KB*OVEN_TEMPERATURE*(L**2)))

    return a1*a2*a3


# --- Continuous to Discret ---
def build_cdf_axial(pdf_func, n_grid=200000): 
    """
    Returns the cumlative distributive functions of the axial PDF on a grid, ranging from V_MIN to V_MAX
    Args:
        - pdf_func (func) : probability density function 
        - n_grid (int) : number of points on the grid used to calculate the CDF
    
    Returns:
        - v (ndarray): Velocity grid spanning [V_MIN_AXIAL, V_MAX_AXIAL].
        - cdf (ndarray): Normalized cumulative distribution function on the grid.
    """
    
    v = np.geomspace(V_MIN_AXIAL, V_MAX_AXIAL, n_grid)
    pdf = pdf_func(v)
    pdf /= np.trapezoid(pdf, v)     # Normalization of the PDF
    dv = np.diff(v)
    cdf = np.cumsum(pdf[:-1] * dv)
    cdf = np.insert(cdf, 0, 0.0)    # Enforce CDF(0) = 0
    cdf /= cdf[-1]                  # Renormalize to ensure CDF[-1] = 1 (guards against numerical drift)

    return v, cdf


def build_cdf_ortho(pdf_func, n_grid=200000): 
    """
    Returns the cumlative distributive functions of the radial PDF on a grid, ranging from V_MIN to V_MAX
    Args:
        - pdf_func (func) : probability density function 
        - n_grid (int) : number of points on the grid used to calculate the CDF
    
    Returns:
        - v (ndarray): Velocity grid spanning [V_MIN_ORTHO, V_MAX_ORTHO].
        - cdf (ndarray): Normalized cumulative distribution function on the grid.

    """
    v = np.geomspace(V_MIN_ORTHO, V_MAX_ORTHO, n_grid)
    pdf = pdf_func(v)
    pdf /= np.trapezoid(pdf, v)
    dv = np.diff(v)
    cdf = np.cumsum(pdf[:-1] * dv)
    cdf = np.insert(cdf, 0, 0.0)
    cdf /= cdf[-1]

    return v, cdf


# --- Sampling Functions ---
def sample_velocities():
    """
    Sample initial atomic velocities from the axial and transverse distributions using inverse transform sampling.

    Returns:
        - velocities (ndarray): Array of shape (N_ATOMS, 3) containing (v_x, v_y, v_z).
    """
    v_r_grid, cdf_r = build_cdf_ortho(orthogonal_distribution)
    v_z_grid, cdf_z = build_cdf_axial(axial_distribution)

    # Draw uniform random numbers and map them through the inverse CDF
    v_r = np.interp(np.random.rand(N_ATOMS), cdf_r, v_r_grid)
    v_z = np.interp(np.random.rand(N_ATOMS), cdf_z, v_z_grid)

    # Uniform azimuthal angle ensures isotropy in the transverse plane
    theta = 2 * np.pi * np.random.rand(N_ATOMS)

    # Convert to Cartesian
    v_x = v_r * np.cos(theta)
    v_y = v_r * np.sin(theta)

    velocities = np.column_stack((v_x, v_y, v_z))

    return velocities



# --- Sampling Initial Positions ---
def sample_positions():
    """
    Sample initial atomic positions at the oven exit.
    Positions are drawn from a 2D Gaussian transverse distribution and truncated by a circular aperture of diameter D.

    Returns:
        - positions (ndarray): Array of shape (N_ATOMS, 3) containing (x, y, 0), the hole being in the (z=0) plan.
    """
    R = D / 2
    positions = np.zeros((N_ATOMS, 3))  

    positions[:, 2] = -0.15

    count = 0
    while count < N_ATOMS:

        n_to_sample = N_ATOMS - count   # Number of atoms still needed

        # Oversample to compensate for rejection at the aperture
        batch_size = int(n_to_sample * 1.5) + 10

        # Draw transverse positions from Gaussian source profile
        x = np.random.normal(0, SIGMA_X, batch_size)
        y = np.random.normal(0, SIGMA_Y, batch_size)

        # Apply circular aperture condition
        mask = x**2 + y**2 <= R**2

        accepted_x = x[mask]
        accepted_y = y[mask]

        # Number of accepted atoms in this batch
        n_accept = min(len(accepted_x), n_to_sample)

        # Store accepted positions
        positions[count:count+n_accept, 0] = accepted_x[:n_accept]
        positions[count:count+n_accept, 1] = accepted_y[:n_accept]
        count += n_accept

    return positions


# --- Setup the Distributions for Plottings ---
v_grid_axial = np.linspace(V_MIN_AXIAL, V_MAX_AXIAL, 1000)
axial_distribution_grid = axial_distribution(v_grid_axial)

v_grid_ortho = np.linspace(V_MIN_ORTHO, V_MAX_ORTHO, 1000)
ortho_distribution_grid = orthogonal_distribution(v_grid_ortho)


# --- Sampling Initial Conditions ---
v0_diffraction = sample_velocities()
r0_diffraction = sample_positions() 
u0_diffraction = np.hstack([r0_diffraction, v0_diffraction])

