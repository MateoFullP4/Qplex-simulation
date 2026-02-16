"""
This file is designed to simulate atomic trajectories and compute the axial
force map for a Zeeman slower detuning in the QPlex setup.

Workflow:

1. Define a grid of initial axial positions and velocities for the atoms.
2. Generate a list of initial conditions `u0_list` with different initial velocities.
3. Integrate the atomic trajectories over time using `RK4` from atomsmltr.
4. Compute the axial force $F_z$ along a (z, v_z) grid using the MOT configuration.
5. Save the resulting trajectories and force map to a `.npz` file in
   `data/force_map_trajectories/`.

Conventions and notes:

- Axial direction is along the z-axis.
- Position coordinates are in meters; velocities in m/s.
- Force is computed in Newtons.
- The force map grid is defined using `Z_GRID` and `VZ_GRID`.
- Initial positions for each atom are set to the origin in x and y,
  with z spanning from -0.15 to 0.35 m.
- `u0_list` defines atoms starting at z = -0.15 m with velocities
  spanning from 0 to 550 m/s.

Purpose:

- Visualize the axial force landscape $F_z(z,v_z)$.
- Provide trajectory data for plotting and analysis of trapped, lost, and
  escaped atoms.
- Facilitate understanding and optimization of Zeeman slower parameters.

Dependencies:

- `configuration_MM` must be defined in `configurations.blue_mot_configuration`
  and compatible with `atomsmltr.RK4`.
- NumPy is used for array handling.
"""


import numpy as np
import sys
import os
from atomsmltr.simulation import RK4
from pathlib import Path


# --- Import Magnetic Configuration ---
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from configurations.blue_mot_configuration import configuration_MM


# --- Result File ---
ROOT = Path(__file__).resolve().parent
DATA_FOLDER = ROOT.parent / "data" / "force_map_trajectories"
SIM_RESULT_FILE = DATA_FOLDER/"zeeman_slower_detuning_12.5.npz"
SAVE = True


# --- Setup Initial Conditions ---
t = np.linspace(0, 0.1, 5000)
v_list = np.linspace(0, 550, 56)
u0_list = [(0, 0, -0.15, 0, 0, vz) for vz in v_list]
z_axis  = np.linspace(-0.15, 0.35, 200)        
vz_axis = np.linspace(-150, 550, 351)    
Z_GRID, VZ_GRID = np.meshgrid(z_axis, vz_axis)


pos = np.column_stack([
    np.zeros(Z_GRID.size),     
    np.zeros(Z_GRID.size),     
    Z_GRID.ravel(),            
    np.zeros(Z_GRID.size),     
    np.zeros(Z_GRID.size),     
    VZ_GRID.ravel(),           
])

# --- Setup Simulation ---
sim_ytterbium = RK4(configuration_MM)
coll = sim_ytterbium._integrate(u0_list, t)
coll = coll.y
z = coll[:, 2, :]
vz = coll[:, 5, :]
force = sim_ytterbium.get_force(pos)     
Fz = force[:, 2].reshape(Z_GRID.shape)


# --- Save Data ---
if SAVE:
    np.savez(SIM_RESULT_FILE, t=t, z=z, vz=vz, Fz=Fz)

