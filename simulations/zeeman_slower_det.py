import numpy as np
import sys
import os
from atomsmltr.simulation import RK4
from pathlib import Path


# --- Import Magnetic Configuration ---
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from configurations.red_mot_configuration import configuration_MM


# --- Result File ---
ROOT = Path(__file__).resolve().parent
DATA_FOLDER = ROOT.parent / "data simulations"
SIM_RESULT_FILE = DATA_FOLDER/"zeeman_slower.npz"
SAVE = False


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

