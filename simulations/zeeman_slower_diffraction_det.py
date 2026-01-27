import numpy as np
import sys
import os
from atomsmltr.simulation import RK4
from pathlib import Path
from atomsmltr.atoms import Ytterbium
from scipy import constants as csts


# --- Import Magnetic Configuration ---
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from configurations.red_mot_configuration import configuration_MM
from configurations.oven_diffraction import u0_diffraction, N_ATOMS


# --- Result File ---
ROOT = Path(__file__).resolve().parent
DATA_FOLDER = ROOT.parent / "data simulations" / "atomic_flux_fraction"
SIM_RESULT_FILE = DATA_FOLDER/"zeeman_slower_diffraction_det.npz"


# --- Simulation Parameters ---
STEPS_NUMBER = 500
TOTAL_DURATION = 0.1
MASS = Ytterbium().mass
KB = csts.Boltzmann


# --- Physical Functions ---
def radial_temperature(atom_array):
    vx = atom_array[:, 3, :]
    vy = atom_array[:, 4, :]
    v_sq = vx**2 + vy**2
    T_radial = MASS * np.mean(v_sq, axis=0) / (2 * KB)
    
    return T_radial


def axial_temperature(atom_array):
    vz = atom_array[:, 5, :]
    v_sq = vz**2    
    T_axial = MASS * np.mean(v_sq, axis=0) / (2 * KB)
    
    return T_axial


# --- Setup Initial Conditions ---
t = np.linspace(0, TOTAL_DURATION, STEPS_NUMBER)  


# --- Setup Simulation ---
sim_ytterbium = RK4(configuration_MM)
coll = sim_ytterbium._integrate(u0_diffraction, t)
trajectories = coll.y


# --- Catching Criteria ---
def catching_rate(traj):
    """
    Check the last state of the trajectories and calculates the proportion of atoms being caught by the MOT. 
    Args:
        - traj (ndarray): Array of shape (n_atoms, 6, len(t)) returned by _integrate.
    
    Returns:
        - rate (float) : proportion of atomic being caught in the MOT.
    """
    z_trajs = traj[:, 2, :]

    #  Step 1: Determine cutoff mask 
    cutoff = 0.35

    # For each atom, find the first index where z > cutoff
    cut_indices = np.searchsorted(z_trajs, cutoff, side='right')

    # Create mask array: True for timesteps to consider
    # Each row is an atom, each column a timestep
    time_indices = np.arrange(z_trajs.shape[1])
    mask = time_indices < cut_indices[:, 'None']    # shape (n_atoms, len(t))
    
    # Step 2: Determine final state for each atom 
    z_final = z_trajs[:, -1]
    z_max = np.max(z_trajs, axis=1)  # maximum z 
    
    # --- Step 3: Apply trapping criteria ---
    escaped = (z_final > 0.05) | (z_final < -0.05) | (~mask[:, -1])
    lost    = (~escaped) & (z_max > 0.01)
    trapped = ~(escaped | lost)


    # --- Step 4: Compute fraction trapped ---
    rate = np.mean(trapped)

    return rate

# --- Simulation ---
print(catching_rate(trajectories))


# --- Save Data ---

# np.savez(SIM_RESULT_FILE, detuning=detuning, rates=rates)


