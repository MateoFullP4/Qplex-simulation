"""

"""

import numpy as np
import sys
import os
from atomsmltr.simulation import RK4
from pathlib import Path
from atomsmltr.atoms import Ytterbium, Strontium
from scipy import constants as csts


# --- Import Magnetic Configuration ---
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from configurations.red_mot_configuration import configuration_MM
from configurations.oven_diffraction import u0_diffraction, N_ATOMS


# --- Result File ---
ROOT = Path(__file__).resolve().parent
DATA_FOLDER = ROOT.parent / "data" / "atomic_flux_fraction"
SAVE = False

# --- Simulation Parameters ---
STEPS_NUMBER = 200
TOTAL_DURATION = 0.1
MASS = Strontium().mass
KB = csts.Boltzmann
TRANS = Strontium().trans["main"]
GAMMA = TRANS.Gamma
DETUNINGS = [(-0.5*i*GAMMA + 0.5) for i in range(2, 40)]


# --- Setup Initial Conditions ---
T = np.linspace(0, TOTAL_DURATION, STEPS_NUMBER)  


# --- Catching Criteria ---
def catching_rate(traj):
    z_trajs = traj[:, 2, :]  # Shape: (n_atoms, n_timesteps)
    cutoff = 0.35

    # --- Step 1: Determine cutoff mask ---
    # Find the first index where z > cutoff for each atom
    # np.argmax returns the first index where the condition is True
    # If it's never True, it returns 0, so we handle that case.
    greater_than_cutoff = z_trajs > cutoff
    cut_indices = np.argmax(greater_than_cutoff, axis=1)
    
    # If argmax is 0 but the value isn't > cutoff, it means it never crossed
    never_crossed = ~np.any(greater_than_cutoff, axis=1)
    cut_indices[never_crossed] = z_trajs.shape[1]

    # Create mask using broadcasting
    time_indices = np.arange(z_trajs.shape[1])
    mask = time_indices < cut_indices[:, np.newaxis] 
    
    # --- Step 2: Determine final state ---
    z_final = z_trajs[:, -1]
    z_max = np.max(z_trajs, axis=1)
    
    # --- Step 3: Apply trapping criteria ---
    # Using the mask to see if the atom was cut off before the end
    escaped = (z_final > 0.05) | (z_final < -0.05) | (cut_indices < z_trajs.shape[1])
    lost    = (~escaped) & (z_max > 0.01)
    trapped = ~(escaped | lost)

    # --- Step 4: Compute fraction trapped ---
    rate = np.mean(trapped)

    return rate


# --- Change Existing Detunings ---
def clear_detunings():
    configuration_MM.rm_atomlight_coupling("main","laser_ZS_1")


def add_detunings(detuning):
    configuration_MM.add_atomlight_coupling("laser_ZS_1", "main", detuning)
    


# --- Running Block ---
def run_simu():
    sim_ytterbium = RK4(configuration_MM)
    coll = sim_ytterbium._integrate(u0_diffraction, T)
    trajectories = coll.y
    return trajectories


# --- Detunings Sweeping ---
def detunings_sweeping():
    rates = []
    for i, det in enumerate(DETUNINGS):
        clear_detunings()
        add_detunings(DETUNINGS[i])
        traj = run_simu()
        rate = catching_rate(traj)
        rates.append(float(rate))
        print(i)

        if SAVE:
            gamma_multiplier = det / GAMMA
            folder_name = f"catching_rate_{gamma_multiplier:.2f}Gamma"
            save_dir = DATA_FOLDER / folder_name
            
            # Create directory if it doesn't exist
            save_dir.mkdir(parents=True, exist_ok=True)
            
            # 3. Saving Data
            file_path = save_dir / "sim_data.npz"
            np.savez(file_path, traj=traj, rate=rate, detuning=det)
            
            print(f"Iteration {i}: Saved to {folder_name}")
    
    return rates

print(detunings_sweeping())