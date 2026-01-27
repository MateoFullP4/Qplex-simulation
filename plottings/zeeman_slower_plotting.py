import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


# --- Path Parameters ---
ROOT = Path(__file__).resolve().parent
DATA_FOLDER = ROOT.parent / "figures"
DATA_FOLDER.mkdir(parents=True, exist_ok=True)
SIM_DATA_FILE = ROOT.parent / "data simulations" / "zeeman_slower.npz"

# --- File parameters ---
SAVE = True
SHOW = True
FIG_NAME = "force_map_trajectories"
SAVE_FMT = "pdf"
SAVE_DIR = DATA_FOLDER


# --- Figure Parameters ---
TEXT_WIDTH = 6
GOLDEN = 1.618
FIGSIZE = (TEXT_WIDTH, TEXT_WIDTH / GOLDEN * 1.2)
COLOR_TRAP = "C0"  # Blue for trapped
COLOR_LOSS = "C3"  # Red for through-put
COLOR_ESCAPED = "k" # Black for escaped


# --- Matplotlib Setup ---
def setup_style():
    """Sets professional plot parameters"""
    plt.rcParams.update({
        "font.family": "serif",
        "font.size": 10,
        "axes.labelsize": 10,
        "legend.fontsize": 8,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "lines.linewidth": 1.0,
        "figure.autolayout": True
    })


def load_data():
    """Loads simulation data from the npz file"""
    if not SIM_DATA_FILE.exists():
        raise FileNotFoundError(f"Data file not found at {SIM_DATA_FILE}. Run simulation first.")
    
    data = np.load(SIM_DATA_FILE)
    return data['z'], data['vz'], data['Fz']

def run_plots():
    setup_style()
    z_trajs, vz_trajs, Fz_grid = load_data()

    # Reconstruct grids for plotting (matching the sim logic)
    z_axis = np.linspace(-0.15, 0.35, 200)
    vz_axis = np.linspace(-150, 550, 351)
    Z_GRID, VZ_GRID = np.meshgrid(z_axis, vz_axis)

    fig, ax = plt.subplots(figsize=FIGSIZE)

    # --- 1) Plot Force Map (Background) ---
    # Fz units are typically ~10^-19 N based on your snippet
    pcm = ax.pcolormesh(
        Z_GRID * 100, # cm
        VZ_GRID,      # m/s
        Fz_grid,
        cmap="bwr",
        shading="auto",
        vmin=-4.5e-19,
        vmax=4.5e-19,
        zorder=0
    )
    cb = fig.colorbar(pcm, ax=ax)
    cb.set_label(r"$F_z$ (N)")

    # --- 2) Plot Trajectories (Foreground) ---
    # z_trajs shape is (N_atoms, Time_Steps)
    for i in range(z_trajs.shape[0]):
        z = z_trajs[i]
        vz = vz_trajs[i]

        # Trajectory masking logic
        cut_idx = np.searchsorted(z, 0.35)
        if cut_idx >= len(z):
            mask = np.ones_like(z, dtype=bool)
        else:
            mask = np.zeros_like(z, dtype=bool)
            mask[:cut_idx] = True

        # Color coding logic
        if z[-1] > 0.05 or z[-1] < -0.05 or not mask[-1]:
            color = COLOR_ESCAPED
        elif np.any(z > 0.01):
            color = COLOR_LOSS
        else:
            color = COLOR_TRAP
        
        ax.plot(z[mask] * 100, vz[mask], color=color, alpha=0.7, linewidth=0.8, zorder=1)

    # --- Final Formatting ---
    ax.set_xlim(-15, 35)
    ax.set_ylim(-150, 550)
    ax.set_xlabel("Position $z$ (cm)")
    ax.set_ylabel("Speed $v_z$ (m/s)")
    ax.set_title("Force Map and Atomic Trajectories")
    ax.grid(True, alpha=0.3, linestyle='--')

    if SAVE:
        save_path = SAVE_DIR / f"{FIG_NAME}.{SAVE_FMT}"
        plt.savefig(save_path, bbox_inches='tight')
        print(f"Figure saved to {save_path}")
    
    if SHOW:
        plt.show()

if __name__ == "__main__":
    run_plots()