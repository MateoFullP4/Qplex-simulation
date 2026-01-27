import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


# --- Path Parameters ---
ROOT = Path(__file__).resolve().parent
DATA_FOLDER = ROOT.parent / "figures" # Folder to save plots
DATA_FOLDER.mkdir(parents=True, exist_ok=True)

# --- File parameters ---
SAVE = True
SHOW = False
FIG_NAME = "magnetic_field_profile"
SAVE_FMT = "pdf"
SAVE_DIR = DATA_FOLDER


# --- Figure Parameters ---
ORANGE = "#ff7f0e"
TEXT_WIDTH = 6
GOLDEN = 1.618
FIG_W = TEXT_WIDTH * 1
FIG_H = FIG_W / GOLDEN * 1.2 
FIGSIZE = (TEXT_WIDTH, TEXT_WIDTH * 1.2) # Adjusted for 3x2 grid
FIG = None
AX = None


# --- Matplotlib Setup ---
def setup_style():
    """Sets basic professional plot parameters"""
    plt.rcParams.update({
        "font.family": "serif",
        "font.size": 10,
        "axes.labelsize": 10,
        "legend.fontsize": 8,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "lines.linewidth": 1.5,
        "figure.autolayout": True
    })

# --- Grid Definitions ---
TX = np.linspace(-0.050, 0.050, 101)
TY = np.linspace(-0.050, 0.050, 101)
TZ = np.linspace(-0.250, 0.400, 651)

gridYZ = np.array([[(0, y, z) for y in TY] for z in TZ])
gridXZ = np.array([[(x, 0, z) for x in TX] for z in TZ])
gridXY = np.array([[(x, y, 0) for x in TX] for y in TY])


# --- Import Magnetic Configuration ---
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from configurations.magnetic_configuration import magnetSet_asymmetric

# --- Data Generation ---
def compute_B_slice(grid):
    """ 
    Compute B-field for a (N1,N2,3) grid for Magpylib.
    Args:
        grid (): 
    """
    flat = grid.reshape(-1, 3)
    Bflat = magnetSet_asymmetric.getB(flat)
    return Bflat.reshape(grid.shape)  


def get_field_data():
    # Resolution
    res_2d = 100
    tx = np.linspace(-0.05, 0.05, res_2d)
    ty = np.linspace(-0.05, 0.05, res_2d)
    tz = np.linspace(-0.20, 0.20, res_2d) # Standardized for slices

    # 2D Grids
    gridXY = np.array([[(x, y, 0) for x in tx] for y in ty])
    gridXZ = np.array([[(x, 0, z) for x in tx] for z in tz])
    gridYZ = np.array([[(0, y, z) for y in ty] for z in tz])

    # 1D Lines
    lx = np.linspace(-0.05, 0.05, 500)
    lz = np.linspace(-0.25, 0.40, 500)
    
    pts_x = np.column_stack([lx, np.zeros_like(lx), np.zeros_like(lx)])
    pts_y = np.column_stack([np.zeros_like(lx), lx, np.zeros_like(lx)])
    pts_z = np.column_stack([np.zeros_like(lz), np.zeros_like(lz), lz])

    # Compute
    B_XY = np.linalg.norm(compute_B_slice(gridXY), axis=2)
    B_XZ = np.linalg.norm(compute_B_slice(gridXZ), axis=2)
    B_YZ = np.linalg.norm(compute_B_slice(gridYZ), axis=2)
    
    B_x = np.linalg.norm(magnetSet_asymmetric.getB(pts_x), axis=1)
    B_y = np.linalg.norm(magnetSet_asymmetric.getB(pts_y), axis=1)
    B_z = np.linalg.norm(magnetSet_asymmetric.getB(pts_z), axis=1)

    return (tx, ty, tz, lx, lz), (B_XY, B_XZ, B_YZ, B_x, B_y, B_z)

def run_plots():
    grids, fields = get_field_data()
    tx, ty, tz, lx, lz = grids
    B_XY, B_XZ, B_YZ, B_x, B_y, B_z = fields

    fig, axes = plt.subplots(3, 2, figsize=FIGSIZE)
    setup_style()

    # --- 2D SLICES (Left Column) ---
    cmap = 'magma'
    
    # XY
    im0 = axes[0, 0].imshow(B_XY, extent=[tx.min()*1e3, tx.max()*1e3, ty.min()*1e3, ty.max()*1e3], origin='lower', cmap=cmap)
    axes[0, 0].set_title(r"$|B(x,y,0)|$")
    plt.colorbar(im0, ax=axes[0, 0])

    # XZ
    im1 = axes[1, 0].imshow(B_XZ, extent=[tx.min()*1e3, tx.max()*1e3, tz.min()*1e3, tz.max()*1e3], origin='lower', cmap=cmap, aspect='auto')
    axes[1, 0].set_title(r"$|B(x,0,z)|$")
    plt.colorbar(im1, ax=axes[1, 0])

    # YZ
    im2 = axes[2, 0].imshow(B_YZ, extent=[ty.min()*1e3, ty.max()*1e3, tz.min()*1e3, tz.max()*1e3], origin='lower', cmap=cmap, aspect='auto')
    axes[2, 0].set_title(r"$|B(0,y,z)|$")
    plt.colorbar(im2, ax=axes[2, 0])

    # --- 1D PROFILES (Right Column) ---
    axes[0, 1].plot(lx*1e3, B_x, color=ORANGE)
    axes[0, 1].set_title("Along X-axis")
    
    axes[1, 1].plot(lx*1e3, B_y, color=ORANGE)
    axes[1, 1].set_title("Along Y-axis")
    
    axes[2, 1].plot(lz*1e3, B_z, color=ORANGE)
    axes[2, 1].set_title("Along Z-axis")

    # Labels and Grid
    for i, row in enumerate(axes):
        for j, ax in enumerate(row):
            ax.grid(True, alpha=0.3, linestyle='--')
            ax.set_xlabel("mm")
            ax.set_ylabel("T")

    plt.tight_layout()
    if SAVE:
        plt.savefig(SAVE_DIR / f"{FIG_NAME}.{SAVE_FMT}")
    if SHOW:
        plt.show()

if __name__ == "__main__":
    run_plots()