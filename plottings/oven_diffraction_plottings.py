"""
This file is designed to plot and compare the initial distributions of atoms just outside of the oven.
The plotted graphs are:
- Distribution of axial velocities                      (1)
- Spread of v_x compared to v_y                         (2)
- Spread of initial position of the atoms at (z=0)      (3)
- Distribution of radial velocities                     (4)

Conventions used in this file:
- Axial direction: z-axis
- Radial (transverse) velocity: v_x = v_r*cos(θ), v_y = v_r*sin(θ)

The sampled initial conditions are compared to the theoretical distributions:

- Axial velocity distribution: Alexandre Dareau, "Manipulation cohérente d’un condensat 
  de Bose-Einstein d’ytterbium sur la transition 'd’horloge'", ENS Paris, 2015.
- Radial velocity distribution: P. T. Greenland et al., J. Phys. D: Appl. Phys., 1985, 18, 1223.
- Initial spatial spread: Gaussian in x and y within the aperture of diameter D centered at (0,0,0).

"""

import os 
import sys
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path


# --- Path Parameters ---
ROOT = Path(__file__).resolve().parent
DATA_FOLDER = ROOT.parent / "figures" # Folder to save plots
DATA_FOLDER.mkdir(parents=True, exist_ok=True)


# --- File parameters ---
SAVE = True
SHOW = False
FIG_NAME = "oven_diffraction_distributions"
SAVE_FMT = "pdf"
SAVE_DIR = DATA_FOLDER


# --- Figure Parameters ---
BLUE = "#1f77b4"
RED = "#d62728"
GREEN = "#2ca02c"
TEXT_WIDTH = 6
GOLDEN = 1.618
FIGSIZE = (TEXT_WIDTH, TEXT_WIDTH / GOLDEN * 1.5) # Adjusted for 2x2 grid
N_BIN_AXIAL = 50
N_BIN_RADIAL = 50



# --- Matplotlib Setup ---
def setup_style():
    """
    Configure matplotlib rcParams for consistent, publication-quality plots.
    """
    plt.rcParams.update({
        "font.family": "serif",
        "font.size": 10,
        "axes.labelsize": 10,
        "legend.fontsize": 8,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "lines.linewidth": 1.2,
        "figure.autolayout": True,
        "text.usetex": False # Set to True if you have LaTeX installed
    })


# --- Import Magnetic Configuration ---
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from configurations.oven_diffraction import u0_diffraction, ortho_distribution_grid, axial_distribution_grid, v_grid_axial, v_grid_ortho, D, N_ATOMS


# --- Main plotting Routine ---
def run_distribution_plots():
    """
    Generate and save a 2×2 figure comparing sampled distributions
    to their corresponding theoretical predictions.

    Panels:
        (0,0) Axial velocity histogram vs theory
        (0,1) Transverse velocity map (v_x, v_y)
        (1,0) Spatial cross-section at oven exit
        (1,1) Radial transverse velocity histogram vs theory
    """
    setup_style()

    # Create 2×2 figure with high DPI for clean PDF rasterization
    fig, axes = plt.subplots(2, 2, figsize=FIGSIZE, dpi=300)

    # Extract velocity components from phase-space array
    v_x = u0_diffraction[:, 3]
    v_y = u0_diffraction[:, 4]
    v_z = u0_diffraction[:, 5]

    # Radial transverse speed
    v_r = np.sqrt(v_x**2 + v_y**2)

    # Extract spatial coordinates
    x = u0_diffraction[:, 0]
    y = u0_diffraction[:, 1]
    z = u0_diffraction[:, 2]
    

    # --------------------------------------------------------------
    # 1. Axial velocity distribution
    # --------------------------------------------------------------
    ax0 = axes[0, 0]
    _, bins_axial, _ = ax0.hist(
        v_z, 
        bins=N_BIN_AXIAL, 
        density=False, 
        alpha=0.6, 
        color='blue', 
        label='Sampled'
        )
    
    # Convert PDF into expected counts per bin
    bin_width_axial = bins_axial[1] - bins_axial[0]
    scale_factor_axial = N_ATOMS * bin_width_axial

    # Normalize theoretical axial PDF
    axial_pdf = axial_distribution_grid / np.trapezoid(
        axial_distribution_grid,
        v_grid_axial
        )
    
    # Overlay theoretical prediction
    ax0.plot(
        v_grid_axial, 
        axial_pdf * scale_factor_axial, 
        color=RED, 
        lw=1.5, 
        label='Theory'
        )
    
    ax0.set_title(r"Axial Velocity $v_z$")
    ax0.set_ylabel("Number of Atoms")
    ax0.set_xlabel("Velocity (m/s)")
    ax0.legend()


    # --------------------------------------------------------------
    # 2. Transverse velocity map (v_x, v_y)
    # --------------------------------------------------------------
    ax1 = axes[0, 1]

    # Scatter plot of transverse velocities
    ax1.scatter(
        v_x, 
        v_y, 
        s=0.5, 
        alpha=0.3, 
        color=BLUE, 
        rasterized=True  # critical for PDF file size
        ) 
    
    ax1.set_title(r"Transverse Map $v_x, v_y$")

    # Force equal scaling to preserve circular symmetry
    ax1.set_box_aspect(1) 


    # --------------------------------------------------------------
    # 3. Spatial cross-section at oven exit
    # --------------------------------------------------------------
    ax2 = axes[1, 0]

    # Scatter plot of transverse positions (converted to mm)
    ax2.scatter(
        x*1e3, 
        y*1e3, 
        s=0.5, 
        color=GREEN, 
        alpha=0.4, 
        rasterized=True
        )
    
    # Overlay physical aperture boundary
    circle = plt.Circle(
        (0, 0), 
        (D/2)*1e3, 
        color=RED, 
        fill=False, 
        lw=1.2, 
        linestyle='--', 
        zorder=5
        )
    ax2.add_patch(circle)

    ax2.set_title("Spatial Cross-section")

    # Preserve circular geometry
    ax2.set_box_aspect(1)


    # --------------------------------------------------------------
    # 4. Radial transverse velocity distribution
    # --------------------------------------------------------------    
    ax3 = axes[1, 1]

    # Histogram of radial transverse speeds
    _, bins_radial, _ = ax3.hist(
        v_r, 
        bins=N_BIN_RADIAL, 
        density=False, 
        alpha=0.6, 
        color='blue', 
        label='Sampled'
        )
    
    # Convert PDF into histogram counts
    bin_width_radial = bins_radial[1] - bins_radial[0]
    scale_factor_radial = N_ATOMS * bin_width_radial

    # Normalize transverse PDF
    ortho_pdf = ortho_distribution_grid / np.trapezoid(
        ortho_distribution_grid, 
        v_grid_ortho
        )
    
    # Overlay radial distribution
    ax3.plot(
        v_grid_ortho, 
        ortho_pdf * scale_factor_radial, 
        color=RED, 
        lw=1.5, 
        zorder=3
        )
    
    ax3.set_title(r"Radial Velocity $v_r$")
    ax3.set_ylabel("Number of Atoms")
    ax3.set_xlabel("Velocity (m/s)")
    ax3.legend()

    # --------------------------------------------------------------
    # Final formatting
    # --------------------------------------------------------------
    for ax in axes.flatten():
        ax.grid(True, alpha=0.2, linestyle=':', zorder=0)
        ax.legend(frameon=False)

        # Axis labeling depends on panel content
        ax.set_xlabel("m/s" if ax in [ax0, ax3] else "mm")

    # Tight layout to prevent clipping in PDF output
    plt.tight_layout(pad=1.2)
    

    # Save or display figure
    if SAVE:
        plt.savefig(
            SAVE_DIR / f"{FIG_NAME}.{SAVE_FMT}"
            , bbox_inches='tight', 
            transparent=False
            )
        
    if SHOW:
        plt.show()


# --- MAIN ---
if __name__ == "__main__":
    run_distribution_plots()
