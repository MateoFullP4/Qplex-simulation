import os
import sys
import importlib.util
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# --- Path Parameters ---
ROOT = Path(__file__).resolve().parent
BASE_DIR = ROOT.parent
# Ensure the folder name matches your OS ("Data_simulations" vs "data simulations")
DATA_DIR = BASE_DIR / "data" / "atomic_flux_fraction"
SAVE_DIR = BASE_DIR / "figures" / "atomic_flux_fraction"
SAVE_DIR.mkdir(parents=True, exist_ok=True)

# --- File parameters ---
SAVE = True
SHOW = True
FIG_NAME = "atomic_flux_fraction"
SAVE_FMT = "pdf"

# --- Figure Parameters ---
ORANGE = "#ff7f0e"
TEXT_WIDTH = 6
GOLDEN = 1.618
FIGSIZE = (TEXT_WIDTH, TEXT_WIDTH / GOLDEN)

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

# --- Data Fetching ---
def fetch_data():
    all_mean_rates = {}
    all_detunings = {}
    
    # Get and sort the .py files numerically
    data_files = sorted(
        [f for f in os.listdir(DATA_DIR) if f.startswith('detunings_') and f.endswith('.py')],
        key=lambda x: int(x.split('_')[1].split('.')[0])
    )

    for filename in data_files:
        file_path = DATA_DIR / filename
        module_name = filename[:-3] 
        
        spec = importlib.util.spec_from_file_location(module_name, str(file_path))
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)

        if hasattr(module, 'detunings') and hasattr(module, 'mean_rates'):
            all_detunings[module_name] = module.detunings
            all_mean_rates[module_name] = module.mean_rates
        else:
            print(f"Warning: Data missing in {filename}")
            
    return all_detunings, all_mean_rates


# --- Data Processing ---
def concat_and_sort(dic):
    combined_list = []
    for key in dic:
        combined_list.extend(dic[key])
    return sorted(combined_list) # Returns a new sorted list

def interleave_mean_rates(dic):
    # Sort keys to ensure we interleave in order: det_0[0], det_1[0], det_2[0]...
    sorted_keys = sorted(dic.keys(), key=lambda x: int(x.split('_')[1]))
    list_of_lists = [dic[k] for k in sorted_keys]
    
    interleaved = []
    for group in zip(*list_of_lists):
        interleaved.extend(group)
    return interleaved


# --- Running Function --- 
def run_plot():
    setup_style()
    
    # 1. Fetch
    all_det, all_rates = fetch_data()
    
    if not all_det:
        print("No data found!")
        return

    # 2. Process
    # We use concat_and_sort for detunings (assuming they define the x-axis)
    # We use interleave for rates to match the x-axis structure
    det_sorted = concat_and_sort(all_det)
    rates_interleaved = interleave_mean_rates(all_rates)

    # 3. Plot
    fig, ax = plt.subplots(figsize=FIGSIZE)
    
    ax.plot(det_sorted, rates_interleaved, color=ORANGE, marker='o', markersize=4, label="Mean Flux Fraction")
    
    ax.set_xlabel("Detuning (Units)")
    ax.set_ylabel("Mean Rates")
    ax.set_title("Atomic Flux Fraction vs Detuning")
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.legend()

    plt.tight_layout()
    
    if SAVE:
        plt.savefig(SAVE_DIR / f"{FIG_NAME}.{SAVE_FMT}")
        print(f"Plot saved to {SAVE_DIR}")
    if SHOW:
        plt.show()

if __name__ == "__main__":
    run_plot()