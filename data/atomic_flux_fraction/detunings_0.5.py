"""
The data files in the "atomic_flux_fraction" subfolder come from running the "Zeeman_slower_diffraction.py" file in the "simulation
folder multiple times. I didn't push the .npz files on the Github (which can be enabled through the SAVE variable for memory size purpose, 
so I just copy-pasted it in these files.

Each file of this subfolder thus contains 10 lists of the finak catching rate depending on the detunings. 
In this file, the detunings are multiple of Gamma. 

You can obviously run the code yourself to validate the data, or use it as it is here. 
"""


import numpy as np
import matplotlib.pyplot as plt

detunings = [(0.5*i) for i in range(2, 40)]

l1 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.001, 0.003, 0.01, 0.016, 0.017, 0.015, 0.02, 0.012, 0.012, 0.01, 0.007, 0.002, 0.009, 0.006, 0.009, 0.013, 0.016, 0.005, 0.001, 0.001, 0.002, 0.005, 0.005, 0.02, 0.03, 0.047, 0.042, 0.046, 0.039, 0.028, 0.023, 0.017, 0.009]


data = np.array([l1])
mean_rates = np.mean(data, axis=0)


