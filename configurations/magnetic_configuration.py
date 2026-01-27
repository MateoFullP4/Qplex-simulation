"""
This file is designed to create the permanent magnet setup used in Qplex.

This configuration can be imported using the following packages : 
    - sys 
    - os
    - pathlib

To use this configuration in an other folder, you can use this syntax : 

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from configurations.magnetic_configuration import magnetSet_asymmetric

"""

import magpylib as magpy
import numpy as np
from atomsmltr.environment.fields.magnetic.magpylib import MagpylibWrapper


# --- Defininition of all the cubes ---
cube1 = magpy.magnet.Cuboid(magnetization=(-8.7*100000,0,0), dimension=(0.004,0.012,0.020))
cube2 = magpy.magnet.Cuboid(magnetization=(-8.7*100000,0,0), dimension=(0.004,0.012,0.020))
cube3 = magpy.magnet.Cuboid(magnetization=(-8.7*100000,0,0), dimension=(0.004,0.012,0.020))
cube4 = magpy.magnet.Cuboid(magnetization=(-8.7*100000,0,0), dimension=(0.004,0.012,0.012))

cube5 = magpy.magnet.Cuboid(magnetization=(-8.7*100000,0,0), dimension=(0.004,0.012,0.020))
cube6 = magpy.magnet.Cuboid(magnetization=(-8.7*100000,0,0), dimension=(0.004,0.012,0.020))
cube7 = magpy.magnet.Cuboid(magnetization=(-8.7*100000,0,0), dimension=(0.004,0.012,0.020))
cube8 = magpy.magnet.Cuboid(magnetization=(-8.7*100000,0,0), dimension=(0.004,0.012,0.012))

cube9 = magpy.magnet.Cuboid(magnetization=(8.7*100000,0,0), dimension=(0.004,0.012,0.020))
cube10 = magpy.magnet.Cuboid(magnetization=(8.7*100000,0,0), dimension=(0.004,0.012,0.020))
cube11 = magpy.magnet.Cuboid(magnetization=(8.7*100000,0,0), dimension=(0.004,0.012,0.020))
cube12 = magpy.magnet.Cuboid(magnetization=(8.7*100000,0,0), dimension=(0.004,0.012,0.012))

cube13 = magpy.magnet.Cuboid(magnetization=(8.7*100000,0,0), dimension=(0.004,0.012,0.020))
cube14 = magpy.magnet.Cuboid(magnetization=(8.7*100000,0,0), dimension=(0.004,0.012,0.020))
cube15 = magpy.magnet.Cuboid(magnetization=(8.7*100000,0,0), dimension=(0.004,0.012,0.020))
cube16 = magpy.magnet.Cuboid(magnetization=(8.7*100000,0,0), dimension=(0.004,0.012,0.012))

cube17 = magpy.magnet.Cuboid(magnetization=(0,0,8.7*100000), dimension=(0.020,0.012,0.004))
cube18 = magpy.magnet.Cuboid(magnetization=(0,0,-8.7*100000), dimension=(0.020,0.012,0.004))

cube19 = magpy.magnet.Cuboid(magnetization=(8.7*100000,0,0), dimension=(0.004,0.012,0.020))
cube20 = magpy.magnet.Cuboid(magnetization=(8.7*100000,0,0), dimension=(0.004,0.012,0.020))
cube21 = magpy.magnet.Cuboid(magnetization=(8.7*100000,0,0), dimension=(0.004,0.012,0.020))
cube22 = magpy.magnet.Cuboid(magnetization=(8.7*100000,0,0), dimension=(0.004,0.012,0.020))

cube23 = magpy.magnet.Cuboid(magnetization=(-8.7*100000,0,0), dimension=(0.004,0.012,0.020))
cube24 = magpy.magnet.Cuboid(magnetization=(-8.7*100000,0,0), dimension=(0.004,0.012,0.020))
cube25 = magpy.magnet.Cuboid(magnetization=(-8.7*100000,0,0), dimension=(0.004,0.012,0.020))
cube26 = magpy.magnet.Cuboid(magnetization=(-8.7*100000,0,0), dimension=(0.004,0.012,0.020))

cube27 = magpy.magnet.Cuboid(magnetization=(0,0,8.7*100000), dimension=(0.012,0.020,0.004))
cube28 = magpy.magnet.Cuboid(magnetization=(0,0,-8.7*100000), dimension=(0.012,0.020,0.004))
cube29 = magpy.magnet.Cuboid(magnetization=(0,0,8.7*100000), dimension=(0.012,0.020,0.004))
cube30 = magpy.magnet.Cuboid(magnetization=(0,0,-8.7*100000), dimension=(0.012,0.020,0.004))

cube31 = magpy.magnet.Cuboid(magnetization=(0,0,8.7*100000), dimension=(0.012,0.020,0.004))
cube32 = magpy.magnet.Cuboid(magnetization=(0,0,-8.7*100000), dimension=(0.012,0.020,0.004))
cube33 = magpy.magnet.Cuboid(magnetization=(0,0,8.7*100000), dimension=(0.012,0.020,0.004))
cube34 = magpy.magnet.Cuboid(magnetization=(0,0,-8.7*100000), dimension=(0.012,0.020,0.004))


# --- Assigning cubes position ---
cube1.position = (0.034,0,-0.0475)
cube2.position = (0.034,0,-0.0875)
cube3.position = (0.034,0,-0.0675)
cube4.position = (0.030,0,-0.0675)

cube5.position = (-0.034,0,-0.0475)
cube6.position = (-0.034,0,-0.0875)
cube7.position = (-0.034,0,-0.0675)
cube8.position = (-0.030,0,-0.0675)

cube9.position = (0.034,0,0.0475)
cube10.position = (0.034,0,0.0875)
cube11.position = (0.034,0,0.0675)
cube12.position = (0.030,0,0.0675)

cube13.position = (-0.034,0,0.0475)
cube14.position = (-0.034,0,0.0875)
cube15.position = (-0.034,0,0.0675)
cube16.position = (-0.030,0,0.0675)

cube17.position = (-0.04,0,-0.000)
cube18.position = (0.04,0,-0.000)

cube19.position = (0.038,0,0.0675)
cube20.position = (-0.038,0,0.0675)
cube21.position = (0.042,0,0.0675)
cube22.position = (-0.042,0,0.0675)

cube23.position = (0.046,0.012,0.0975)
cube24.position = (-0.046,0.012,0.0975)
cube25.position = (0.046,-0.012,0.0975)
cube26.position = (-0.046,-0.012,0.0975)

cube27.position = (-0.09,0.055,-0.01)
cube28.position = (0.09,0.055,-0.01)
cube29.position = (-0.09,-0.055,-0.01)
cube30.position = (0.09,-0.055,-0.01)

cube31.position = (-0.09,0.08,0.001)
cube32.position = (0.09,0.08,0.001)
cube33.position = (-0.09,-0.08,0.001)
cube34.position = (0.09,-0.08,0.001)


# --- Add all the magnets to the field ---
magnetSet_asymmetric =	magpy.Collection(cube1 + cube2 + cube3 + cube4 + cube5 + cube6 + cube7 + cube8 + cube9 + cube10 + cube11 + cube12 + cube13 + cube14 + cube15 + cube16 + cube17 + cube18 + cube19 + cube20 + cube21 + cube22 + cube23 + cube24 + cube25 + cube26 + cube27 + cube28 + cube29 + cube30 + cube31 + cube32 + cube33 + cube34)
mag_field_asymmetric = MagpylibWrapper(magnetSet_asymmetric)
mag_field_asymmetric.tag = "Real Permanent Magnet Configuration"


# --- MAIN ---
if __name__ == "__main__":
    magnetSet_asymmetric.show()