# Purpose

This code uses symmetry and framework atom positions to speed up adsorption energy calculation on a gemmi grid associated to a nanoporous material.

# Compile and run

c++ -std=c++11 -I./include -O3 src/main.cpp -o graes.out <br>
./graes.out structure/KAXQIL_clean_14.cif forcefield/UFF.def 298.0 12.0 Xe 0.12 100 0.8

# Results format

KAXQIL_clean_14,-44.626,2.21797,2.14037,10.0172,15.0287,73.1689,0.0300725,0.330829
Structure name, Enthalpy (kJ/mol), enthalpy_std, enthalpy_skew, enthalpy_kurt, mean_grid, std_grid, Henry coeff (mol/kg/Pa), Time (s)

# Acknowledgement

This code has been developped during a PhD thesis co-financed by the CEA and Orano under the supervision of François-Xavier Coudert: https://github.com/fxcoudert

This code includes the library developped in Gemmi:
https://github.com/project-gemmi/gemmi.git
