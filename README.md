**About vdwp:**

The vdwp package is writen in python which is designed for efficient prediction of the relative stability of a protein due to a single amino acid point mutation.In this approach, we calculate the free energy change due to an arbitrary point mutation of a protein from a single MD trajectory of the wild type protein. The method is tested on 27 diverse protein systems with a total of 853 mutations and the calculated relative free energies show a generally good correlation with the experimental values (a correlation coefficient of 0.63).
Mutation of samples used the fixbb program of Rosetta. Minimization of method is based on Amber14SB forced filed.

**Code structure:**

The code is organized as follows:
vdwid.py: main program for calculating Van der Waals interaction energy of wild type and mutation samples.<br>
change_filt_type.py: changeing PDB file type between Amber and Rosetta.<br>
tools.py: tools for read Amber's parmtop and crd files.<br>

inputfile: Necessary input files

examples: example amber PDB parm file, crd file, mutation file

**Reference:**

If you use this code in any future publications, please cite this using Wang, B., et al. (2020). "A method for efficient calculation of thermal stability of proteins upon point mutations." Physical Chemistry Chemical Physics 22(16): 8461-8466.

**Use vdwp:**

Necessary third-party software:
Amber (amber16,18 were tested)
Rosetta2017
Python3

calculating protein stability:
python vdwid.py parm.top md.crd frames muta.txt
parm.top: Parm file of Amber
md.crd: Coordinate file of Amber
frames: how many frames you want used for calculating energy
muta.txt: mutation information file
