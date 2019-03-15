rotation function
For rotation of ligand matrices and computation of docking energy between the ligand and protein.
Note: This function requires cal_energy function to calculate energies and cal_energy function does not require additional arguments.
Arguments:
 Mandatory:
•	coeff <-  Input format - Matrix
  Matrix of Lennard-Jones parameters.
•	ligand - Matrix of ligand pdbm format
•	protein  - Matrix of protein pdbm format
Optional:
•	iter – (By default = 100) Number of iterations required to loop through 0 to 360 degrees.

cal_energy function
If cal_energy function is used alone to calculate the docking energy, the arguments required are the same as rotation function except optional argument “iter”. 
Docking algorithm Description:
Please run the R script. It saves the rotation function and its dependent cal_energy function. Input the Lennard-Jones parameters, ligand and protein matrices. If needed, please specify the number of iterations required to rotate the ligand from 0 to 360 degrees. The rotation function concurrently performs rotations from 0 to 2pi and calculates energies of the rotated matrices. At the end it prints the rotation matrix associated with the lowest energy. It also writes the iterated degrees and their corresponding energies to a .csv file in the directory. 


