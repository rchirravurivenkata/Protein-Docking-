Molecular docking is a computer simulation procedure to predict the conformation of a receptor-ligand complex. Each docking program makes use of one or more specific search algorithms, which are the methods used to predict the possible conformations of a binary complex. However, in practice with current computational resources, it is impossible to exhaustively explore the search space—this would involve enumerating all possible distortions of each molecule (molecules are dynamic and exist in an ensemble of conformational states) and all possible rotational and translational orientations of the ligand relative to the protein at a given level of granularity. 

This algorithm attempts to find the best (lowest energy) docking pose of a given ligand for a particular molecule/protein by searching in only one angle. Thus, this algorithm is no way near to more complex and comprehensive algorithms. However, this program is still time-intensive. One should expect the program to output results in approximately 20-30 minutes. 


rotation (function)

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

cal_energy (function)

If cal_energy function is used alone to calculate the docking energy, the arguments required are the same as rotation function except optional argument “iter”.


Docking algorithm expected results:

Please run the R script. It saves the rotation function and its dependent cal_energy function. Input the Lennard-Jones parameters, ligand and protein matrices. If needed, please specify the number of iterations required to rotate the ligand from 0 to 360 degrees. The rotation function concurrently performs rotations from 0 to 2pi and calculates energies of the rotated matrices. At the end it prints the rotation matrix associated with the lowest energy. It also writes the iterated degrees and their corresponding energies to a .csv file in the directory. 



References:

- Dias, R., de Azevedo, J., & Walter, F. (2008). Molecular docking algorithms. Current Drug Targets, 9(12), 1040-1047.
- https://en.wikipedia.org/wiki/Docking_(molecular)#Search_algorithm
