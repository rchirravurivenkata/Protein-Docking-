###############################################################################################################################################
###### Docking Algorithm ######################################################################################################################
##### By Ramakanth Chirravuri Venkata #########################################################################################################
###############################################################################################################################################

## Function for Calculating Energies
## Input the lookup table of LJ terms, ligand and protein matrices

cal_energy <- function(coeff, ligand, protein){
  energies <- c()
  for(i in 1:nrow(ligand)){
    for(j in 1:nrow(protein)){

      ## Calculates Euclidean distance between the protein and the ligand
      
      euclidean <- (sqrt(((protein[j, 7]-ligand[i, 7])**2) + ((protein[j, 8]-ligand[i, 8])**2)+ ((protein[j, 9]-ligand[i, 9])**2)))/10
      idx1 <- which( coeff[, 1] == protein[j, 3]) 
      idx2 <- which(coeff[, 2] == ligand[i, 3])
      idx <- intersect(idx1, idx2)
      if ((length(idx)) < 1){
      idx1 <- which( coeff[, 1] == ligand[i, 3]) 
      idx2 <- which( coeff[, 2] == protein[j, 3]) 
       idx <- intersect(idx1, idx2)}
      
      
      ## Calculates Lennard-Jones term 
      
      v_lj <- (coeff[idx, 4]/(euclidean**12)) - (coeff[idx, 3]/(euclidean**6))
      
      ## Calculates Electrostatics term
      
      if (euclidean < 1.32){
        epsilon = 8}
      else{
      epsilon <- 6.02944 + (72.37056/(1+(213.5782* (2.718281**(-0.018733345*72.37056*euclidean)))))}
      
    
      v_e <- (138.935485 * ((ligand[i, 10] * protein[j, 10])/(epsilon* euclidean)))
      vij <- v_lj + v_e
      
      energies <- append(vij, energies)
    }}
  return(sum(energies))}

###########################################################################################################################

#### Rotation function to simulataneously perform rotations and energy calculations

rotation <- function(coeff, ligand, protein, iter = 100){
  energy_vals <- cal_energy(coeff, ligand, protein)
  positions <- c(0)
  coord <- as.matrix(ligand[, c(7,8,9)])
  coord[, 1] <- as.numeric(as.character(coord[, 1]))
  coord[, 2] <- as.numeric(as.character(coord[, 2]))
  coord[, 3] <- as.numeric(as.character(coord[, 3]))
  
  
  ## Loop to make 100 iterations from 0 to 2pi and simultaneously calculate the energy 
  
  for(i in seq(1:iter)){
    print(i)
    
    degree <- as.numeric((6.2831853072/iter) * i)    
    mat <- matrix(0, 3, 3)
    diag(mat) <- 1
    sin_value <- sin(degree)
    cos_value <- cos(degree)
    mat[1,1] <- cos_value
    mat[3,3] <- cos_value
    mat[1, 3] <- sin_value
    mat[3,1] <- -(sin_value)
    interim <- coord %*% mat
    new_ligand <- ligand
    new_ligand [, 7:9] <- interim
    set_energy_vals <- cal_energy(coeff, new_ligand, protein)
    if(as.numeric(set_energy_vals) < min(as.numeric(energy_vals))){
      final_matrix <- mat}
    energy_vals <- append(energy_vals, set_energy_vals)
    positions <- append(positions, (as.numeric(degree) * 57.295779513))}
  
  
  ## Index position of lowest energy
  
  pos_low_energy <- which(energy_vals == min(energy_vals, na.rm = T))
  print(noquote("***********************************************************************************************"))
  print(noquote("                                       Results                                                  "))
  print(noquote("***********************************************************************************************"))
  print(noquote("                                                                                                "))
  print(paste("Lowest energy value is:", min(energy_vals, na.rm = T), sep = " "))
  print(paste("Degree associated with lowest energy value is:", positions[pos_low_energy], sep = " "))
  deg_n_energy <- data.frame(positions, energy_vals)
  colnames(deg_n_energy) <- c("Degrees", "Energy Values")
  write.csv(deg_n_energy, "degrees_and_energy.csv", quote = F, row.names = F)
  print(noquote("                                                                                                "))
  print(noquote("                                                                                                "))
  print(noquote("The rotation matrix that produces best docking pose is returned"))
  print(noquote("                                                                                                "))
  return(final_matrix)
}
        
##############################################################################################################################################################     

### Root Mean Square Deviation Function


## Provide structure numeric coordinate matrices of same size as x and y
## If only few rows of the matrices have to considered, provide number of rows (n).

rmsd <- function(x,y, n = 0){
  if(n == 0){
    n <- nrow(x)}
  sd <- sqrt(1/n) 
  ligand <- as.matrix(x[1:n, 7:9])
  protein <- as.matrix(y[1:n, 7:9])
  ligand <- apply(ligand, 2, as.numeric)
  protein <- apply(protein, 2, as.numeric)
  a <- ligand - protein
  a <- as.matrix(a)
  for(i in 1:length(a)){
    a[i] <- a[i] **2
  }
  c <- sqrt(sum(a))
  ls <- sd * c
  return(ls)}

########################################################################################################################################