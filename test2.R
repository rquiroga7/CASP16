# Define the range of ligands
ligands <- paste0("L", 1001:1016)

# Define the number of models
num_models <- 5

# Define the gas constant and temperature
R <- 1.987
T <- 298

# Ensure the output directory exists
output_dir <- "L1000_models"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Loop through each ligand
for (ligand in ligands) {
  # Print header
  header <- paste(
    "PFRMAT     LG",
    "TARGET     L1000",
    "AUTHOR     3419-7309-1619",
    "METHOD     Beta version of 2Vinardo, uses empirical docking scoring function F36 to predict pose and ligand binging affinity",
    sep = "\n"
  )

  ###VER#
    #Now get ligand coordinates for corresponding model from C:/Users/Usuario/Documents/GitHub/CASP16/L1000_results/DOCK.F36.P2hvx/L1001_ligand.sdf
ligand_sdf <- readLines(paste0("C:/Users/Usuario/Documents/GitHub/CASP16/L1000_results/DOCK.F36.P2hvx/", ligand, "_ligand.sdf"))
# Split ligands by lines that contain "$$$$"
ligand_coords <- unlist(strsplit(paste(ligand_sdf, collapse = "\n"), "\\$\\$\\$\\$"))
#######    
  
  # Get protein pdb file coordinates
  protein_pdb <- readLines(paste0("C:/Users/Usuario/Documents/GitHub/CASP16/L1000_results/Prot/2hvx_protein.pdb"))
  
  # Get ligand ID and name
  ligand_info <- read.table(paste0("C:/Users/Usuario/Documents/GitHub/CASP16/L1000_ligands/", ligand, ".tsv"), header = FALSE)
  ligand_id <- ligand_info[1, 1]
  ligand_name <- ligand_info[1, 2]
  
  # Get delta_G for each model
  delta_G_lines <- grep("REMARK 921   NORMALIZED FREE ENERGY", readLines(paste0("C:/Users/Usuario/Documents/GitHub/CASP16/L1000_results/DOCK.F36.P2hvx/", ligand, "_ligand.pdbt")), value = TRUE)
  delta_G_values <- as.numeric(sapply(delta_G_lines, function(line) strsplit(line, " ")[[1]][10]))
  
  # Loop through each model
  for (model in 1:num_models) {
    # Calculate Kd in nM
    delta_G <- delta_G_values[model]
    Kd_nM <- exp(delta_G * 1000 / (R * T)) * 1e9
    
    # Print model information
    model_info <- paste(
      paste("MODEL", model),
      "PARENT 2hvx",
      paste("LIGAND", ligand_id, ligand_name),
      paste("AFFNTY", format(Kd_nM, scientific = FALSE)),
      sep = "\n"
    )


    
    # Combine all parts
    output <- c(header, protein_pdb, model_info)
    
    # Define the output file path
    output_file <- paste0(output_dir, "/", ligand, "LG363_", model, ".txt")
    
    # Write to file
    writeLines(output, output_file)
    
    # Debugging print statement
    cat("Written to file:", output_file, "\n")
  }
}

# Create tar
# tar("L1000LG363.tgz", files = list.files(output_dir, full.names = TRUE))