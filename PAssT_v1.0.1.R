# Load library
if (!require(openxlsx)) install.packages('openxlsx')
if (!require(dplyr)) install.packages('dplyr')

# Features
Maldi_mass <- read.csv("~/Maldi_mass.csv", sep = ";")
Maldi_mass$m.z <- Maldi_mass$m.z-1.0078 #convert from dataframe to list and substract hydrogen

# Load reference
ESI_mass <- read.csv("~/ESI_mass.csv", sep = ";")

# Function to compare masses
dis2Maldi_mass <- c()
list.m <- c()
for (i in Maldi_mass$m.z){

  # Calc distance to feature
  diff <- abs(ESI_mass$Mass - i)

  # Set high pseudo distance if na
  diff[which(is.na(diff))] <- 999999

  # Remove distance >= 1
  if(min(diff) >= 1.0) next

  # Add min distance
  dis2Maldi_mass <- c(dis2Maldi_mass, min(diff))

  # Extract entries that have min distance
  sub <- ESI_mass[which(diff == min(diff)),]

  # Add feature to extraction
  sub <- cbind(MALDI_mass = rep(i,nrow(sub)), sub)

  # Select max score
  sub <- sub[which(sub$X.10lgP == max(sub$X.10lgP)),]

  # Add extraction to output
  list.m <- as.data.frame(rbind(list.m, sub))
}

# Add distance
dist <- round(abs(list.m$MALDI_mass-list.m$Mass), digits = 4)
list.m <- as.data.frame(cbind(Distance = dist, list.m))

# Remove duplicates
list.m <- list.m[order(list.m[,'Mass'],list.m[,"Distance"]),]
list.m <- list.m[!duplicated(list.m$Mass),]

# Count Accession
acc <- unlist(lapply(as.character(list.m$Accession), function(x) unlist(strsplit(x,":"))))

# Save
wb <-createWorkbook(creator = "BFH" ,title = "PeptideASSignmentTool")
addWorksheet(wb, "Full_Peptide_Library")
writeData(wb, "Full_Peptide_Library", list.m, rowNames = F)
saveWorkbook(wb, "~/PAssT_LIbrary.xlsx", overwrite = T)