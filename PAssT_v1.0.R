# Load library
library(openxlsx)

# Preparation of MALDI mass features
features <- read.csv("~/MALDI_mass.csv") #Load mass feature list - Add your directory!!
features$m.z <- features$m.z-1.0078 #Substract hydrogen mass

# Load reference database
ref <- read.csv("~/ESI_mass.csv") #Load database - Add your directory!!

#  Loop to match each mass feature to a database hit
dis2feature <- c()
list.mh <- c()
for (i in features){

  # Calculate distance to mass feature
  diff <- abs(ref$Mass - i)

  # Set high pseudo distance if NA
  diff[which(is.na(diff))] <- 999999

  # Remove distance >= 1
  if(min(diff) >= 1.0) next

  # Apply for minimum distance
  dis2feature <- c(dis2feature, min(diff))

  # Extract entries that have minimum distance
  sub <- ref[which(diff == min(diff)),]

  # Add mass feature to extraction
  sub <- cbind(MALDI.mass = rep(i,nrow(sub)), sub)

  # Select maximum logP score
  sub <- sub[which(sub$X.10lgP == max(sub$X.10lgP)),]

  # Add extraction to output
  list.mh <- as.data.frame(rbind(list.mh, sub))
}

# Remove duplicates
list.mh <- list.mh[!duplicated(list.mh),]

# Add distance
dist <- round(abs(list.mh[,1]-list.mh[,4]), digits = 4)
list.mh <- as.data.frame(cbind(Distance = dist, list.mh))

# Count Accession
acc <- unlist(lapply(as.character(list.mh$Accession), function(x) unlist(strsplit(x,":"))))

# Save data frame as xlsx table
wb <-createWorkbook()
addWorksheet(wb, "MH + .calc.")
addWorksheet(wb, "MH + .calc. - peptide")
writeData(wb, "MH + .calc.", list.mh, rowNames = F)

writeData(wb, "MH + .calc. - peptide", table(acc), rowNames = F)
saveWorkbook(wb, "C:/Downloads/", overwrite = T) #Add your directory!!