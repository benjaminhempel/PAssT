# PAssT
Peptide (mass) Assignment Tool (PAssT) allow to identify mass features from MALDI-MSI. A peptide reference library correctly assign MALDI-MSI mass features to peptide signatures and corresponding proteins. Peptide signatures with the lowest mass difference and highest -logP value compared to the reference list were assumed to be a match.

## Getting Started: Quick-Run Protocol
Here, we introduce a *step-by-step* introduction to our **PAssT**, allowing to follow the script.

### Installation
First, packages have to be loaded for further processing (Installation will performed automatically when packages are used for the first time).

```{library}
library(dplyr)
library(openxlsx)
```

### Data Upload
Next, load the peptide reference library (ESI_Mass) and MSI mass feature list (Maldi_mass).

```{csv files}
MALDI_mass <- read.csv("C:/~/MALDI_mass.csv", sep = ";")
ESI_mass <- read.csv("C:/~/ESI_mass.csv", sep = ";")
```
**Note** Output of the SCiLS Lab csv file represent <span style="color:red">**[M+H]<sup>+</sup>**</span> mass values and therefore have to substract the hydrogen adduct ion!

```{hydrogen adduct}
ESI_mass$m.z <- ESI_mass$m.z-1.0078
```
### Mass Feature Loop
After data upload, run the mass feature loop:

1. Construct empty vectors with minimal distances for all mass features (*dis2feature*) and respectives masses with minimal distances (*list.mh*)
2. Loop for calculation of distances for all mass features against each others and setting some parameters:
    
    a) Set a high pseudo distance if NA<br />
    b) Remove distance of values >= 1<br />
    c) Select mass features with highest -10logP score
    
```{mass feature loop}
dis2Maldi_mass <- c()
list.m <- c()

for (i in MALDI_mass$m.z){

  diff <- abs(ESI_mass$Mass - i)
  diff[which(is.na(diff))] <- 999999
  
  if(min(diff) >= 1.0) next
  dis2Maldi_mass <- c(dis2Maldi_mass, min(diff))
  
  sub <- ESI_mass[which(diff == min(diff)),]
  sub <- cbind(MALDI_mass = rep(i,nrow(sub)), sub)
  sub <- sub[which(sub$X.10lgP == max(sub$X.10lgP)),]
  
  list.m <- as.data.frame(rbind(list.m, sub))
}
```

### Data Processing
Further, remove of duplicates, add distances, and count entry accessions in the data frame (*list.mh*)

```{data processing}
dist <- round(abs(list.m$MALDI_mass-list.m$Mass), digits = 4)
list.m <- as.data.frame(cbind(Distance = dist, list.m))

list.m <- list.m[order(list.m[,'Mass'],list.m[,"Distance"]),]
list.m <- list.m[!duplicated(list.m$Mass),]

acc <- unlist(lapply(as.character(list.m$Accession), function(x) unlist(strsplit(x,":"))))
```

### Data Download
Finally, save **PAssT** output as an excel work book in the folder **MassAssignment** on the desktop.

```{save workbook}
wb <-createWorkbook(creator = "BFH" ,title = "PeptideAssignmentTool")
addWorksheet(wb, "Full_Peptide_Library")
writeData(wb, "Full_Peptide_Library", list.m, rowNames = F)
saveWorkbook(wb, "~/PAssT_LIbrary.xlsx", overwrite = T)
```

## How to Contribute
To allow an uncomplicated and quick identification of MALDI mass features
PAssT has been developed to allow for a qualitative and fast identification of MALDI-MSI mass features and is now publicily available for the scientific community. Your help is very valuable to make it better for everyone. The most effective way in regards to bugs or contributions is to open an issue on the github issue tracker. Github allows you to classify your issues into bug report, feature request or feedback to the authors.

   + Check out to see what can be improved.
   + Check out to add new useful features.
   + Open issue for any problems.

## Citation
BF. Hempel, M. Damm, D. Petras, TD. Kazandjian, NR. Casewell, CA. Szentiks, G. Fritsch, G. Nebrich, NR. Casewell, O. Klein, RD. Süssmuth  “Spatial Venomics – Cobra Venom System Reveals Spatial Differentiation of Snake Toxins by Mass Spectrometry Imaging” biorxiv 2022, https://doi.org/10.1101/2022.01.31.478453.

## Update

### PAssT-v1.0.1
   + Small bug fixes to ensure correct match process
   + Update README.md file
