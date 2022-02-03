# PAssT
Peptide (mass) Assignment Tool (PAssT) allow to identify mass features from MALDI-MSI. A peptide reference library correctly assign MALDI-MSI mass features to peptide signatures and corresponding proteins. Peptide signatures with the lowest mass difference and highest -logP value compared to the reference list were assumed to be a match.

## Getting Started: Quick-Run Protocol
Here, we introduce a *step-by-step* introduction to our **PAssT**, allowing to follow the script.

### Data Upload
First, load the peptide reference library and MALDI-MSI mass feature list.

```{csv files}
Maldi_mass <- read.csv("C:/~/Maldi_mass.csv")
ESI_mass <- read.csv("C:/~/ESI_mass.csv")
```
**Note** that the output of the SCiLS Lab csv file represent <span style="color:red">**[M+H<sup>+</sup>]<sup>+</sup>**</span> mass values and therefore have to substract the hydrogen adduct ion!

```{hydrogen adduct}
ESI_mass$m.z <- ESI_mass$m.z-1.0078
```
### Mass Feature Loop
After data upload, run the mass feature loop:

1. Construct empty vectors with minimal distances for all mass features (*dis2feature*) and respectives masses with minimal distances (*list.mh*)
2. Loop for calculation of distances for all mass features against each others and setting some parameters:
    
    a) Set a high pseudo distance if NA
    b) Remove distance of values >= 1
    c) Select mass features with highest -10logP score
    
```{loop mass assignment}
dis2feature <- c()
list.mh <- c()

for (i in features){

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
```

## How to contribute
