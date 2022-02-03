# PAssT
Peptide (mass) Assignment Tool (PAssT) allow to identify mass features from MALDI-MSI. A peptide reference library correctly assign MALDI-MSI mass features to peptide signatures and corresponding proteins. Peptide signatures with the lowest mass difference and highest -logP value compared to the reference list were assumed to be a match.

## Getting Started: Quick-Run Protocol
Here, we introduce a *step-by-step* introduction to our **PAssT", allowing to follow the script.

### Data Upload
First, load the peptide reference library and MALDI-MSI mass feature list.

```{csv files}
Maldi_mass <- read.csv("C:/~/Maldi_mass.csv")
ESI_mass <- read.csv("C:/~/ESI_mass.csv")
```
Please keep in mind to substract a **hydrogen adduct ion** from your MALDI-MSI mass feature list.

```{hydrogen adduct}
ESI_mass$m.z <- ESI_mass$m.z-1.0078
```


## How to contribute
