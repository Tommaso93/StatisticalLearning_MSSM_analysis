# Statistical Learning Ph.D. Course - Tommaso Diotalevi

This repository contains the code written and the plots produced for the TMVA multivariate analysis.

**NB**: In order to reproduce the results, the ROOT framework with TMVA must be installed: https://root.cern.ch/downloading-root contains the binaries for the latest version (v6.18/04) and the instruction for the installation.

## List of the files cointained in this repository

* NewAnalyzer.C : This script contains the analyzer required to read the ROOT files coming from the Monte Carlo simulations and extract only the events that fulfill certain conditions, as shown in the "Event selection" section;
* NewAnalyzer.h : Header file of the NewAnalyser.C script, containing the paths for the MC dataset (available only in a CERN environment);
* analysis_least1bjet.root : ROOT file produced by the NewAnalyzer.C script, containing the "good" events after the selection (containing at least 1 bjet) and divided into signal and background events.

To open a ROOT file, simply:
```
root -l <filename.ROOT>
[1] new TBrowser() 
```
In the **TMVA folder**:

* dataset/ : folder containing the weights and plots produced by the classifier and required for the production of the significance;
* TMVAClassification.C : this script contains the TMVA code for the training and testing of the various BDT algorithms adopted for this project, the output with all the information is then saved in a ROOT file;
* TMVA.root : ROOT file containing all the plots and information concerning the classification. When a new plot is produced, it is automatically saved into the dataset folder.

To open a TMVA ROOT file, simply:

```
root -l 
[1] TMVA::TMVAGui("TMVA.root")
```

* ApplicationTMVA.C : this script applies the classifier trained in the previous file and then - using the BDT response that maximizes the significance - is able to produce the significance in the Higgs search region of 200-400 GeV.
