# evaluateCutpoints

evaluateCutpoints is a R package for multiple cutpoints determination in biomedical research. It is an extension of evaluateCutpoints application (http://wnbikp.umed.lodz.pl/Evaluate-Cutpoints/).
It allows to calculate optimal cutoff value for each biomarker in a dataset and produces a set of statistics for each cutpoint (histograms, Kaplan-Meier curves and estimates - hazard ratios, confidence intervals, p-values).
Additionally, the user can pick a cutpoint value for each biomarker manually.

Publication is available at: www.sciencedirect.com/science/article/pii/S0169260718312252

# Installation

Please use following command to install the package.

```
devtools::install_github("mogluszka/evaluateCutpoints")
```

# Data

Data should contain columns with the results for each biomarker, time and outcome variables.

You can download data example from `data` folder. Then use the following:

```
table <- read.csv("~/pathToYourFile/sample-data.csv", header = TRUE)
df <- data.frame(table)
```

# Package usage example

To test the package, you can use the code below. The algorithm will create a new folder in the main folder and produce statistics for each biomarker.

Required variables:

- **mainDir:** main folder,
- **resultsDirName:** name of the folder where the analysis results will be produced,
- **df:** data frame object,
- **biomarkerList:** a list of analyzed biomarkers,
- **time:** time variable,
- **event:** event variable,
- **setCutpoint:** should be NULL if cutpoints need be calculated; if the user wants to set one cutoff value and produce statistics for two resulting groups, setCutpoint should be a list of cutpoint values for each biomarker in biomarkerList,
- **setCutpoint2:** should be NULL if cutpoints need be calculated; if the user wants to set two cutoff values per biomarker and produce statistics for three resulting groups, setCutpoint2 should be a second list of cutpoint values (additionally to setCutpoint variable) for each biomarker in biomarkerList.

```
library(evaluateCutpoints)

table <- read.csv("~/Documents/ev-cutpoints-multiple-analyses/sample-data.csv", header = TRUE)
df <- data.frame(table)

time <- c("time")
event <- c("event")
biomarkers <- c("biomarker1", "biomarker2", "biomarker3")
setCutpoint <- c(3, 4, 2)
setCutpoint2 <- c(9, 14, 8)

evaluateCutpoints(
  mainDir="~/Documents/folder-name/",
  resultsDirName="analysis-results",
  df=df,
  biomarkerList=biomarkers,
  time=time,
  event=event,
  setCutpoint=NULL,
  setCutpoint2=NULL
)

```
