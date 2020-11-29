#' A function used when both sets of cutpoints are NULL
#' This function calculates the best cutpoint using cutp, maxstat, ROC and Youden methods, generates ROC Plots, Kaplan-Meier plots and histograms for 2 and 3 groups.
#' @param df input data frame
#' @param time Time variable
#' @param event Event variable
#' @param biomarker Biomarker variable
#' @return 2 folders for each biomarker with results

defaultAnalyses <- function(df, time, event, biomarker) {
  numberOfGroups = c("2","3")
  setwd(biomarker)
  for (n in numberOfGroups) {
    dir.create(n)
    setwd(n)
    if (n == 2) {
      methods = c("Youden", "ROC01", "cutp", "maxstat")
      for (m in methods) {
        dir.create(m)
        setwd(m)
        if (m == "Youden") {
          youdenMethod(df, time, event, biomarker)
        } else if (m == "ROC01") {
          ROCMethod(df, time, event, biomarker)
        } else if (m == "cutp") {
          cutPMethod(df, time, event, biomarker)
        } else if (m == "maxstat") {
          maxstatMethod(df, time, event, biomarker)
        }
        setwd('..')
      }
    }
    if (n == 3) {
      dir.create("rolr")
      setwd("rolr")
      rolrMethod(setCutpoint, df, time, event, biomarker)
      setwd('..')
    }
    setwd('..')
  }
  setwd('..')
}

#' A function used when one set of cutpoints is defined
#' This function generates ROC Plots, Kaplan-Meier plots and histograms for a selected cutpoint value.
#' @param df input data frame
#' @param time Time variable
#' @param event Event variable
#' @param biomarker Biomarker variable
#' @param setCutpoint cutpoint variable
#' @return a folders for each biomarker with results

adaptCutpoint2groups <- function(setCutpoint, df, time, event, biomarker) {
  setwd(biomarker)
  dirname <- paste(as.character(setCutpoint), "_2groups")
  dir.create(dirname)
  setwd(dirname)
  tryCatch({
    oneCutpointAdapt(setCutpoint, df, time, event, biomarker)
  }, error = function(e) {
    print("Please select cutpoints between minimum and maximum biomarker value")
  })
  setwd('../..')
}

#' A function used when two sets of cutpoints are defined
#' This function generates ROC Plots, Kaplan-Meier plots and histograms for two selected cutpoint values.
#' @param df input data frame
#' @param time Time variable
#' @param event Event variable
#' @param biomarker Biomarker variable
#' @param setCutpoint cutpoint variable
#' @return a folders for each biomarker with results

adaptCutpoint3groups <- function(setCutpoint, df, time, event, biomarker) {
  setwd(biomarker)
  dirname <- paste(as.character(setCutpoint[1]), as.character(setCutpoint[2]), "_3groups")
  dir.create(dirname)
  setwd(dirname)
  
  tryCatch({
    rolrMethod(setCutpoint, df, time, event, biomarker)
  }, error = function(e) {
    print("Please select cutpoints between minimum and maximum biomarker value")
  })
  
  setwd('../..')
}

#' A function that decides which type of analysis to use.
#' 
#' 
#' @param df input data frame
#' @param time Time variable
#' @param event Event variable
#' @param biomarker Biomarker variable
#' @param cutpoints cutpoints variable
#' @return results from a selected analysis

mainFunction <- function(cutpoints, df, time, event, biomarker) {
  if (is.null(cutpoints)) {
    defaultAnalyses(df, time, event, biomarker)
  } else if (length(cutpoints) == 1) {
    adaptCutpoint2groups(cutpoints, df, time, event, biomarker)
  } else if (length(cutpoints) == 2) {
    adaptCutpoint3groups(cutpoints, df, time, event, biomarker)
  } else {
    print("Please select a correct cutpoint value")
  }
}

#' A function creating folders for each biomarker in a data table with subfolders that include results of the analysis.
#' @param df input data frame
#' @param time Time variable
#' @param event Event variable
#' @param biomarker Biomarker variable
#' @param cutpoints cutpoints variable
#' @return results from a selected analysis
#' @export

evaluateCutpoints <- function (mainDir, file, resultsDirName, biomarkerList, time, event, setCutpoint, setCutpoint2, df) {
  
  table <- read.table(file, sep = "\t", header = TRUE)
  df <- data.frame(table)
  
  setwd(mainDir)
  
  if ( dir.exists(resultsDirName)) {
    print("Directory already exists")
  } else {
    dir.create(resultsDirName)
  }
  
  currentDate <- format(Sys.time(), "%a %b %d %Y %X")
  setwd(resultsDirName)
  dir.create(currentDate)
  
  colNames <- colnames(df)
  colsCommaSeparated <- dput(colNames)
  
  setwd(currentDate)
  
  prepareCutpointsDf <- function(biomarker, biomarkerList, cutpoint = NULL, cutpoint2 = NULL) {
    
    if (is.null(cutpoint) && is.null(cutpoint2)) {
      return (NULL)
    } else {
      cutpointsFunc(biomarker, biomarkerList, cutpoint, cutpoint2)
    }
  }
  
  cutpointsFunc <- function(biomarker, biomarkerList, cutpoint = NULL, cutpoint2 = NULL) {
    
    if (is.null(cutpoint2)) {
      cutpointDf <- data.frame(rbind(cutpoint))
    } else {
      cutpointDf <- data.frame(rbind(cutpoint, cutpoint2))
    }
    colnames(cutpointDf) <- biomarkerList
    rownames(cutpointDf) <- NULL
    
    result <- cutpointDf[[biomarker]]
    
  }
  
  if(is.null(setCutpoint) || is.null(setCutpoint2) || (length(setCutpoint) = length(biomarkerList)) || (length(setCutpoint2) = length(biomarkerList))) {
    for (b in biomarkerList) {
      dir.create(b)
      timeLen <- length(time)
      eventLen <- length(event)
      if ((timeLen == 0) || (eventLen == 0)) {
        print("Please correct time or event variables")
      } else if ((timeLen == 1) && (eventLen == 1)) {
        print("One time and event variable - results will be produced in biomarker folder")
        cutpointResult <- R.utils::doCall(prepareCutpointsDf, args=list(biomarker = b, biomarkerList=biomarkerList, cutpoint = setCutpoint, cutpoint2 = setCutpoint2))
        mainFunction(cutpoints = cutpointResult, df, time, event, b)
      } else if ((timeLen > 1) | (eventLen > 1)) {
        print("Please provide one time and one event variable")
      }
    }
  } else {
    print("Please select correct biomarkers, time, event and cutpoint values")
  }
}
