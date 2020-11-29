#' Calculate cutpoints with Youden method. 
#' 
#' This function calculates the best cutpoint using Youden method, generates a ROC Plot and histogram.
#' @param df input data frame
#' @param time Time variable
#' @param event Event variable
#' @param biomarker Biomarker variable
#' @return ROC plot, youden Plot, a csv result table

youdenMethod <- function(df, time, event, biomarker) {
  
  method <- "Youden"
  
  optimal.cutpoint.Youden <- optimal.cutpoints(X = biomarker, status = event, tag.healthy = 0, methods = method, data = df, pop.prev = NULL, control = control.cutpoints(), ci.fit = FALSE, conf.level = 0.95, trace = FALSE)
  Youdenptable <- summary(optimal.cutpoint.Youden)$p.table$Global
  
  summaryDf <- data.frame(Youdenptable[[1]])
  write.csv(summaryDf, file="YoudenTableOutput.csv")
  
  cutpoint <- Youdenptable[[method]][[1]][[1]]
  
  histogramPlotYouden <- ggplot(df, aes(x=df[, biomarker])) +
    geom_histogram(fill="#2c3e50") +
    geom_vline(aes(xintercept=cutpoint))
  
  histogramYoudenWithLabs <- histogramPlotYouden + labs(x=paste(biomarker))
  ggsave(paste(biomarker, "histogramYouden.png"))
  
  png(paste(biomarker, "ROCYouden.png"))
  youdenPlot <- plot(optimal.cutpoint.Youden, which = 1)
  dev.off()
  
}

#' Calculate cutpoints with ROC method. 
#' 
#' This function calculates the best cutpoint using ROC method, generates a ROC Plot and histogram.
#' @param df input data frame
#' @param time Time variable
#' @param event Event variable
#' @param biomarker Biomarker variable
#' @return ROC plot, youden Plot, a csv result table

ROCMethod <- function(df, time, event, biomarker) {
  
  method <- "ROC01"
  
  optimal.cutpoint.roc <- optimal.cutpoints(X = biomarker, status = event, tag.healthy = 0, methods = method, data = df, pop.prev = NULL, control = control.cutpoints(), ci.fit = FALSE, conf.level = 0.95, trace = FALSE)
  Rocptable <- summary(optimal.cutpoint.roc)$p.table$Global
  
  summaryDf <- data.frame(Rocptable[[1]])
  write.csv(summaryDf, file="ROCTableOutput.csv")
  
  
  cutpoint <- Rocptable[[method]][[1]][[1]]
  
  histogramPlotROC <- ggplot(df, aes(x=df[, biomarker])) +
    geom_histogram(fill="#2c3e50") +
    geom_vline(aes(xintercept=cutpoint))
  
  histogramROCWithLabs <- histogramPlotROC + labs(x=paste(biomarker))
  ggsave(paste(biomarker, "histogramROC.png"))
  
  png(paste(biomarker, "ROCPlot.png"))
  rocPlot <- plot(optimal.cutpoint.roc, which = 1)
  dev.off()
  
}

biomarkerCounts <- function(biomarker, cutpoint) {
  length.total <- length(biomarker)
  length.lower <- length(biomarker[biomarker < cutpoint ])
  length.upper <- length(biomarker[biomarker > cutpoint ])
  length.lower.ratio <- length.lower/length.total
  length.upper.ratio <- length.upper/length.total
  length.lower.result <- paste(length.lower, "(", percent(length.lower.ratio), ")")
  length.upper.result <- paste(length.upper, "(", percent(length.upper.ratio), ")")
  length.results <- c(length.total, length.lower.result, length.upper.result)
  return(length.results)
}

biomarkerCountsThree <- function(biomarker, low.cutoff, high.cutoff) {
  length.total <- length(biomarker)
  print(low.cutoff)
  print(length.total)
  length.lower <- length(biomarker[biomarker < low.cutoff])
  length.upper <- length(biomarker[biomarker > high.cutoff])
  length.medium <- length.total - length.lower - length.upper
  length.lower.ratio <- length.lower / length.total
  length.upper.ratio <- length.upper / length.total
  length.medium.ratio <- length.medium / length.total
  length.lower.result <- paste(length.lower, "(", percent(length.lower.ratio), ")")
  length.upper.result <- paste(length.upper, "(", percent(length.upper.ratio), ")")
  length.medium.result <- paste(length.medium, "(", percent(length.medium.ratio), ")")
  length.results <- c(length.total, length.lower.result, length.medium.result, length.upper.result)
  return(length.results)
}

resultTableSurvival <- function(biomarker, cutpoint, survival, event, df, x) {
  
  length.results <- biomarkerCounts(biomarker, cutpoint)
  
  res.cox <- do.call(
    coxph,
    list(formula = Surv(survival, event) ~ x, data = df)
  )
  
  model <- summary( res.cox )
  coef <- model$coefficients
  HR <- signif( coef[ 2 ], digits = 3 )
  q <- 1-(1-95/100)/2
  z <- qnorm( q )
  HR.lower <- signif( exp( coef[1] - z * coef[3] ), digits = 3 )
  HR.upper <- signif( exp( coef[1] + z * coef[3] ), digits = 3 )
  CI <- paste(" ( ", HR.lower, "-", HR.upper, " ) " )
  p <- signif( model$sctest[ "pvalue" ], digits = 3 )
  
  result.table.col.names <- c( "Cutpoint", "Biomarker < Cutpoint", "Biomarker > Cutpoint", "HR", "CI", "P-value" )
  result.table.row.names <- c( cutpoint, as.character(length.results[2]), as.character(length.results[3]), HR, CI, p )
  result.table <- data.frame(result.table.col.names, result.table.row.names)
  colnames(result.table) <- c("Estimate", "Result")
  return(result.table)
}

#' Calculate cutpoints with cutp method. 
#' 
#' This function calculates the best cutpoint for a continuous variable with coxph and survfit model.
#' @param df input data frame
#' @param time Time variable
#' @param event Event variable
#' @param biomarker Biomarker variable
#' @return Csv result table, histogram, Kaplan-Meier plot, heatmap

cutPMethod <- function(df, time, event, biomarker) {
  
  method <- "cutp"
  
  vector.biomarker <- df[, biomarker ]
  vector.survival <- df[, time ]
  vector.event <- df[, event ]
  
  cph1 <- do.call(
    coxph,
    list( formula = Surv( vector.survival, vector.event ) ~ vector.biomarker, data = df)
  )
  
  
  allCutpointsDf <- data.frame(cutp(cph1))
  colnames(allCutpointsDf) <- c(biomarker, "U", "Q", "pvalue")
  write.csv(allCutpointsDf, file="allCutpointsOutput.csv")
  biomarkerPValueDf <- data.frame( allCutpointsDf[biomarker], allCutpointsDf$pvalue )
  colnames(biomarkerPValueDf) <- c("biomarker", "pvalue")
  biomarkerPValueDfSorted <- biomarkerPValueDf[order(biomarkerPValueDf$biomarker),]
  colnames(biomarkerPValueDfSorted) <- c(biomarker, "pvalue")
  write.csv(biomarkerPValueDfSorted, file="outputSortedbyBiomarkerValue.csv")
  
  named.biomarkerPValueDfSorted <- data.frame(biomarkerPValueDfSorted)
  colnames(named.biomarkerPValueDfSorted) <- c("biomarker", "pvalue")
  bind.cutpointsDf <- cbind(named.biomarkerPValueDfSorted$pvalue, named.biomarkerPValueDfSorted$pvalue)
  rownames(bind.cutpointsDf) <- named.biomarkerPValueDfSorted$biomarker
  rev.cutpointsDf <- t(bind.cutpointsDf)
  colnames(rev.cutpointsDf) <- rownames(bind.cutpointsDf)
  rownames(rev.cutpointsDf) <- c( ""," ")
  
  y.axisSettings <- list(
    title = "",
    zeroline = FALSE,
    showline = FALSE,
    showticklabels = FALSE,
    showgrid = FALSE,
    ticks = ""
  )
  x.axisSettings <- list(
    title = biomarker
  )
  
  heatmap <- plot_ly(
    z = data.matrix(rev.cutpointsDf),
    x = colnames(rev.cutpointsDf),
    y = rownames(rev.cutpointsDf),
    type = "heatmap",
    colorbar = list(
      xanchor = "left",
      yanchor = "top",
      ypad = 0,
      xpad = 0,
      lenmode = "pixels",
      len= 150,
      nticks=3
    )) %>% layout(
      yaxis = y.axisSettings,
      xaxis = x.axisSettings
    )
  
  htmlwidgets::saveWidget(heatmap, file = "heatmap.html")
  
  cutpoint <- signif(allCutpointsDf[[ 1, 1 ]], digits = 4 )
  p.value <- signif(allCutpointsDf[[ 1, 4 ]], digits = 4 )
  
  category <- ifelse( vector.biomarker < cutpoint, "low", "high" )
  x <- ifelse( vector.biomarker < cutpoint, 0, 1 )
  category.df <- data.frame( df, category, x )
  
  result.table.cutp <- resultTableSurvival(vector.biomarker, cutpoint, vector.survival, vector.event, category.df, x)
  write.csv(result.table.cutp, file="CutPResults.csv")
  
  fit <- do.call(
    survfit,
    list( formula = Surv( vector.survival, vector.event ) ~ category, data = category.df )
  )
  
  histogramPlotCutp <- ggplot(df, aes(x=vector.biomarker)) +
    geom_histogram(fill="#2c3e50") +
    geom_vline(aes(xintercept=cutpoint))
  
  histogramCutpWithLabs <- histogramPlotCutp + labs(x=paste(biomarker))
  ggsave(paste(biomarker, "histogramCutP.png"))
  
  ggsurvfit <- ggsurvplot(fit, surv.col = c( "#2c3e50", "Red" ))
  ggsave(paste(biomarker, "KaplanMeierPlotCutp.png"))
}

#' Calculate cutpoints with maxstat method. 
#' 
#' This function calculates the best cutpoint for a continuous variable using maximally selected rank statistics.
#' @param df input data frame
#' @param time Time variable
#' @param event Event variable
#' @param biomarker Biomarker variable
#' @return Csv result table, histogram, Kaplan-Meier plot, Standarized log-rank statistics

maxstatMethod <- function(df, time, event, biomarker) {
  
  method <- "maxstat"
  
  vector.biomarker <- df[, biomarker ]
  vector.survival <- df[, time ]
  vector.event <- df[, event ]
  
  mod <- maxstat.test( 
    Surv( vector.survival, vector.event ) ~ vector.biomarker, 
    data=df, smethod="LogRank", pmethod="Lau92", iscores=TRUE
  )
  
  cutpoint <- signif(mod$estimate[[ 1 ]], digits = 4)
  p.value <- mod$p.value
  
  modstats <- mod$stats
  modcuts <-mod$cuts
  dataxy <- data.frame(modcuts, modstats)
  
  maxstatplot <- ggplot( dataxy, aes( x=modcuts, y=modstats ) ) +
    geom_line(colour= "#2c3e50") +
    geom_point(colour= "#2c3e50") +
    geom_vline(aes(xintercept=cutpoint)) +
    labs(x=biomarker,y="Standarized log-rank statistics")
  ggsave(paste(biomarker, "MaxstatPlot.png"))
  
  category <- ifelse( vector.biomarker < cutpoint, "low", "high" )
  x <- ifelse(vector.biomarker < cutpoint, 0, 1)
  category.df <- data.frame( df, category, x )
  
  fit <- do.call(
    survfit,
    list( formula = Surv( vector.survival, vector.event ) ~ category, data = category.df )
  )
  
  result.table.cutp <- resultTableSurvival(vector.biomarker, cutpoint, vector.survival, vector.event, category.df, x)
  write.csv(result.table.cutp, file="MaxStatResults.csv")
  
  histogramPlotMaxstat <- ggplot(df, aes(x=vector.biomarker)) +
    geom_histogram(fill="#2c3e50") +
    geom_vline(aes(xintercept=cutpoint))
  
  histogramMaxstatWithLabs <- histogramPlotMaxstat + labs(x=paste(biomarker))
  ggsave(paste(biomarker, "histogramMaxstat.png"))
  
  ggsurvfit <- ggsurvplot(fit, surv.col = c( "#2c3e50", "Red" ))
  ggsave(paste(biomarker, "KaplanMeierPlotMaxstat.png"))
  
}

#' Calculate cutpoints with maxstat method. 
#' 
#' This function calculates two cutpoints for a continuous variable using hierarchical method.
#' @param df input data frame
#' @param time Time variable
#' @param event Event variable
#' @param biomarker Biomarker variable
#' @param setCutpoint cutpoint variable
#' @return Csv result table, histogram, Kaplan-Meier plots

rolrMethod <- function(setCutpoint, df, time, event, biomarker) {
  
  vector.biomarker <- df[, biomarker ]
  vector.survival <- df[, time ]
  vector.event <- df[, event ]
  
  func <- function(df) {
    df <- df[, c(biomarker, time, event)]
    df <- na.omit(df)
    res <- rhier(times = df[, time ], status = df[, event ], x = df[, biomarker ], ns = 15, alt = 'increase')
  }
  
  rolr <- func( df )
  
  if (is.null(setCutpoint)) {
    low.cutoff = rolr[[ 1 ]][[ 1 ]]
    high.cutoff = rolr[[ 1 ]][[ 2 ]]
  } else {
    low.cutoff = setCutpoint[1]
    high.cutoff = setCutpoint[2]
  }
  
  
  result.table.col.names <- c("Cutpoint 1", "Cutpoint 2", "Low", "Medium", "High")
  length.results <- biomarkerCountsThree(vector.biomarker, low.cutoff, high.cutoff)
  result.table.row.names <- c(low.cutoff, high.cutoff, as.character(length.results[2]), as.character(length.results[3]), as.character(length.results[4]))
  result.table.rolr <- data.frame(result.table.col.names, result.table.row.names)
  colnames(result.table.rolr) <- c("Estimate", "Result")
  write.csv(result.table.rolr, file="GroupCounts.csv")
  
  category <- ifelse(vector.biomarker < low.cutoff, "low", ifelse(vector.biomarker > high.cutoff, "high", "medium" ))
  category.df <- data.frame(vector.biomarker, vector.survival, vector.event, category )
  
  category.df.lowmed <- data.frame( subset( category.df, subset= category == "medium" | category == "low" ))
  category.df.lowhigh <- data.frame( subset( category.df, subset= category == "high" | category == "low" ))
  category.df.medhigh <- data.frame( subset( category.df, subset= category == "high" | category == "medium" ))
  
  
  fit <- do.call(
    survfit,
    list( formula = Surv( vector.survival, vector.event ) ~ category, data = category.df )
  )
  
  if (low.cutoff != high.cutoff) {
    
    lowhigh <- do.call(
      survfit,
      list( formula = Surv( vector.survival, vector.event ) ~ category, data = category.df.lowhigh)
    )
    
    lowmedium <- do.call(
      survfit,
      list( formula = Surv( vector.survival, vector.event ) ~ category, data = category.df.lowmed)
    )
    
    mediumhigh <- do.call(
      survfit,
      list( formula = Surv( vector.survival, vector.event ) ~ category, data = category.df.medhigh)
    )
    
    cat.df <- data.frame( vector.biomarker, vector.survival, vector.event )
    
    categorize.low <- subset( cat.df, subset = vector.biomarker < low.cutoff | biomarker > high.cutoff )
    category.low = ifelse( categorize.low$vector.biomarker < low.cutoff, "low", "high" )
    cat.low <- data.frame( categorize.low, category.low )
    
    lowhigh.surv <- do.call(
      coxph,
      list( formula = Surv( vector.survival, vector.event ) ~ category.low, data = cat.low)
    )
    
    categorize.med <- subset( cat.df, subset = vector.biomarker < high.cutoff )
    category.med = ifelse( categorize.med$vector.biomarker < low.cutoff, "medium", "high" )
    cat.med <- data.frame(categorize.med, category.med)
    
    lowmedium.surv <- do.call(
      coxph,
      list( formula = Surv( vector.survival, vector.event ) ~ category.med, data = cat.med)
    )
    
    categorize.high <- subset( cat.df, subset = vector.biomarker > low.cutoff )
    category.high = ifelse( categorize.high$vector.biomarker < high.cutoff, "medium", "high" )
    cat.high <- data.frame( categorize.high, category.high )
    
    mediumhigh.surv <- do.call(
      coxph,
      list( formula = Surv( vector.survival, vector.event ) ~ category.high, data = cat.high)
    )
    
    groupEstimates <- function (group, summaryModel) {
      model <- summary( summaryModel )
      coef <- model$coefficients
      HR <- signif(coef[ 2 ], digits = 3 )
      q <- 1-(1-95/100)/2
      z <- qnorm(q)
      HR.lower <- signif( exp( coef[ 1 ] - z * coef[ 3 ]), digits = 3 )
      HR.upper <- signif( exp( coef[ 1 ] + z * coef[ 3 ]), digits = 3 )
      CI <- paste( " (", HR.lower, "-", HR.upper, ") " )
      p <- signif( model$sctest[ "pvalue" ], digits = 3 )
      estimates <- c(group,HR,CI, p)
      return(estimates)
    }
    
    estimates.total.table <- rbind(
      groupEstimates("Low vs High", lowhigh.surv),
      groupEstimates("Low vs Medium", lowmedium.surv),
      groupEstimates("Medium vs High", mediumhigh.surv)
    )
    colnames(estimates.total.table) <- c("Group", "HR", "CI", "p-value")
    write.csv(estimates.total.table, file="EstimatesTable.csv")
    
    lowmediumfit.plot <- ggsurvplot(lowmedium, surv.col = c( "Red", "Blue" ))
    ggsave(paste(biomarker, "LowVsMedium.png"))
    lowhighfit.plot <- ggsurvplot(lowhigh, surv.col = c( "#2c3e50", "Red" ))
    ggsave(paste(biomarker, "LowVsHigh.png"))
    mediumhighfit.plot <- ggsurvplot(lowhigh, surv.col = c("#2c3e50", "Red"))
    ggsave(paste(biomarker, "MediumVsHigh.png"))
    lowmediumhighfit.plot <- ggsurvplot(fit, surv.col = c("#2c3e50", "Red", "Blue"))
    ggsave(paste(biomarker, "LowVsMediumVsHigh.png"))
    
    histogramPlotRolr <- ggplot(df, aes(x=vector.biomarker)) +
      geom_histogram(fill="#2c3e50") +
      geom_vline(aes(xintercept=low.cutoff)) +
      geom_vline(aes(xintercept=high.cutoff))
    
    histogramRolrWithLabs <- histogramPlotRolr + labs(x=paste(biomarker))
    ggsave(paste(biomarker, "histogramRolr.png"))
    
    
  } else {
    print("Please select two distinct cutpoints.")
  }
}

#' Calculate statistics for one selected cutpoint. 
#' 
#' This function calculates statistics if one set of cutpoints is selected.
#' @param df input data frame
#' @param time Time variable
#' @param event Event variable
#' @param biomarker Biomarker variable
#' @param setCutpoint cutpoint variable
#' @return Csv result table, histogram, Kaplan-Meier, Standarized log-rank statistics plots

oneCutpointAdapt <- function(setCutpoint, df, time, event, biomarker) {
  
  vector.biomarker <- df[, biomarker ]
  vector.survival <- df[, time ]
  vector.event <- df[, event ]
  
  cph1 <- do.call(
    coxph,
    list( formula = Surv( vector.survival, vector.event ) ~ vector.biomarker, data = df)
  )
  
  
  allCutpointsDf <- data.frame(cutp(cph1))
  colnames(allCutpointsDf) <- c(biomarker, "U", "Q", "pvalue")
  
  category <- ifelse( vector.biomarker < setCutpoint, "low", "high" )
  x <- ifelse( vector.biomarker < setCutpoint, 0, 1 )
  category.df <- data.frame( df, category, x )
  
  result.table.cutp <- resultTableSurvival(vector.biomarker, setCutpoint, vector.survival, vector.event, category.df, x)
  write.csv(result.table.cutp, file="CutPResults.csv")
  
  fit <- do.call(
    survfit,
    list( formula = Surv( vector.survival, vector.event ) ~ category, data = category.df )
  )
  
  histogramPlotCutp <- ggplot(df, aes(x=vector.biomarker)) +
    geom_histogram(fill="#2c3e50") +
    geom_vline(aes(xintercept=setCutpoint))
  
  histogramCutpWithLabs <- histogramPlotCutp + labs(x=paste(biomarker))
  ggsave(paste(biomarker, "histogramCutP.png"))
  
  ggsurvfit <- ggsurvplot(fit, surv.col = c( "#2c3e50", "Red" ))
  ggsave(paste(biomarker, "KaplanMeierPlotCutp.png"))
  
  mod <- maxstat.test( 
    Surv( vector.survival, vector.event ) ~ vector.biomarker, 
    data=df, smethod="LogRank", pmethod="Lau92", iscores=TRUE
  )
  
  modstats <- mod$stats
  modcuts <-mod$cuts
  dataxy <- data.frame(modcuts, modstats)
  
  maxstatplot <- ggplot( dataxy, aes( x=modcuts, y=modstats ) ) +
    geom_line(colour= "#2c3e50") +
    geom_point(colour= "#2c3e50") +
    geom_vline(aes(xintercept=setCutpoint)) +
    labs(x=biomarker,y="Standarized log-rank statistics")
  ggsave(paste(biomarker, "MaxstatPlot.png"))
  
}