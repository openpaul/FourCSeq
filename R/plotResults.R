##' plotScatter
##' 
##' \code{plotScatter} Plot scatter plots for all experiments with replicates.
##'  
##' @param object A \code{FourC} object.
##' @param assay Character vector selecting the assay of the \code{FourC} object
##'        for which the correlation between the replicates should be plotted.
##' @param ... Additional arguments passed to \code{smoothScatter}.       
##' 
##' @seealso \code{\link[graphics]{smoothScatter}}
##' 
##' @examples
##' 
##'                        
##' data(fc, package="FourCSeq")
##' 
##' fc <- combineFragEnds(fc)
##' fc
##' 
##' plotScatter(fc)
##' 
##' @export
plotScatter <- function(object,
                        assay="counts",
                        ...){    
  
  stopifnot(class(object)=="FourC")
  
  if(!assay %in% names(assays(object))) 
    stop("Select a valid assay of the FourC object.")
  
  ## get the viewpoint data
  vpData = as.data.frame(colData(object))
  
  ## get the count data of the fragment ends
  data <- assays(object)[[assay]]
  
  for (vp in unique(vpData$viewpoint)){
    experiments <- paste(vp, 
                         vpData$condition[vp == vpData$viewpoint], 
                         vpData$replicate[vp == vpData$viewpoint], 
                         sep="_")
    if(length(experiments) < 2) next
    comb <- combinations(length(experiments), 2)
    
    for (i in seq_len(nrow(comb))){
      cols = c(experiments[comb[i,1]], experiments[comb[i,2]])  
      if(!exists("xlab") & !exists("ylab")){
        plotSmoothScatter(vp, cols,
                          data[,cols], 
                          xlab=cols[2],
                          ylab=cols[1],
                          ...)
      } else {
        plotSmoothScatter(vp, cols,
                          data[,cols], ...) 
      }
    }
  }
}


## plot the correlation for any given viewpoint and condition
plotSmoothScatter <- function(viewpoint, rep, ftc, ...){
  filter = apply(ftc, 1, sum) == 0

  heatscatter(log10(ftc[!filter, 2] + 0.1), 
              log10(ftc[!filter, 1] + 0.1), 
              xaxt = "n", yaxt = "n",
              colpal=blues9[-c(1:2)],
              cor=FALSE,
              ...)
  x = c(0, 10^0, 10^1, 10^2, 10^3, 10^4, 10^5, 10^6, 10^7) + 0.1
  labels = c(0, 1, 10, 100, expression(10^3), expression(10^4),
             expression(10^5), expression(10^6), expression(10^7))
  axis(1, at=log10(x),labels=labels, las=2, cex.axis=.8)
  axis(2, at=log10(x),labels=labels, las=2, cex.axis=.8)
}

##' Plot differences
##' 
##' \code{plotDifferences} generate plots to investigate the results 
##' of \code{getDifferences}
##' 
##'                
##' @param object A \code{FourC} object.
##' @param plotWindows Window sizes around the viewpoint for which plots are generated.
##' @param textsize Adjust text size.
##' @param diffThresh Threshold on adjusted p-values calculated in the differential test.
##' @param controls Auxiliary \code{GRanges} object that contains information about the 
##' viewpoint and closest gene TSS.
##' @param txdb Auxiliary \code{TranscriptDb} object that contains gene models for the 
##' investigated organism.
##' @param conditionAB Condition A and B, for which the comparison is plotted.
##' 
##' @return NULL
##' @author Felix A. Klein, \email{felix.klein@@embl.de}
##' @seealso \code{\link{FourC}}, \code{\link{getDifferences}}
##' @examples
##' 
##'                        
##' data(fcf, package="FourCSeq")
##' 
##' fcf <- getDifferences(fcf, referenceCondition="WE_68h")
##' 
##' plotDifferences(fcf)
##' 
##' @export
plotDifferences <- function(object,
                            plotWindows = c(1e+05, 1e+06),
                            textsize = 20,
                            diffThresh = 0.01,
                            controls=NULL,
                            txdb=NULL,
                            conditionAB=NULL){    
  
  stopifnot(class(object)=="FourC")
  
  if(!"results" %in% mcols(mcols(object))$type)
    stop("No results available, call getDifferences first.")  
  if( !(is.null(conditionAB) | length(conditionAB) == 2) )
      stop("conditionAB has to contain exactly two conditions.")
  
  ## set theme for ggbio
  theme_new = theme_bw()
  theme_new$text$size = textsize
  theme_new$axis.text$size = rel(0.8)
  theme_set(theme_new)
  
  ## get the viewpoint data
  vpData = colData(object)
  
  lvls <- levels(vpData$condition)
  
  if(length(lvls) > 1){
  
    ## get fragment data
    fragData <- rowRanges(object)
    
    ####################################
    ## go through all viewpoints
    ####################################
    viewpoint = unique(vpData$viewpoint)
    print(viewpoint)
    
    dse <- object[mcols(object)$selectedForFit,]
    
    trafo <- assay(dse, "trafo")
    fit <- assay(dse, "fit")
    zScore <- assay(dse, "zScore")
    pValue <- assay(dse, "pValue")
    fdr <- assay(dse, "pAdjusted")
    sd <- colData(dse)$sd
    frags <- rowRanges(dse)

    if("peaks" %in% names(assays(dse))){
      peaks = assay(dse, "peaks")
    } else {
      dse <- addPeaks(dse)
      peaks = assay(dse, "peaks")
    }
    
    zScoreThresh = exptData(dse)$peakParameter$zScoreThresh
        
    if(!is.null(txdb)){
      exons = reduce(exons(txdb, columns="gene_id"))
      cds = reduce(cds(txdb, columns="gene_id"))
      transcripts = reduce(transcripts(txdb, columns="gene_id"))
      
      which = frags[abs(frags$dist) < max(plotWindows),]
      
      genes = ggplot() + 
        geom_rect(subsetByOverlaps(exons, which), rect.height = 25) + 
        geom_rect(subsetByOverlaps(transcripts, which), rect.height = 2) + 
        geom_rect(subsetByOverlaps(cds, which), rect.height = 50)
    }
    
    ## go through all condition pairs
    if(is.null(conditionAB)) conditionAB <- lvls[1:2]
    vpConditions <- paste(viewpoint, paste(conditionAB, collapse="_"), sep="_")
    
    deseqResults = results(dse[mcols(dse)$selectedForFit,], 
                           contrast=c("condition", conditionAB))
    
    differential <- rep(FALSE, nrow(deseqResults))
    differential[deseqResults$padj < diffThresh] = TRUE
    
    vpFragMid = unique(colData(dse)$start + colData(dse)$width %/% 2)
    
    colsA <-  rownames(vpData[which(vpData$condition == conditionAB[1]),])
    colsB <-  rownames(vpData[which(vpData$condition == conditionAB[2]),])
    cols = c(colsA, colsB)
                              
    ## plot first condition
    plotsA <- lapply(colsA, function(col, distMax=NULL){
      tmp <- as.data.frame(rowRanges(dse))
      tmp <- tmp[,c(1:9, match("dist",names(tmp)))]
      tmp$mid = tmp$start + (tmp$end - tmp$start)%/%2 
      tmp$count <- trafo[,col]
      tmp$fit <- fit[,col]
      tmp$fitUp <- tmp$fit + zScoreThresh*sd[col]
      tmp$fitDown <- tmp$fit - zScoreThresh*sd[col]    
      tmp$peak <- peaks[,col]
      tmp$differentialInteraction = differential
      
      if(is.null(distMax)) distMax = max(abs(tmp$dist)) + 1 
      plotVp(tmp[abs(tmp$dist) < distMax & tmp$dist!=0,], 
             col, 
             time="")}, distMax=max(plotWindows))
    names(plotsA) <- colsA
    
    ## plot second condition
    plotsB <- lapply(colsB, function(col, distMax=NULL){
      tmp <- as.data.frame(rowRanges(dse))
      tmp <- tmp[,c(1:9, match("dist",names(tmp)))]
      tmp$mid = tmp$start + (tmp$end - tmp$start)%/%2 
      tmp$count <- trafo[,col]
      tmp$fit <- fit[,col]
      tmp$fitUp <- tmp$fit + zScoreThresh*sd[col]
      tmp$fitDown <- tmp$fit - zScoreThresh*sd[col]    
      tmp$peak <- peaks[,col]
      tmp$differentialInteraction = differential
      
      if(is.null(distMax)) distMax = max(abs(tmp$dist)) + 1 
      plotVp(tmp[abs(tmp$dist) < distMax & tmp$dist!=0,], 
             col, 
             time="")}, distMax=max(plotWindows))
    names(plotsB) <- colsB
    
    forPlot <- as.data.frame(rowRanges(dse))
    forPlot <- forPlot[,c(1:9, match("dist",names(forPlot)))]
    forPlot$mid = forPlot$start + (forPlot$end - forPlot$start)%/%2 
    forPlot$peak <- apply(peaks, 1, any)
    forPlot$differentialInteraction = differential
    forPlot <- cbind(forPlot, as.data.frame(deseqResults))
    forPlot$logFDR <- pmax(log10(forPlot$padj), -10)
    
    forPlot$change <- with(forPlot, 
                           ifelse(differential, 
                                  ifelse(log2FoldChange > 0, 1, -1),
                                  0))
    
    colorBar <- ggplot(forPlot, aes(xmax=end, xmin=start, ymax=1, ymin=0, fill=change)) + 
      geom_rect() +
      scale_fill_gradient2(low="red", high="green", mid="white", guide=FALSE) + 
      theme(axis.text.y = element_blank(), axis.ticks.y= element_blank())
    
    plots <- c(plotsA, plotsB, list(colorBar))
    heights <- c(rep(3, length(plotsA) + length(plotsB)), 1)
    
    log2FoldChange <- ggplot(forPlot, 
                             aes(mid, log2FoldChange)) + 
      labs(x="genomic position", y="log2 fold change") +
      geom_path(colour="darkgrey", size=0.5, na.rm = TRUE) +
      geom_point(size=3, na.rm = TRUE) +
      scale_y_continuous(limits=c(min(forPlot$log2FoldChange, na.rm=TRUE), 
                                  max(forPlot$log2FoldChange, na.rm=TRUE)))
    
    log2FoldChange <- log2FoldChange +
      geom_point(aes(mid, log2FoldChange), 
                 colour=2, 
                 data=subset(forPlot, peak & !differentialInteraction),
                 size=4,
                 na.rm = TRUE) +
      geom_point(aes(mid, log2FoldChange), 
                 colour=4, 
                 data=subset(forPlot, !peak & differentialInteraction),
                 size=5,
                 shape=18,
                 na.rm = TRUE) +
      geom_point(aes(mid, log2FoldChange), 
                 colour="darkorange", 
                 data=subset(forPlot, peak & differentialInteraction),
                 size=5,
                 shape=18,
                 na.rm = TRUE)
    plots <- c(plots, list(log2FoldChange))    
    heights <- c(heights, 3)
    
    if(exists("genes")){
      plots[["genes"]] <- genes
      heights <- c(heights, 1.2)
    }      
    
    if(!is.null(controls) & viewpoint %in% controls$vp){
      plots[["VP/TSS"]] <- autoplot(controls[controls$vp == viewpoint,])
      heights = c(heights, 1)
    }
    
    for(windowSize in plotWindows){
      print(tracks(plots, 
                   heights=heights,
                   xlim=frags[abs(frags$dist) < windowSize],
                   xlab="Genomic position",
                   label.text.cex = 1.1))
    }
    
    cat("Successfully plotted results.\n")
    
  } else stop("Only one level and no differential results to plot.")
}


##' Plot z-score results.
##' 
##' \code{plotZScores} generates plots to check the z-score calculation.
##' 
##' 
##' @param object A \code{FourC} object.
##' @param cols A character or numeric vector used to subset the columns of \code{object}.
##' @param plotWindows Window sizes around the viewpoint for which plots are generated.
##' @param textsize Plot parameter passed to the plotting function if
##' @param controls Auxiliary \code{GRanges} object that contains information about the 
##' viewpoint and closest gene TSS.
##' @param txdb Auxiliary \code{TranscriptDb} object that contains gene models for the 
##' investigated organism.
##' @param plotSingle If set to true each selected column will be plotted in a single plot.
##' @return NULL
##' 
##' @details
##' Plots are generated to visualize the results of the getZScores function
##' 
##' @author Felix A. Klein, \email{felix.klein@@embl.de}
##' @seealso \code{\link{FourC}}, \code{\link{getZScores}},
##' \code{\link{distFitMonotone}}, \code{\link{distFitMonotoneSymmetric}}
##' @examples
##' 
##'                        
##' data(fcf, package="FourCSeq")
##' 
##' plotZScores(fcf)
##'  
##' @export
plotZScores <- function(object,
                        cols=NULL,
                        plotWindows = c(1e+05, 1e+06),
                        controls=NULL,
                        textsize = 20,
                        txdb=NULL,
                        plotSingle=FALSE){    
  
  stopifnot(class(object)=="FourC")
  
  if( ! c("zScore") %in% names(assays(object)))
    stop("No z-scores available, call getZScores first.")
  
  ## set theme for ggbio
  theme_new = theme_bw()
  theme_new$text$size = textsize
  theme_new$axis.text$size = rel(0.8)
  theme_set(theme_new)
  
  ## get the viewpoint data
  vpData = colData(object)
  
  ## prepare for DESeqSE object
  countData = assays(object)[["counts"]]

  ####################################
  ## go through all viewpoints
  ####################################
  for(viewpoint in unique(vpData$viewpoint)){
    print(viewpoint)
    
    ## select the data columns of the current viewpoint
    if(is.null(cols)) cols = rownames(vpData[which(vpData$viewpoint==viewpoint),])
    
    ## replace rowRanges with fragData containing distance information and convert fragData to data.frame
    fragData = as.data.frame(rowRanges(object))
    
    ## select only fragments that passed the filter and selected columns
    dse <- object[mcols(object)$selectedForFit,cols]
    
    trafo <- assay(dse, "trafo")
    fit <- assay(dse, "fit")
    zScore <- assay(dse, "zScore")
    pValue <- assay(dse, "pValue")
    pAdjusted <- assay(dse, "pAdjusted")
    sd <- colData(dse)$sd
    
    if("peaks" %in% names(assays(dse))){
      peaks = assay(dse, "peaks")
    } else {
      dse <- addPeaks(dse)
      peaks = assay(dse, "peaks")
    }
    
    zScoreThresh = exptData(dse)$peakParameter$zScoreThresh
    
    if(plotSingle){
      for(col in cols){
        tmp <- fragData
        tmp$count <- trafo[,col]
        tmp$fit <- fit[,col]
        tmp$fitUp <- tmp$fit + zScoreThresh*sd[col]
        tmp$fitDown <- tmp$fit - zScoreThresh*sd[col]    
        tmp$peak <- peaks[,col]
        
        for(windowSize in plotWindows){
          print(plotVp(tmp[abs(tmp$dist) < windowSize & tmp$dist!=0,], 
                       col, 
                       time=""))        
        }
      }
    } else {
      frags <- rowRanges(dse)
      
      viewpoint <- gsub("CRM_", "", unique(colData(dse)$viewpoint))
      
      vpFragMid = unique(colData(dse)$start + colData(dse)$width %/% 2)
      
      if(!is.null(txdb)){
        exons = reduce(exons(txdb, columns="gene_id"))
        cds = reduce(cds(txdb, columns="gene_id"))
        transcripts = reduce(transcripts(txdb, columns="gene_id"))
        
        which = frags[abs(frags$dist) < max(plotWindows),]
        
        genes = ggplot() + 
          geom_rect(subsetByOverlaps(exons, which), rect.height = 25) + 
          geom_rect(subsetByOverlaps(transcripts, which), rect.height = 2) + 
          geom_rect(subsetByOverlaps(cds, which), rect.height = 50)
      }
      
      vpData=colData(dse)

      plots <- lapply(cols, function(col, distMax=NULL){
        tmp <- as.data.frame(rowRanges(dse))
        tmp <- tmp[,c(1:9, match("dist",names(tmp)))]
        tmp$mid = tmp$start + (tmp$end - tmp$start)%/%2 
        tmp$count <- trafo[,col]
        tmp$fit <- fit[,col]
        tmp$fitUp <- tmp$fit + zScoreThresh*sd[col]
        tmp$fitDown <- tmp$fit - zScoreThresh*sd[col]
        tmp$peak <- peaks[,col]
        
        if(is.null(distMax)) distMax = max(abs(tmp$dist)) + 1 
        plotVp(tmp[abs(tmp$dist) < distMax & tmp$dist!=0,], 
               col, 
               time="")}, distMax=max(plotWindows))
      names(plots) <- cols
      
      heights <- rep(3, length(plots))
      
      if(exists("genes")){
        plots[["genes"]] <- genes
        heights <- c(heights, 1)
      }
      
      if(!is.null(controls) & viewpoint %in% controls$vp){
        plots[["VP/TSS"]] <- autoplot(controls[controls$vp == viewpoint,])
        heights = c(heights, 1)
      }
      
      for(windowSize in plotWindows){
        print(tracks(plots, 
                     heights=heights,
                     xlim=frags[abs(frags$dist) < windowSize],
                     xlab="Genomic position",
                     label.text.cex = 1.1))
      }
    }
  }
  cat("Successfully plotted results.\n")
}


plotVp <- function(forPlot, col, time=""){
  plot <- ggplot(forPlot, aes(mid, count)) + 
    labs(title = paste(col, time), x="Genomic coordinate", y="Count (vst)") +
    geom_path(colour="darkgrey", size=0.5, na.rm = TRUE) +
    geom_path(aes(mid, fit), colour=3, na.rm = TRUE) +
    geom_path(aes(mid, fitUp), colour=4, linetype=2, na.rm = TRUE) +
    geom_path(aes(mid, fitDown), colour=4, linetype=2, na.rm = TRUE)+
    geom_point(size=3, na.rm = TRUE) +
    scale_y_continuous(limits=c(min(forPlot$count, na.rm=TRUE), 
                                max(forPlot$count, na.rm=TRUE)))
  
  if("differentialInteraction" %in% names(forPlot)){
    plot <- plot +
      geom_point(aes(mid, count), 
                 colour=2, 
                 data=subset(forPlot, peak & !differentialInteraction),
                 size=4,
                 na.rm = TRUE) +
      geom_point(aes(mid, count), 
                 colour=4, 
                 data=subset(forPlot, !peak & differentialInteraction),
                 size=5,
                 shape=18,
                 na.rm = TRUE) +
      geom_point(aes(mid, count), 
                 colour="darkorange", 
                 data=subset(forPlot, peak & differentialInteraction),
                 size=5,
                 shape=18,
                 na.rm = TRUE)
  } else {
    plot <- plot +
      geom_point(aes(mid, count), 
                 colour=2, 
                 data=subset(forPlot, peak),
                 size=4,
                 na.rm = TRUE)
  }
  plot
}


##' Plot fit results.
##' 
##' \code{plotFits} generates plots of the fits used to calculate the z-scores.
##' 
##' @param object A \code{FourC} object, after successfully calling \code{getZScores}.
##' @param viewpoint A character vector of the viewpoint for which the plots are generated. 
##'        If set to \code{NULL} the first viewpoint in the \code{FourC} object is used.
##' @param main Main text for the plots. If set to \code{NULL} the column names are printed.     
##' @return NULL
##' 
##' @details
##' Plots are generated to visualize the results of the fits used to calculate the z-scores.
##' 
##' @author Felix A. Klein, \email{felix.klein@@embl.de}
##' @seealso \code{\link{FourC}}, \code{\link{getZScores}},
##' \code{\link{distFitMonotone}}, \code{\link{distFitMonotoneSymmetric}}
##' @examples
##' 
##'                        
##' data(fcf, package="FourCSeq")
##' 
##' plotFits(fcf)
##'  
##' @export
plotFits <- function(object,viewpoint=NULL, main=NULL){
  stopifnot(class(object) == "FourC")
  if( ! c("fit") %in% names(assays(object)))
    stop("No fit values available, call getZScores first.")
  
  vpData <- colData(object)
  
  if(is.null(viewpoint)) viewpoint <- unique(vpData$viewpoint)[1]
  
  cols = rownames(vpData[which(vpData$viewpoint==viewpoint),])
  
  ## select vp chromosome and columns
  sel <- mcols(object)$selectedForFit
  object <- object[mcols(object)$vpChr,cols]
  
  ## order by absolute Distance
  object <- object[order(abs(mcols(object)$dist)),]
  
  trafo <- assay(object, "trafo")
  fit <- assay(object, "fit")
  countData <- counts(object)
  sd <- colData(object)$sd
  sdFactor <- exptData(object)$peakParameter$zScoreThresh
  if(is.null(sdFactor)) sdFactor <- 2
  
  tickValue <- 0:(ceiling(log10(max(countData))))
  tickPos <- getVST(object)(10^tickValue)
  tickLabel <- c(1, 10, 100, expression(10^3), expression(10^4), 
                 expression(10^5), expression(10^6), 
                 expression(10^7), expression(10^8))[1:length(tickValue)]
  
  xPos = 0:8
  labels = c(1, 10, 100, expression(10^3), expression(10^4),
             expression(10^5), expression(10^6), expression(10^7),
             expression(10^8))
  
  for(col in cols){
    leftSide <- rowRanges(object)$dist < 0 & sel
    def.par <- par(no.readonly = TRUE)
    par(mar=c(5,4,4,5)+.1)
    plot(-log10(-rowRanges(object)$dist[leftSide]), 
         trafo[leftSide,col], 
         main = ifelse(is.null(main), paste0(col, "_left_side"), main),
         pch=20,
         cex=0.5,
         xlab="Distance from viewpoint [nt]",
         ylab="Counts",
         xaxt="n",
         yaxt="n",         
         ylim=c(min(tickPos), max(tickPos)))
    axis(1, at=-xPos,labels=labels, las=2)
    axis(2, 
         at = tickPos, 
         labels = tickLabel, 
         las=2)
    axis(4)
    mtext("Variance stabilized counts",side=4,line=3)
    points(-log10(-rowRanges(object)$dist[leftSide]), 
           fit[leftSide,col], 
           type="l", 
           col=2)
    points(-log10(-rowRanges(object)$dist[leftSide]), 
           fit[leftSide,col] +sdFactor*sd[col], 
           type="l",
           lty=2,
           col=4)
    
    rightSide <- rowRanges(object)$dist > 0 & sel
    plot(log10(rowRanges(object)$dist[rightSide]), 
         trafo[rightSide,col], 
         main = ifelse(is.null(main), paste0(col, "_right_side"), main),
         pch=20,
         cex=0.5,
         xlab="Distance from viewpoint [nt]",
         ylab="Counts",
         xaxt="n",
         yaxt="n",
         ylim=c(min(tickPos), max(tickPos)))
    axis(1, at=xPos,labels=labels, las=2)
    axis(2, 
         at = tickPos, 
         labels = tickLabel, 
         las=2)
    axis(4)
    mtext("Variance stabilized counts",side=4,line=3)
    points(log10(rowRanges(object)$dist[rightSide]), 
           fit[rightSide,col], 
           type="l", 
           col=2)
    points(log10(rowRanges(object)$dist[rightSide]), 
           fit[rightSide,col] +sdFactor*sd[col], 
           type="l",
           lty=2,
           col=4)
    
    
    plot(log10(abs(rowRanges(object)$dist[sel])), 
         trafo[sel,col], 
         main = ifelse(is.null(main), paste0(col, "_combined"), main),
         pch=20,
         cex=0.5,
         xlab="Distance from viewpoint [nt]",
         ylab="Counts",
         xaxt="n", 
         yaxt="n",
         ylim=c(min(tickPos), max(tickPos)))
    axis(1, at=xPos,labels=labels, las=2)
    axis(2, 
         at = tickPos, 
         labels = tickLabel, 
         las=2)
    axis(4)
    mtext("Variance stabilized counts",side=4,line=3)
      points(log10(abs(rowRanges(object)$dist[sel])), 
           fit[sel,col], 
           type="l", 
           col=2)
    points(log10(abs(rowRanges(object)$dist[sel])), 
           fit[sel,col] +sdFactor*sd[col], 
           type="l",
           lty=2,
           col=4)
    par(def.par)
  }
}

##' Plot the estimated normalization factors.
##' 
##' \code{plotNormalizationFactors} generates plots of the fits used to calculate the z-scores.
##' 
##' @param object A \code{FourC} object, after successfully calling \code{getZScores}.
##' @param dist If \code{TRUE} (default), the normalizationFactors are plotted as a 
##'        function for distance from the viewpoint. Otherwise they are plotted
##'        as a function of ranks.
##' @param sizeFactors If \code{TRUE}, sizeFactors are calculated and plotted as
##'        horizontal bars.
##' @return NULL
##' 
##' @details
##' Plots are generated to visualize the results of the fits used to calculate the z-scores.
##' 
##' @author Felix A. Klein, \email{felix.klein@@embl.de}
##' @seealso \code{\link{FourC}}, \code{\link{getDifferences}}, 
##'         \code{\link[DESeq2]{estimateSizeFactors}}
##' @examples
##' 
##'                        
##' data(fcf, package="FourCSeq")
##' 
##' fcf <- getDifferences(fcf)
##' 
##' plotNormalizationFactors(fcf)
##'  
##' @export
plotNormalizationFactors <- function(object, dist=TRUE, sizeFactors=TRUE){
  
  stopifnot(class(object) == "FourC")
  if( ! c("normalizationFactors") %in% names(assays(object)))
    stop("No normalizationFactors available, call getDifferences first.")
  
  normFactors <- normalizationFactors(object)
  xlab = ifelse(is.null(dist), "Index", "Distance from viewpoint")
  if(dist){
    dist = mcols(object)$dist
  } else {
    dist = 1:nrow(normFactors)
  }
  plot(dist,
       normFactors[,1], 
       type="l", 
       ylab="Normalization Factor",
       xlab= xlab,
       log = "y",
       ylim=c(max(min(normFactors)-0.1, 0.001), max(normFactors)+0.1))
  for(i in 2:ncol(normFactors)) points(dist, normFactors[,i], type="l", col=i)
  legend("topright", 
         pch = 20,
         legend=colnames(normFactors), 
         col=1:ncol(normFactors),
         cex=0.5)
  if(sizeFactors){
    sizeFactors <- sizeFactors(estimateSizeFactors(object))
    for(i in 1:ncol(normFactors)) abline(h=sizeFactors[i], col=i)
  }
}
