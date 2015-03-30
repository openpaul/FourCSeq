##' Combine the counts of both fragment ends.
##' 
##' 
##' \code{combineFragEnds} combines the counts of both fragment ends. A
##' multiplication factor can be used for fragments that only have counts for
##' one valid fragment end.
##' 
##' @param object A \code{FourC} object.
##' @param multFactor Multiplication factor that can be used to multiply the
##' counts of fragments which only have one valid end. Default is \code{1}.
##' \code{filter} is automatically set to \code{TRUE} if \code{multFactor} 
##' is different from 1.
##' @param filter If filter is \code{TRUE}, only reads from valid fragment
##' ends are summed up for each fragment. This means if only one fragment end
##' is valid only counts from this end are considered.
##' If filter is \code{FALSE}, counts from both fragment ends are summed up
##' without any filtering.
##' @return Returns an updated \code{FourC} object with a new assay \code{counts}
##' containing the combined count data of both fragment ends for all viewpoints.
##' @author Felix A. Klein, \email{felix.klein@@embl.de}
##' @seealso \code{\link{FourC}}, \code{\link{countFragmentOverlaps}}
##' @examples
##' 
##' 
##' 
##'                        
##' data(fc, package="FourCSeq")
##' 
##' fc <- combineFragEnds(fc)
##' fc
##' 
##' @export
combineFragEnds <- function(object,
                            multFactor=1,
                            filter=FALSE){

  stopifnot(class(object)=="FourC")
  
  ## 
  if(any( ! c("countsLeftFragmentEnd","countsRightFragmentEnd") %in% names(assays(object))))
      stop("Use countFragmentOverlaps to import the counts for the left and right fragment ends.")
  
  if(multFactor != 1 | filter){
    ## create empty matrix for the combined counts
    counts=matrix(NA, nrow=length(rowRanges(object)), 
      ncol=nrow(colData(object)))

    ## get fragments for which both ends are valid
    bothValid = mcols(object)$leftValid & mcols(object)$rightValid

    counts[bothValid,] <- assay(object, "countsLeftFragmentEnd")[bothValid,] +
      assay(object, "countsRightFragmentEnd")[bothValid,]
  
    counts[!bothValid & rowRanges(object)$leftValid] <- multFactor *
      assay(object, "countsLeftFragmentEnd")[!bothValid & rowRanges(object)$leftValid,]
  
    counts[!bothValid & rowRanges(object)$rightValid] <- multFactor * 
      assay(object, "countsRightFragmentEnd")[!bothValid & rowRanges(object)$rightValid,]

  } else {
    counts <- assay(object, "countsLeftFragmentEnd") + assay(object, "countsRightFragmentEnd")
  }
  if (any(round(counts) != counts)) {
    stop("Some values in counts are not integers.
          Make sure that countsRight/LeftFragmentEnd and multFactor
          only contain integers.")
  } 
  mode(counts) <- "integer"

  assays(object)[["counts"]] <- counts
  ## reorder that counts is the first assay as required by DESeq2
  assays(object) <- assays(object)[c("counts", setdiff(names(assays(object)), "counts"))]

  object
}


##' Write track files of an selected \code{assay}
##' 
##' The files are saved in the specified folder. The filenames are the 
##' combination of the assay name, the selected column name and the corresponding
##' file extension.
##'
##' \code{writeTrackFiles} 
##' 
##' @param object A \code{FourC} object.
##' @param assay Character vector selecting the assay of the \code{FourC} object
##'         that should be saved as track file.
##' @param folder Path relative to the project folder, where the results are track files should be saved.
##' @param format Character vector specifying the format of the output.
##' Can either be 'bedGraph' or 'bw'. 'bw' is the default.
##' @param removeZeros Define whether fragments with zero counts should be included with
##'        value 0 or not. On default zeros are removed.        
##' @return Message whether the track export of assay was successful.
##' @author Felix A. Klein, \email{felix.klein@@embl.de}
##' @examples  
##' 
##' 
##' 
##'                        
##' data(fc, package="FourCSeq")
##' exptData(fc)$projectPath = tempdir()
##' 
##' fc <- combineFragEnds(fc)
##' fc
##' 
##' writeTrackFiles(fc)
##' 
##' @export
writeTrackFiles <- function(object, 
                          assay="counts", 
                          folder="tracks",
                          format="bw",
                          removeZeros=TRUE){    

  stopifnot(class(object)=="FourC")
  
  if(!format %in% c("bw", "bedGraph"))
     stop("Format has to be 'bw' or 'bedGraph'")  
  
  if(!assay %in% names(assays(object))) 
    stop("Select a valid assay of the FourC object.")

  ## get the fragment intervals
  ivCols <- c("seqnames", "start", "end")
  fragData = rowRanges(object)
  mcols(fragData) = NULL

  ## get the count data of the fragment ends
  data <- assays(object)[[assay]]
  
  outputDir <- file.path(exptData(object)$projectPath, folder)
  if(!file.exists(outputDir)) dir.create(outputDir)
  
  tmp <- options(scipen=12, digits=2)
  
  for(name in colnames(data)){
    outputFile = file.path(outputDir, paste0(assay, "_", name, ".", format))

    notNA <- !is.na(data[, name])
    if(removeZeros) notNA[notNA][(data[notNA,name] == 0)] <- FALSE

    gr = fragData[notNA,]
    gr$score = data[notNA, name]
    gr=keepSeqlevels(gr, unique(as.character(seqnames(gr))))
    export(gr, outputFile, format=format)
    rm(gr)
  }
  options(tmp)
  return(paste0("Successfully created ", format, 
                " files of the ", assay, " data."))
}

##' Normalize count data to rpm
##' 
##' \code{normalizeRPM} Normalizes the counts of each experiment to RPM. 
##' 
##' @param object A \code{FourC} object.
##' @param assay Character vector selecting the assay of the \code{FourC} object
##'        that should be normalized.
##' @param normalized Character vector specifying the assay name of the normalized
##'        data.
##' @param removeHighestValues Defines the number of highest count values that should
##'        be removed before the library sizes are calculated in the denominator.
##' \code{normalized} containing the normalized data.
##' @author Felix A. Klein, \email{felix.klein@@embl.de}
##' @seealso \code{\link{FourC}}, \code{\link{combineFragEnds}}
##' @examples  
##' 
##'                        
##' data(fc, package="FourCSeq")
##' 
##' fc <- combineFragEnds(fc)
##' fc
##' 
##' fc <- normalizeRPM(fc)
##' fc
##' 
##' @export
normalizeRPM <- function(object, 
                         assay="counts",
                         normalized="rpm",
                         removeHighestValues=0){    
  
  stopifnot(class(object)=="FourC")
  
  if(!assay %in% names(assays(object))) 
    stop("Select a valid assay of the FourC object.")
  
  assays(object)[[normalized]] <- apply(assays(object)[[assay]], 
                                        2, 
                                        function(x, remove){
                                          if(remove > 0) 
                                            1e+6*x/sum(x[x < sort(x, decreasing=TRUE)[remove]])
                                          else 
                                            1e+6*x/sum(x)
                                          }, 
                                        remove=removeHighestValues)
  object  
}

##' getDistAroundVp
##' 
##' \code{getDistAroundVp} Add the distance information for a given viewpoint 
##' to the fragment data.
##' 
##' @param vp A single character string with the selected viewpoint.
##' @param vpData \code{DataFrame} containing the column data of the \code{FourC} object.
##' @param fragData \code{GRanges} object containing the row data (fragment data) of 
##' the \code{FourC} object.
##' 
##' @return Updated \code{GRanges} object that contains the distance from the viewpoint 
##'         for the fragments on the viewpoint chromosome. 
##'         
##' @details Mid positions of the fragments are defined as start + (end - start) %/% 2
##' 
##' @examples
##' 
##' 
##' data(fc, package="FourCSeq")
##' 
##' fragmentDataWithDistance <- getDistAroundVp("ap", colData(fc), rowRanges(fc))
##' 
##' @export
getDistAroundVp <- function(vp, vpData, fragData){
  ## calculate the distance from the view point for the same chromosome using the fragment mid points
  ## get the genomic coordinates of the current viewpoint
  curVp = unique(vpData[vpData$viewpoint==vp,c("chr", "start","end")])
  curVp$mid = curVp$start + (curVp$end - curVp$start)%/%2
  
  fragData$mid = start(fragData) + (end(fragData) - start(fragData))%/%2
    
  fragData$dist = NA

  sameChr = as.logical(seqnames(fragData) == curVp$chr)

  fragData$dist[sameChr] = fragData$mid[sameChr] - curVp$mid

  fragData$posLeft = fragData$dist<0
  fragData$posLeft[is.na(fragData$posLeft)] = FALSE
  fragData$posRight = fragData$dist>0
  fragData$posRight[is.na(fragData$posRight)] = FALSE

  ## add information about blind fragments
  fragData$blind = fragData$leftSize < 0  & fragData$rightSize < 0
  
  ## add vp chromosome
  fragData$vpChr = sameChr

  newCols <- c("mid", "dist", "posLeft", "posRight", "blind", "vpChr")
  mcolsRows <- DataFrame(type=rep("viewpointInfo", length(newCols)),
                         description=rep("",length(newCols)))
  mcols(mcols(fragData))[colnames(mcols(fragData)) %in% newCols,] <- mcolsRows
  
  fragData
}

##' Detect differences
##' 
##' \code{getDifferences} detects differences in the interaction frequencies between 
##' different conditions. For each viewpoint all possible combinations of 
##' conditions are tested using the \code{DESeq2} package.
##' 
##'                
##' @param object A \code{FourC} object.
##' @param referenceCondition Reference condition to be used for testing.
##' @param fitNormFactors Defines if the distance dependency should be used to 
##' estimate normalizationFactors. Default is \code{TRUE}.
##' @return Returns an updated \code{FourC} object that contains the results and
##'         information from the differential testing. Normalization factors
##'         are added to the assays slot.
##' @author Felix A. Klein, \email{felix.klein@@embl.de}
##' @seealso \code{\link{FourC}}, \code{\link{getZScores}}, 
##'          \code{\link[DESeq2]{estimateDispersions}},
##'          \code{\link[DESeq2]{nbinomWaldTest}}
##' @examples
##' 
##'                        
##' data(fcf, package="FourCSeq")
##' 
##' fcf <- getDifferences(fcf, referenceCondition="WE_68h")
##' 
##' results <- getAllResults(fcf)
##' results
##' 
##' @export
getDifferences <- function(object,
                           referenceCondition=NULL,
                           fitNormFactors=TRUE){      
  stopifnot(class(object)=="FourC")
  
  if( ! c("counts") %in% names(assays(object)))
    stop("No assay 'counts' found. Use combineFragEnds first.")
  if ("results" %in% mcols(mcols(object))$type) {
    stop("Differences are already calculated. To recalculate use the object returned by 'getZScores'.")
  }
  
  ## get the viewpoint data
  viewpoint = unique(colData(object)$viewpoint)
  if(length(viewpoint) != 1)
    stop("None or more than one viewpoint are contained in the 'FourC' object.
         Use a 'FourC' object that contains only one viewpoint.")
  print(viewpoint)
  
  if(!fitNormFactors) viewpoint = paste0(viewpoint, "_noNormFactors")

  colData(object)$condition <- factor(colData(object)$condition)
  
  if(is.null(referenceCondition)){
    referenceCondition <- levels(colData(object)$condition)[1]
    cat("No reference condition provided. ", referenceCondition, 
        " was chosen by default.")
  }
  
  if(referenceCondition %in% levels(colData(object)$condition)){
    colData(object)$condition <- relevel(colData(object)$condition, referenceCondition)
  } else {
    cat("Provided reference condition ", referenceCondition, " not found in levels.\n")
    cat("Using ", levels(colData(object)$condition)[1], " as reference condition instead.")  
  }
  object <- estimateSizeFactors(object)
  
  ## fit normalizationFactors/estimate dispersions
  if(fitNormFactors) {
    object <- try(getNormalizationFactors(object))
    if(class(object) == "try-error") break     
  }
  
  ## test for differences  
  design(object) <- formula(~ condition)
  object <- try(estimateDispersions(object))
  object <- nbinomWaldTest(object)
  
  invisible(object)
}


##' Calculate z-scores using the residuals of the general trend fit.
##' 
##' \code{getZScores} calculates the z-score for each fragment within a
##' distance from the viewpoint defined by \code{distAroundVP}. For the
##' calculation the count data is first transformed with a variance stabilizing
##' transformation from the \code{DESeq2} package. The decay trend of this transformed
##' data with distance from the viewpoint is fitted either with local regression
##' model using \code{locfit} or a monotone decay fit using the \code{fda} package. 
##' z-scores are finally calculated from the residuals of the fit values.
##' 
##' @param object A \code{FourC} object.
##' @param removeZeros Parameter to define if fragments with zero counts across 
##' all replicates and conditions should be removed from the calculation. 
##' Default is \code{TRUE}.
##' @param minCount The minimum median count value across replicates, a fragment 
##' has to have to be used for the calculation of z-scores. Default is \code{40}
##' @param minDist Define a minimum distance to both sides of the viewpoint that is
##' discarded for analysis. If undefined (\code{NULL}), the borders of the 
##' viewpoint peak are estimated as the minimum values next to the viewpoint before
##' the signal rises again.
##' @param fitFun Fit function to be used. Either \code{distFitMonotoneSymmetric} (default), 
##' \code{distFitMonotone} or any self-defined function to fit the distance dependency.
##' @param sdFun Function to calculate the variation of the data. Default is \code{mad}.
##' @param ... Additional Parameters passed to the distance fit function specified in
##' fitFun  (e.g. \code{distFitMonotone}, \code{distFitMonotoneSymmetric}).
##' @return Returns a \code{FourC} object for the selected viewpoint
##' with z-score values for all fragments on the viewpoint chromosome that 
##' passed the \code{minCount} threshold and that were not too close to the 
##' viewpoint. All additional required data is saved in the object. Especially
##' the following information is added to the \code{FourC} object:
##' \enumerate{
##'      \item{exptData:} parameters passed to the getZScore function.
##'      \item{colData:} sdFun values calculated from the fit residuals
##'      \item{rowRanges:} the distance information to the viewpoint and information
##'                      calculated for the variance stabilizing data are added.
##'      \item{assays:} variance stabilized count values (trafo), fit values (fit),
##'                     z-score values (zScore), associated p-values (pValue) and
##'                     adjusted p-values (pAdjusted)
##' }
##' 
##' 
##' @author Felix A. Klein, \email{felix.klein@@embl.de}
##' @seealso \code{\link{FourC}}, \code{\link{countFragmentOverlaps}},
##' \code{\link{distFitMonotone}}, \code{\link{distFitMonotoneSymmetric}}
##' @examples 
##'                        
##' data(fc, package="FourCSeq")
##' 
##' fcf <- getZScores(fc)
##' 
##' fcf
##' 
##'  
##' @export
getZScores <- function(object,
                       removeZeros = TRUE,                             
                       minCount=40,
                       minDist=NULL,
                       fitFun="distFitMonotoneSymmetric",
                       sdFun=mad,
                       ...){  
  stopifnot(class(object)=="FourC")
  
  if( ! c("counts") %in% names(assays(object)))
    stop("No assay 'counts' found. Use 'combineFragEnds' first.")

  if(any(!c("chr", "start", "end") %in% names(colData(object))))
    stop("No information about viewpoint position provided. Add this information as described in the vignette.")
    
  if( c("zScore") %in% names(assays(object)))
    stop("z-scores are already calculated. To recalculate z-scores use the object returned by 'combineFragEnds'.")
  
  ## get the viewpoint data
  viewpoint = unique(colData(object)$viewpoint)
  if(length(viewpoint) != 1)
    stop("None or more than one viewpoint are contained in the 'FourC' object.
         Use a 'FourC' object that contains only one viewpoint.")
  print(viewpoint)
  
  ## just keep viewpoint chromosome
  object <- object[seqnames(object) == unique(colData(object)$chr),]
    
  #######################################################
  ## calculated distances around the current viewpoint
  #######################################################  
  fragData = getDistAroundVp(viewpoint, colData(object), rowRanges(object))
  
  ## calculate the mean counts for each fragments for filtering
  medianCounts <- apply(counts(object), 1, median)
    
  ## remove counts of fragments to close to the viewpoint
  if(!is.null(minDist)){
    tooClose = abs(fragData$dist) <= minDist
  } else {
    ## find min to left
    toLeft <- fragData$dist > -20000 & fragData$dist < 0 & !is.na(fragData$dist)
    afterMin <- 1:sum(toLeft) > tail(which(sign(diff(medianCounts[toLeft])) < 0), 1) + 1
    toExclude <- which(toLeft)[afterMin]
    
    ## find min to right
    toRight <- fragData$dist < 20000 & fragData$dist > 0 & !is.na(fragData$dist)
    beforeMin <- 1:sum(toRight) < which(sign(diff(medianCounts[toRight])) > 0)[1]
    toExclude <- c(toExclude, which(toRight)[beforeMin])
    
    ## make sure minimum distance from the viewpoint is 1000
    toExclude = union(toExclude, which(abs(fragData$dist) < 1000))
    
    tooClose = rep(FALSE, length(fragData$dist))
    tooClose[toExclude] = TRUE
  }
  
  ## remove points to close to the viewpoint and store the info for fitting 
  fragData$posLeft[tooClose] = FALSE
  fragData$posRight[tooClose] = FALSE
  fragData$tooClose = tooClose
  
  ## remove low count values and store the info for fitting
  lowCounts = medianCounts < minCount
  fragData$posLeft[lowCounts] = FALSE
  fragData$posRight[lowCounts] = FALSE
  fragData$lowCounts = lowCounts
  
  ## select only values that fullfill the filter criterion & have loggeomeans finite
  fragData$selectedForFit <- (fragData$posLeft | fragData$posRight) 
  
  ## update metadata 
  newCols <- c("tooClose", "lowCounts", "selectedForFit")
  mcolsRows <- DataFrame(type=rep("fragmentSelection", length(newCols)),
                         description=rep("",length(newCols)))
  mcols(mcols(fragData))[colnames(mcols(fragData)) %in% newCols,] <- mcolsRows
  
  ## replace rowRanges with fragData containing distance information and convert fragData to data.frame
  rowRanges(object) <- fragData
  fragData = as.data.frame(fragData)
  
  ## factor condition
  colData(object)$condition <- factor(colData(object)$condition)
  
  ###############################
  ## calculate z-scores
  ###############################
  
  ## select only fragments that passed the filter
  dds <- object[mcols(object)$selectedForFit,]
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersions(dds)

  if (attr(dispersionFunction(dds), "fitType") != "parametric") {
    stop("Failed to estimate the parameters of the Variance stabilizing transformation.")
  } else {
    coefs <- attr(dispersionFunction(dds), "coefficients")
    attr(dds, "vst") <- function(q) {
      log((1 + coefs["extraPois"] + 2 * coefs["asymptDisp"] * 
             q + 2 * sqrt(coefs["asymptDisp"] * q * (1 + coefs["extraPois"] + 
                                                       coefs["asymptDisp"] * q)))/(4 * coefs["asymptDisp"]))/log(2)}
    attr(dds, "inverse-vst") <- function(q){
      (4*coefs["asymptDisp"]*2^q-(1+coefs["extraPois"]))^2 / 
        (4*coefs["asymptDisp"]*(1+coefs["extraPois"] + 
                                  (4*coefs["asymptDisp"]*2^q-(1+coefs["extraPois"]))))}
  }

  ## get transformed values
  trafo <- getVST(dds)(counts(dds))

  fit <- apply(trafo,
               2,
               fitFun,
               fragData=as.data.frame(rowRanges(dds)),
               removeZeros=removeZeros,
               ...)
  residuals <- trafo - fit
  sd <- apply(residuals, 2, sdFun)

  zScore <- sweep(residuals, 2, sd, "/")
  
  ## calculate p-values
  pValue <- apply(zScore, 2, pnorm, lower.tail=FALSE)
  pAdjusted <- apply(pValue,2, p.adjust, method="BH")
    
  ## update colData with sd values
  colData(dds)$sd <- sd
  idx <- which(colnames(colData(dds)) == "sd")
  metaDataFrame <- DataFrame(type="intermediate",
                             description="sd/mad calculated from the residuals")
  mcols(colData(dds))[idx,] <- metaDataFrame
  
  ## save the used parameters of the function call in exptData
  exptData(dds)$parameter <- DataFrame(fitFun=fitFun,
                                       removeZeros=removeZeros,
                                       minCount=minCount,
                                       sdFun=sdFun,
                                       ...)
  if(!is.null(minDist)) exptData(dds)$parameter$minDist = minDist
  
  ## add assays
  assays(dds) <- c(assays(dds), SimpleList(trafo=trafo,
                                           fit=fit,
                                           zScore=zScore,
                                           pValue=pValue,
                                           pAdjusted=pAdjusted))
  
  invisible(dds)
}

##' Get normalization factors for each fragment
##' 
##' \code{getNormalizationFactors} uses the distance dependency to calculate 
##' normalization factors for each fragment.
##' 
##'                
##' @param object A \code{FourC} object.
##' @return Returns an updated \code{FourC} object that contains normalization
##' factors.
##' @author Felix A. Klein, \email{felix.klein@@embl.de}
##' @seealso \code{\link{getDifferences}}
##' @examples
##'                        
##' data(fcf, package="FourCSeq")
##' 
##' normalizationFactors <- getNormalizationFactors(fcf)
##' head(normalizationFactors)
##' 
##' @export
getNormalizationFactors <- function(object){
  stopifnot(class(object)=="FourC") 
  if(! ("fit" %in% names(assays(object)))) stop("Call getZScores first!")
   
  
  fit <- getVST(object, inverse=TRUE)(assays(object)[["fit"]])
  
  ## estimate size factors for each fragment
  loggeomeans <- rowMeans(log(fit))

  normFactors <- exp(log(fit) - loggeomeans)

  normalizationFactors(object) <- normFactors
  if("sizeFactor" %in% names(colData(object))) colData(object) <- colData(object)[!names(colData(object)) == "sizeFactor"]
  
  invisible(object)
}



##' Fit the distance dependency
##' 
##' \code{distFitMonotone} takes the variance stabilized count values and 
##' calculates a monotone fit for the distance dependency. The position 
##' information about the viewpoint is stored in \code{fragData}. 
##' The signal trend is fitted for the data left and right to the viewpoint 
##' separately using the fda package. 
##' 
##' @param count Matrix of variance stabilized count data in the defined interval around
##' the viewpoint.
##' @param fragData Data frame with all the information on restriction
##' fragments and the interval around the viewpoint.
##' @param alpha Approximate number of fragments between two breaks of the B-spline
##' basis. \code{floor((number of fragments))/alpha)} is passed to \code{create.bspline.basis} 
##' as \code{nbasis} argument (if this value is smaller than 4, 4 is the default value for \code{nbasis}). 
##' Default is \code{20}.
##' @param penalty Penalty term passed to  \code{fdPar} as \code{lambda} argument.
##' Default is \code{0.1}.
##' @param removeZeros Defines if zero values should be removed from the fit. 
##' Default is \code{FALSE}.
##' @param ... Additional parameters.
##' @return Returns a matrix with z-score values of all fragments around the given
##' viewpoint.
##' @author Felix A. Klein, \email{felix.klein@@embl.de}
##' 
##' @export
distFitMonotone <- function (count, fragData, alpha = 20, penalty=0.1, removeZeros=FALSE, ...) 
{

  if(removeZeros){
    fragData$posLeft <- fragData$posLeft & count > min(count)
    fragData$posRight <- fragData$posRight & count > min(count)
  }
  
  fragData$logDist = log(abs(fragData$dist))

  ## left side
  basisobj <- create.bspline.basis(c(min(fragData$logDist[fragData$posLeft]),max(fragData$logDist[fragData$posLeft])), 
                                   nbasis=max(floor(sum(fragData$posLeft)/alpha), 4))
  fdParobj = fdPar(basisobj, 2, penalty)
  monoLeft = smooth.monotone(argvals=fragData$logDist[fragData$posLeft], y=count[fragData$posLeft], WfdParobj=fdParobj)
  betaLeft <- monoLeft$beta

  ## right side 
  basisobj <- create.bspline.basis(c(min(fragData$logDist[fragData$posRight]),max(fragData$logDist[fragData$posRight])), 
                                   nbasis=max(floor(sum(fragData$posRight)/alpha), 4))
  fdParobj = fdPar(basisobj, 2, penalty)
  monoRight = smooth.monotone(argvals=fragData$logDist[fragData$posRight], y=count[fragData$posRight], WfdParobj=fdParobj)
  betaRight <- monoRight$beta

  fit = rep(NA, length(count))
  selLeft <- fragData$dist < 0
  if ("tooClose" %in% names(fragData)) 
    selLeft <- selLeft & !fragData$tooClose
  fit[selLeft] <- predict(monoLeft, fragData$logDist[selLeft])
  ## extend the minimum value to the end
  fit[selLeft][!is.finite(fit[selLeft])] <- min(fit[selLeft], na.rm=TRUE)

  selRight <- fragData$dist > 0
  if ("tooClose" %in% names(fragData)) 
    selRight <- selRight & !fragData$tooClose
  fit[selRight] <- predict(monoRight, fragData$logDist[selRight])
  ## extend the minimum value to the end  
  fit[selRight][!is.finite(fit[selRight])] <- min(fit[selRight], na.rm=TRUE)
  
  invisible(fit)
}


##' Fit the distance dependency
##' 
##' \code{distFitMonotoneSymmetric} takes the variance stabilized count values and 
##' calculates a symmetric monotone fit for the distance dependency. The position 
##' information about the viewpoint is stored in \code{fragData}. 
##' The signal trend is fitted for the combined data left and right to the viewpoint
##' using the fda package. 
##' 
##' @param count Matrix of variance stabilized count data in the defined interval around
##' the viewpoint.
##' @param fragData Data frame with all the information on restriction
##' fragments and the interval around the viewpoint.
##' @param alpha Approximate number of fragments between two breaks of the B-spline
##' basis. \code{floor((max(number of fragments left or right))/alpha)} is passed to 
##' \code{create.bspline.basis} as \code{nbasis} argument (if this value is smaller 
##' than 4, 4 is the default value). Default is \code{20}.
##' @param penalty Penalty term passed to  \code{fdPar} as \code{lambda} argument.
##' Default is \code{0.1}.
##' @param removeZeros Defines if zero values should be removed from the fit. 
##' Default is \code{FALSE}.
##' @param ... Additional parameters.
##' @return Returns a matrix with z-score values of all fragments around the given
##' viewpoint.
##' @author Felix A. Klein, \email{felix.klein@@embl.de}
##' 
##' @export
distFitMonotoneSymmetric <- function (count, fragData, alpha = 20, penalty=0.1, removeZeros=FALSE, ...) 
{
  
  if(removeZeros){
    fragData$posLeft <- fragData$posLeft & count > min(count)
    fragData$posRight <- fragData$posRight & count > min(count)
  }
  
  fragData$logDist = log(abs(fragData$dist))
  
  ## left side
  basisobj <- create.bspline.basis(c(min(fragData$logDist[fragData$posLeft | fragData$posRight]),
                                     max(fragData$logDist[fragData$posLeft | fragData$posRight])), 
                                   nbasis=max(floor(max(sum(fragData$posLeft), sum(fragData$posRight))/alpha), 4))
  fdParobj = fdPar(basisobj, 2, penalty)
  mono = smooth.monotone(argvals=fragData$logDist[fragData$posLeft | fragData$posRight], 
                             y=count[fragData$posLeft | fragData$posRight], 
                             WfdParobj=fdParobj)
  beta <- mono$beta

  
  fit = rep(NA, length(count))
  sel <- fragData$dist != 0
  if ("tooClose" %in% names(fragData)) 
    sel <- sel & !fragData$tooClose
  fit[sel] <- predict(mono, fragData$logDist[sel])
  ## extend the minimum value to the end
  fit[sel][!is.finite(fit[sel])] <- min(fit[sel], na.rm=TRUE)
  
  invisible(fit)
}


getVST <- function(object, inverse=FALSE){
  if(!is.null(attr(object, "vst")) & !is.null(attr(object, "inverse-vst"))){
    if(!inverse) return(attr(object, "vst"))
    return(attr(object, "inverse-vst"))
  } else {
    stop("No variance stabilizing transformation available. Call getZScores first.")
  }
}


##' FourCSeq analysis results
##'
##' Using the \code{DESEeq2} function \code{results},
##' \code{getAllResults} extracts results from a FourCSeq analysis
##' giving base means across samples,
##' log2 fold changes, standard errors, test statistics,
##' p-values and adjusted p-values for all pair-wise
##' conditions tested.
##' 
##' @param object \code{FourC} object for which z-scores and differences
##' have been calculated.
##' @param ... Additional parameters that can be passed to \code{results}.
##' @return DataFrame of results columns for all pairwise tests, with metadata
##' columns of coefficient and test information
##' 
##' @author Felix A. Klein, \email{felix.klein@@embl.de}
##' 
##' @seealso \code{\link[DESeq2]{results}}
##' 
##' @examples
##' 
##'                        
##' data(fcf, package="FourCSeq")
##' 
##' fcf <- getDifferences(fcf, referenceCondition="WE_68h")
##' 
##' results <- getAllResults(fcf)
##' results
##' 
##' @export
getAllResults <- function(object, ...){
  stopifnot(class(object)=="FourC") 
  if(!"results" %in% mcols(mcols(object))$type) stop("Call getDifferences first!")
  
  object <- object[rowRanges(object)$selectedForFit,]
  lvls <- levels(colData(object)$condition)
  
  if(length(lvls) > 1){
    combi <- combinations(length(lvls), 2)
    
    res <- do.call(cbind, lapply(seq_len(nrow(combi)), function(i, ...){
      tmp <- results(object, contrast=c("condition", lvls[combi[i,]]), ...)[,-1]
      colnames(tmp) <- paste(colnames(tmp),
                             paste(lvls[combi[i,]], collapse="_"),
                             sep="_")
      tmp
    }, ...))
    res = cbind(mcols(object)[,"baseMean", drop=FALSE], res)
  } else stop("Only one level and no differential results")
  invisible(res)
}

##' Add peaks based on z-scores and adjusted p-values
##' 
##' Fragments are annotated to be a z-score peak (TRUE or FALSE) based 
##' on z-scores and fdr values for the given fragment
##' 
##' @param object \code{FourC} object for which z-scores have been calculated.
##' @param zScoreThresh Threshold on z-score that has to met by all replicates 
##' of a condition.
##' @param fdrThresh Threshold on adjusted p-value that has to be met in at 
##' least one replicated of a condition.
##' @return \code{FourC} object with an additional assay 'peaks' of boolean 
##' values. The parameters of the call to addPeaks are stored as DataFrame 
##' \code{peakParameter} in the exptData slot.
##' 
##' @author Felix A. Klein, \email{felix.klein@@embl.de}
##' 
##' data(fcf, package="FourCSeq")
##' 
##' fcf <- addPeaks(fcf)
##' 
##' @export
addPeaks <- function(object, zScoreThresh = 2, fdrThresh = 0.1){
  stopifnot(class(object)=="FourC") 
  if(! ("zScore" %in% names(assays(object)))) stop("Call getZScores first!")
  
  zScore <- assay(object, "zScore")
  fdr <- assay(object, "pAdjusted")
  
  peaks = array(FALSE, dim=dim(fdr), dimnames=dimnames(fdr))
  
  for(lvl in levels(colData(object)$condition)){
    condCols <- colData(object)$condition == lvl 
    tmp <- (apply(fdr[,condCols, drop=FALSE] < fdrThresh, 1, any, na.rm=TRUE) & 
               apply(zScore[,condCols, drop=FALSE] > zScoreThresh, 1, all, na.rm=TRUE))    
    ## annotate peaks
    peaks[tmp,condCols] = TRUE
  }

  assays(object)[["peaks"]] <- peaks
  
  exptData(object)$peakParameter <- DataFrame(zScoreThresh = zScoreThresh, 
                                              fdrThresh = fdrThresh)
  
  invisible(object)
}


##' Smooth the counts of neighboring fragments
##' 
##' Counts are smoothed using the number of fragments provided by binWidth.
##' binWidth has to be an odd number so that an equal number of fragments
##' to each side of the current fragment are used for smoothing.
##' 
##' @param object A \code{FourC} object.
##' @param assay Assay name that will be smoothed.
##' @param binWidth Integer vector of odd numbers.
##' @return Returns an updated \code{FourC} object with smoothed counts
##' for each binWidth as new assays.
##' @author Felix A. Klein, \email{felix.klein@@embl.de}
##' @seealso \code{\link{FourC}}
##' @examples 
##' 
##' data(fc, package="FourCSeq")
##' 
##' fc <- smoothCounts(fc)
##' fc
##' 
##' @export
smoothCounts <- function(object,
                         assay="counts",
                         binWidth = 5){
  
  stopifnot(class(object)=="FourC")
  
  if(any( ! c("counts") %in% names(assays(object))))
    stop("Use combineFragEnds to combine the counts for left and right fragment ends.")
  if(any( binWidth %% 2 == 0))
    stop("binWidth must only contain odd numbers.")
  
  toSmooth = assay(object, assay)
  
  smoothed = lapply((binWidth - 1)/2, function(w){
    smooth = array(0, dim=dim(toSmooth), dimnames=dimnames(toSmooth))
    for(chr in unique(seqnames(rowRanges(object)))){
      cat(2*w + 1, chr, "\n")
      currentChr = as.logical(seqnames(rowRanges(object)) == chr)
      if(sum(currentChr) > 2*w + 1){
        band = bandSparse(n=sum(currentChr), m=sum(currentChr), k=(-w):w)
        band = band / rowSums(band)
        smooth[currentChr,] = matrix(band %*% toSmooth[currentChr,])
      } else {
        cat("Smoothing window to larg than chromosome size.\n")
        cat("Values are replaces by mean values.\n")
        smooth[currentChr,] <- t(apply(toSmooth[currentChr,,drop=FALSE], 2, mean))
      }
    }
    invisible(smooth)
  })
  
  names(smoothed) <- paste0(assay, "_", binWidth)  
  smoothed <- SimpleList(smoothed)  
  assays(object) <- c(assays(object), smoothed)
  
  invisible(object)
}


##' Smooth the hits of neighboring fragments
##' 
##' Counts are transformed to a hit if they exceed the given threshold.
##' Hits are then smoothed using the number of fragments provided by binWidth.
##' binWidth has to be an odd number so that an equal number of fragments
##' to each side of the current fragment are used for smoothing.
##' 
##' @param object A \code{FourC} object.
##' @param binWidth Integer vector of odd numbers.
##' @param thresh Single integer defining the threshold for calling a fragment a hit.
##' @return Returns an updated \code{FourC} object with smoothed counts
##' for each binWidth as new assays.
##' @author Felix A. Klein, \email{felix.klein@@embl.de}
##' @seealso \code{\link{FourC}}
##' @examples 
##' 
##' data(fc, package="FourCSeq")
##' 
##' fc <- smoothHitPerCent(fc)
##' fc
##' 
##' @export
smoothHitPerCent <- function(object,
                             binWidth = 101,
                             thresh=1){
  
  stopifnot(class(object)=="FourC")
  
  if(any( ! c("counts") %in% names(assays(object))))
    stop("Use combineFragEnds to combine the counts for left and right fragment ends.")
  if(any( binWidth %% 2 == 0))
    stop("binWidth must only contain odd numbers.")
  
  hits = abs(counts(object) > thresh)
  
  smoothed = lapply((binWidth - 1)/2, function(w){
    smooth = array(0, dim=dim(hits), dimnames=dimnames(hits))
    for(chr in unique(seqnames(rowRanges(object)))){
      cat(2*w + 1, chr, "\n")
      currentChr = as.logical(seqnames(rowRanges(object)) == chr)
      if(sum(currentChr) > 2*w + 1){
        band = bandSparse(n=sum(currentChr), m=sum(currentChr), k=(-w):w)
        band = band / rowSums(band)
        smooth[currentChr,] = matrix(band %*% hits[currentChr,])
      } else {
        cat("Smoothing window larger than chromosome size.\n")
        cat("Values are replaces by mean values.\n")
        tmp <- t(smooth[currentChr,])
        tmp[,] <- apply(hits[currentChr,,drop=FALSE], 2, mean)
        smooth[currentChr,] <- t(tmp)
      }
    }
    invisible(smooth)
  })
  
  names(smoothed) <- paste0("hitPerCent_", binWidth, "_", thresh)  
  smoothed <- SimpleList(smoothed)  
  assays(object) <- c(assays(object), smoothed)
  
  invisible(object)
}
