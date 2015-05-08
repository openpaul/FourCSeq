##' Function to read reference sequences
##' 
##' This functions allows to retrieve the reference sequence.
##' 
##' 
##' @param object A \code{FourC} object.
##' @return A \code{DNAStringSet} object containing the sequences of the
##' reference genome from the \code{FaFile} or \code{BSgenome} object.
##' @author Felix A. Klein
##' @examples
##' 
##' 
##' 
##' metadata <- list(projectPath=tempdir(),
##'                  fragmentDir="re_fragments",
##'                  referenceGenomeFile=system.file("extdata/dm3_chr2L_1-6900.fa", 
##'                                                  package="FourCSeq"),
##'                  reSequence1="GATC",
##'                  reSequence2="CATG",
##'                  primerFile=system.file("extdata/primer.fa", 
##'                                         package="FourCSeq"),
##'                  bamFilePath=system.file("extdata/bam", package="FourCSeq"))
##' 
##' colData <- DataFrame(viewpoint = "testdata", 
##'                      condition = factor(rep(c("WE_68h", "MESO_68h", "WE_34h"),                    
##'                                             each=2),
##'                                         levels = c("WE_68h", "MESO_68h", "WE_34h")),
##'                      replicate = rep(c(1, 2), 
##'                                      3),
##'                      bamFile = c("CRM_ap_ApME680_WE_6-8h_1_testdata.bam", 
##'                                  "CRM_ap_ApME680_WE_6-8h_2_testdata.bam",       
##'                                  "CRM_ap_ApME680_MESO_6-8h_1_testdata.bam", 
##'                                  "CRM_ap_ApME680_MESO_6-8h_2_testdata.bam", 
##'                                  "CRM_ap_ApME680_WE_3-4h_1_testdata.bam",
##'                                  "CRM_ap_ApME680_WE_3-4h_2_testdata.bam"),
##'                      sequencingPrimer="first")                     
##'                     
##'
##' fc <- FourC(colData, metadata)
##' fc
##' 
##' refSeq <- getReferenceSeq(fc)
##' 
##' 
##' @export
getReferenceSeq <- function (object){
  stopifnot(class(object)=="FourC")
  refGenFile = metadata(object)$referenceGenomeFile
  
  if(class(refGenFile) == "BSgenome"){
    ref = getSeq(refGenFile)
    names(ref) = seqnames(refGenFile)
  } else if(class(refGenFile) == "character"){
    ref = getSeq(FaFile(refGenFile))
    names(ref) = sapply(strsplit(names(ref), " ", fixed=TRUE), 
                        "[", 1)
  } else {
    stop("referenceGenomeFile hast to be of class 'BSgenome' or 'character'.")
  }
  if(any(duplicated(names(ref)))) 
    warning("There are duplicated chromosome names in the reference.")
  ref
}

########################################################


##' Find the fragments to which the viewpoint primers map.
##' 
##' \code{findViewpointFragments} finds the position of the viewpoint primer in
##' the reference genome and on which restriction fragment they fall.
##' It saves the results in these files: 
##'  - projectPath/fragmentDir/primerFragments.rda
##'  - projectPath/fragmentDir/primerFragments.txt
##' 
##' @usage findViewpointFragments(object)
##' 
##' @param object A \code{FourC} object.
##' @return NULL
##' @author Felix A. Klein, \email{felix.klein@@embl.de}
##' @seealso \code{\link{FourC}}
##' \code{\link{getReferenceSeq}}
##' @examples
##' 
##' 
##' 
##' metadata <- list(projectPath=tempdir(),
##'                  fragmentDir="re_fragments",
##'                  referenceGenomeFile=system.file("extdata/dm3_chr2L_1-6900.fa", 
##'                                                  package="FourCSeq"),
##'                  reSequence1="GATC",
##'                  reSequence2="CATG",
##'                  primerFile=system.file("extdata/primer.fa", 
##'                                         package="FourCSeq"),
##'                  bamFilePath=system.file("extdata/bam", package="FourCSeq"))
##' 
##' colData <- DataFrame(viewpoint = "testdata", 
##'                      condition = factor(rep(c("WE_68h", "MESO_68h", "WE_34h"),                    
##'                                             each=2),
##'                                         levels = c("WE_68h", "MESO_68h", "WE_34h")),
##'                      replicate = rep(c(1, 2), 
##'                                      3),
##'                      bamFile = c("CRM_ap_ApME680_WE_6-8h_1_testdata.bam", 
##'                                  "CRM_ap_ApME680_WE_6-8h_2_testdata.bam",       
##'                                  "CRM_ap_ApME680_MESO_6-8h_1_testdata.bam", 
##'                                  "CRM_ap_ApME680_MESO_6-8h_2_testdata.bam", 
##'                                  "CRM_ap_ApME680_WE_3-4h_1_testdata.bam",
##'                                  "CRM_ap_ApME680_WE_3-4h_2_testdata.bam"),
##'                      sequencingPrimer="first")
##'
##' fc <- FourC(colData, metadata)
##' fc
##' 
##' fc <- addFragments(fc)
##' 
##' findViewpointFragments(fc) 
##' 
##' fc <- addViewpointFrags(fc)
##' fc
##' 
##' 
##' @export
findViewpointFragments <- function(object){
  stopifnot(class(object)=="FourC")
  
  if(length(rowRanges(object)) == 0) 
      stop("Add fragments before calling 'findViewpointFragments'")

  projectPath = metadata(object)$projectPath
  primerFile = metadata(object)$primerFile

  frag <- rowRanges(object)
  
  ## read primer sequences
  primer = getSeq(FaFile(primerFile))

  ## read reference genome
  ref = getReferenceSeq(object)
    
  ## match primers to reference genome
  m = unlist(GRangesList(lapply(primer, 
                                function(i){
                                  strand="+"
                                  tmp=unlist(vmatchPattern(toString(i), ref))
                                  if(!length(start(tmp))){
                                    tmp=unlist(vmatchPattern(
                                      toString(reverseComplement(i)), ref))
                                    strand="-"
                                  } 
                                  gr <- GRanges(seqnames=Rle(names(tmp)),
                                                tmp, 
                                                strand=strand)      
                                  names(gr) = names(i)
                                  gr
                                })))
  
  ## find overlaps with fragments
  ov = findOverlaps(m, frag)

  ## format into data frame for output 
  vpFrag <- cbind(data.frame(viewpoint=names(m)), 
                  as.data.frame(frag[subjectHits(ov)]))
  names(vpFrag)[names(vpFrag) == "seqnames"] = "chr"
  vpFrag$strand=NULL
  
  vpFrag$primerStart = start(m)
  vpFrag$primerEnd = end(m)
  startDiff = abs(vpFrag$start - vpFrag$primerStart)
  endDiff = abs(vpFrag$end - vpFrag$primerEnd)
  vpFrag$primerSide = c("left", "right")[ifelse(startDiff < endDiff, 1, 2)]
  
  outputDir = file.path(metadata(object)$projectPath,
                        metadata(object)$fragmentDir)
  if(!file.exists(outputDir)) dir.create(outputDir, recursive=TRUE)
  
  save(vpFrag, file=file.path(outputDir, "primerFragments.rda"))

  write.table(vpFrag, 
              file = file.path(outputDir, "primerFragments.txt"), 
              quote = FALSE, 
              sep = "\t", 
              row.names = FALSE)
}
########################################################


##' Add the information of the viewpoint fragments
##' 
##' \code{addViewpointFrags} adds the genomic position and information of the
##' viewpoint primer fragments to \code{colData} of the \code{FourC} object.
##' 
##' @usage addViewpointFrags(object, primerFragFile)
##' @param object A \code{FourC} object.
##' @param primerFragFile character string defining the file that contains the
##' information about the primer fragments. Defaults to
##' \code{"re_fragments/primerFragments.rda"} which is the default output file
##' of \code{findViewpointFragments}.
##' @return Updated \code{FourC} object with information about the viewpoint
##' fragments added to \code{colData}.
##' @author Felix A. Klein, \email{felix.klein@@embl.de}
##' @seealso \code{\link{FourC}}, \code{\link{findViewpointFragments}}
##' @examples
##' 
##' 
##' 
##' metadata <- list(projectPath=tempdir(),
##'                  fragmentDir="re_fragments",
##'                  referenceGenomeFile=system.file("extdata/dm3_chr2L_1-6900.fa", 
##'                                                  package="FourCSeq"),
##'                  reSequence1="GATC",
##'                  reSequence2="CATG",
##'                  primerFile=system.file("extdata/primer.fa", 
##'                                         package="FourCSeq"),
##'                  bamFilePath=system.file("extdata/bam", package="FourCSeq"))
##' 
##' colData <- DataFrame(viewpoint = "testdata", 
##'                      condition = factor(rep(c("WE_68h", "MESO_68h", "WE_34h"),                    
##'                                             each=2),
##'                                         levels = c("WE_68h", "MESO_68h", "WE_34h")),
##'                      replicate = rep(c(1, 2), 
##'                                      3),
##'                      bamFile = c("CRM_ap_ApME680_WE_6-8h_1_testdata.bam", 
##'                                  "CRM_ap_ApME680_WE_6-8h_2_testdata.bam",       
##'                                  "CRM_ap_ApME680_MESO_6-8h_1_testdata.bam", 
##'                                  "CRM_ap_ApME680_MESO_6-8h_2_testdata.bam", 
##'                                  "CRM_ap_ApME680_WE_3-4h_1_testdata.bam",
##'                                  "CRM_ap_ApME680_WE_3-4h_2_testdata.bam"),
##'                      sequencingPrimer="first")
##'
##' fc <- FourC(colData, metadata)
##' fc
##' 
##' fc <- addFragments(fc)
##' 
##' findViewpointFragments(fc) 
##' 
##' fc <- addViewpointFrags(fc)
##' fc
##' 
##' 
##' @export
addViewpointFrags <- function(object,
                              primerFragFile="re_fragments/primerFragments.rda"){

  ## read viewpoint fragment data
  if(!file.exists(file.path(metadata(object)$projectPath, primerFragFile)))
    stop("Call 'findViewpointFragments' first.")

  load(file.path(metadata(object)$projectPath, primerFragFile))
  
  if(sum(dim(colData(object))) == 0 | !("viewpoint" %in% names(colData(object))))
    stop("colData must contain the column 'viewpoint'
       with the names of the viewpoints")
  if(!("viewpoint" %in% names(vpFrag)))
    stop("colData must contain the column 'viewpoint'
       with the names of the viewpoints")  

  ## add mcols and convert to DataFrame
  vpFrag <- DataFrame(vpFrag)
  mcols(vpFrag) <- DataFrame(type=rep("viewpointInfo",ncol(vpFrag)),
                             description=c("viewpoint", 
                                           "viewpoint chromosome",
                                           "viewpoint fragment start position",
                                           "viewpoint fragment end position",
                                           "viewpoint fragment size",
                                           "viewpoint fragment end size until cutting site of second cutter",
                                           "viewpoint fragment end size until cutting site of second cutter",
                                           "valid fragment end",
                                           "valid fragment end",
                                           "primer start position",
                                           "primer end position",
                                           "side of the primer on the viewpoint fragment"))
  
  ## merge viewpoint data with the data about the viewpoint fragments 
  ## and prevent character to factor conversion
  matchedVp <- match(colData(object)$viewpoint, vpFrag$viewpoint)
  if(any(is.na(matchedVp))) stop(paste(paste(colData(object)$viewpoint, collapse=", "),
                                       "contains viewpoint names not present in the primer names:\n",
                                       paste(vpFrag$viewpoint, collapse=", ")))
                                       
  matchedVpFrag <- do.call(rbind, lapply(matchedVp, function(row) vpFrag[row,]))
  row.names(matchedVpFrag)=row.names(colData(object))
  
  colData(object) <- c(colData(object), matchedVpFrag[,-1])

  object
}


##########################################################
getFragments <- function(rePattern, 
                         ref){
  sites = getSites(rePattern, ref)
  wholeRef = GRanges(seqnames=Rle(seqnames(seqinfo(sites))),
                     IRanges(1,
                             seqlengths(sites)))
  setdiff(wholeRef, sites)
}


########################################################
getSites <- function(rePattern, 
                     ref){
  if(class(ref)=="BSgenome"){
    sites = vmatchPattern(rePattern, ref, 
                          exclude = seqnames(ref)[grep("random", seqnames(ref))])
    strand(sites) = "*"
    sites <- unique(sites)
    seqlevels(sites) = as.character(unique(seqnames(sites)))
    seqinfo(sites) <- seqinfo(sites)[seqlevels(sites)]    
    names(sites) = 1:length(sites)
  } else {
    refinfo =  Seqinfo(seqnames=names(ref),
                       seqlengths=width(ref))
    sites = unlist(vmatchPattern(rePattern, ref))
    sites = GRanges(seqnames=Rle(names(sites)),
                    sites,
                    seqinfo=refinfo)
    names(sites) = 1:length(sites)
  }
  sites
}
###########################################

##' Add the restriction fragment information
##' 
##' \code{addFragments} adds the genomic position and information of the
##' restriction fragments as \code{GRanges} object to \code{rowRanges} of 
##' the \code{FourC} object.
##' 
##' @param object A \code{FourC} object.
##' @param minSize Minimum size of a restriction fragment end to be valid.
##' Default is \code{20} bases.
##' @param filter Defines whether fragments that do not contain a cutting site
##' of the second restriction enzyme or are smaller than \code{minSize} should
##' be filtered out.
##' @param save Defines if the fragment information should be saved as txt and bed files 
##' in the fragmentDir folder of the projectPath.
##' @return Updated \code{FourC} object that contains the  information about 
##' the restriction fragments in \code{rowRanges}.
##' @author Felix A. Klein, \email{felix.klein@@embl.de}
##' @seealso \code{\link{FourC}}, \code{\link{findViewpointFragments}}
##' @examples
##' 
##' 
##' 
##' metadata <- list(projectPath=tempdir(),
##'                  fragmentDir="re_fragments",
##'                  referenceGenomeFile=system.file("extdata/dm3_chr2L_1-6900.fa", 
##'                                                  package="FourCSeq"),
##'                  reSequence1="GATC",
##'                  reSequence2="CATG",
##'                  primerFile=system.file("extdata/primer.fa", 
##'                                         package="FourCSeq"),
##'                  bamFilePath=system.file("extdata/bam", package="FourCSeq"))
##' 
##' colData <- DataFrame(viewpoint = "testdata", 
##'                      condition = factor(rep(c("WE_68h", "MESO_68h", "WE_34h"),                    
##'                                             each=2),
##'                                         levels = c("WE_68h", "MESO_68h", "WE_34h")),
##'                      replicate = rep(c(1, 2), 
##'                                      3),
##'                      bamFile = c("CRM_ap_ApME680_WE_6-8h_1_testdata.bam", 
##'                                  "CRM_ap_ApME680_WE_6-8h_2_testdata.bam",       
##'                                  "CRM_ap_ApME680_MESO_6-8h_1_testdata.bam", 
##'                                  "CRM_ap_ApME680_MESO_6-8h_2_testdata.bam", 
##'                                  "CRM_ap_ApME680_WE_3-4h_1_testdata.bam",
##'                                  "CRM_ap_ApME680_WE_3-4h_2_testdata.bam"),
##'                      sequencingPrimer="first")
##'
##' fc <- FourC(colData, metadata)
##' fc
##' 
##' fc <- addFragments(fc)
##' fc
##' 
##' @export
addFragments <- function(object,
                         minSize=20,
                         filter=TRUE,
                         save=FALSE){
  
  stopifnot(class(object)=="FourC")
  
  ref = getReferenceSeq(object)

  frag = getFragments(metadata(object)$reSequence1, ref)
  site = getSites(metadata(object)$reSequence2, ref)  
  
  ## get distance of second cutter to fragment end
  fragWithSite = subsetByOverlaps(frag, site, minoverlap = nchar(metadata(object)$reSequence2))
  
  firstSites = site[findOverlaps(fragWithSite, site, select="first")]
  lastSites = site[findOverlaps(fragWithSite, site, select="last")]
  
  ## get size of fragments after second digestion
  mcols(fragWithSite)$leftSize = start(firstSites) - start(fragWithSite)
  mcols(fragWithSite)$rightSize = end(fragWithSite) - end(lastSites)
  
  ## check for possible overlaps of cutters and set them to -1
  mcols(fragWithSite)$leftSize[start(firstSites) < start(fragWithSite)] = -1
  mcols(fragWithSite)$rightSize[end(fragWithSite) < end(lastSites)] = -1
  
  mcols(fragWithSite)$leftValid = mcols(fragWithSite)$leftSize >= minSize
  mcols(fragWithSite)$rightValid = mcols(fragWithSite)$rightSize >= minSize  
  
  ## get invalid fragments that do not contain any full cutting site of the second cutter
  fragNoSite = setdiff(frag, fragWithSite)
  mcols(fragNoSite)$leftSize = -1
  mcols(fragNoSite)$rightSize = -1
  mcols(fragNoSite)$leftValid = FALSE
  mcols(fragNoSite)$rightValid = FALSE
  
  ## combine all fragments
  frag = sort(c(fragNoSite, fragWithSite))
  
  if(save){
    tmp <- options(scipen = 12)    
    ## save re fragment information in the specified fragmentDir 
    outputDir = file.path(metadata(object)$projectPath,
                          metadata(object)$fragmentDir)
    if(!file.exists(outputDir)) dir.create(outputDir)
    ## save frags in 0 based half open format for inspection / mapping with python
    output = file.path(outputDir, 
                       "valid_fragments.txt")
    saveGR(object, 
           frag,
           output)
    
    ## save bed tracks of the RE1/2 cutting sites in the bedgraph folder
    site1 <- getSites(metadata(object)$reSequence1, ref)

    ## add score=1 for rtracklayer export
    site1$score = 1
    site$score = 1
    
    reSites = list(site1, site)
    names(reSites) = paste0("re_sites_", 
                            c(metadata(object)$reSequence1,
                              metadata(object)$reSequence2),
                            ".bedGraph")
    
    ivCols <- c("seqnames", "start", "end")
    for (name in names(reSites)) {
      outputFile = file.path(outputDir, name)
      export(reSites[[name]], outputFile, "bedGraph")
    }
    options(tmp)
  }
  
  if(filter) frag = frag[frag$leftValid | frag$rightValid]
  
  mcolsRows <- DataFrame(type=rep("referenceFragments",ncol(mcols(frag))),
                         description=rep("",ncol(mcols(frag))))
  mcols(mcols(frag)) <- mcolsRows
  
  rowRanges(object) <- frag
  object
}


## Function to save genomic ranges object in a tab separated text file.
saveGR <- function(object, 
                   gr, 
                   output,
                   keepCols = c("leftSize", 
                                "rightSize", 
                                "leftValid", 
                                "rightValid")){
  tmp = options(scipen=12)
  cat(paste0("@RE\tre1 pattern: ", 
             metadata(object)$reSequence1,
             "\tre2 pattern: ", 
             metadata(object)$reSequence2),
      file=output,
      sep="\n")
  
  chromLength=seqlengths(gr)
  suppressWarnings(
    write.table(cbind(paste0("@CHROM\t", names(chromLength)), 
                      chromLength), 
                file = output, 
                append = TRUE,
                quote = FALSE, 
                sep = "\t",
                row.names = FALSE,
                col.names = FALSE)
  )
  
  grData=as.data.frame(gr)[,c("seqnames", "start", "end", keepCols)]
  ## shift to 0 based half open system for python
  grData$start=grData$start-1
  names(grData) = gsub("seqnames", "chr", names(grData))
  suppressWarnings(
    write.table(grData, 
                file=output,
                append=TRUE,
                row.names=FALSE, 
                sep="\t", 
                quote=FALSE)
  )
  options(tmp)
}

##' Count fragment overlaps
##' 
##' \code{countFragmentOverlaps} counts the number of reads mapping to each
##' fragment end in \code{rowRanges} of the \code{FourC} object.
##' 
##' @usage countFragmentOverlaps(object, trim=0, minMapq=0, shift=0)
##' @param object A \code{FourC} object.
##' @param trim Number of bases that should be trimmed at the start of a read. 
##' This is necessary if the read still contains the restriction enzyme sequence.
##' Default is \code{0} bases.
##' @param minMapq Minimum mapping quality required for counting the read. 
##' Default is \code{0}. If set to negative values the filtering step is skipped.
##' @param shift Maximum difference in starts or ends between read and fragment
##' positions. Default is \code{0}.
##' @return Updated \code{FourC} object that contains two new \code{assays} 
##' \code{countsLeftFragmentEnd} and \code{countsRightFragmentEnd} with the 
##' count values at the respective fragment end.
##' @author Felix A. Klein, \email{felix.klein@@embl.de}
##' @seealso \code{\link{FourC}}, \code{\link{findViewpointFragments}}, 
##' \code{\link{countFragmentOverlapsSecondCutter}}
##' @examples
##' 
##' 
##' 
##' metadata <- list(projectPath=tempdir(),
##'                  fragmentDir="re_fragments",
##'                  referenceGenomeFile=system.file("extdata/dm3_chr2L_1-6900.fa", 
##'                                                  package="FourCSeq"),
##'                  reSequence1="GATC",
##'                  reSequence2="CATG",
##'                  primerFile=system.file("extdata/primer.fa", 
##'                                         package="FourCSeq"),
##'                  bamFilePath=system.file("extdata/bam", package="FourCSeq"))
##' 
##' colData <- DataFrame(viewpoint = "testdata", 
##'                      condition = factor(rep(c("WE_68h", "MESO_68h", "WE_34h"),                    
##'                                             each=2),
##'                                         levels = c("WE_68h", "MESO_68h", "WE_34h")),
##'                      replicate = rep(c(1, 2), 
##'                                      3),
##'                      bamFile = c("CRM_ap_ApME680_WE_6-8h_1_testdata.bam", 
##'                                  "CRM_ap_ApME680_WE_6-8h_2_testdata.bam",       
##'                                  "CRM_ap_ApME680_MESO_6-8h_1_testdata.bam", 
##'                                  "CRM_ap_ApME680_MESO_6-8h_2_testdata.bam", 
##'                                  "CRM_ap_ApME680_WE_3-4h_1_testdata.bam",
##'                                  "CRM_ap_ApME680_WE_3-4h_2_testdata.bam"),
##'                      sequencingPrimer="first")
##'
##' fc <- FourC(colData, metadata)
##' fc
##' 
##' fc <- addFragments(fc)
##' 
##' findViewpointFragments(fc) 
##' 
##' fc <- addViewpointFrags(fc)
##' fc
##' 
##' fc <- countFragmentOverlaps(fc, trim=4, minMapq=30)
##' fc
##' 
##' @export
countFragmentOverlaps <- function(object, trim=0, minMapq=0, shift=0){
  stopifnot(class(object)=="FourC")
  if(length(rowRanges(object)) == 0) 
    stop("Add fragments before calling 'findViewpointFragments'")
  
  cat("reading bam files\n")
  bamFiles = file.path(metadata(object)$bamFilePath, colData(object)$bamFile)
  colData(object)$originalReads = sapply(bamFiles, function(bamFile) countBam(bamFile)$records)
  reads = lapply(bamFiles, function(bamfile){  
    what <- c("mapq")
    flag <- scanBamFlag(isUnmappedQuery = FALSE)
    param <- ScanBamParam(what=what,flag = flag)
    
    readGAlignments(bamfile, param=param)
  })  
  colData(object)$rawReads = sapply(reads, length)
  
  ## filter bad quality reads
  if(minMapq >= 0){
    reads = lapply(reads, 
                   function(rga, minMapq){                 
                     if(!any(is.na(mcols(rga)$mapq))){
                       return(rga[mcols(rga)$mapq > minMapq])
                     }
                     warning("mapq quality filter could not be applied due to NA values in the mapq quality scores. All reads are used for counting overlaps.")
                     rga
                   },                     
                   minMapq=minMapq)
    colData(object)$lowQualityReads = colData(object)$rawReads - sapply(reads, length)
  }
  
  ## convert to GRanges to certify strand specific counting
  reads = lapply(reads, granges)
  
  ## trim reads if trim > 0
  if(trim>0){
    reads = lapply(reads, function(gr, trim){
      start(gr[strand(gr)=="+"]) = start(gr[strand(gr)=="+"]) + trim
      end(gr[strand(gr)=="-"]) = end(gr[strand(gr)=="-"]) - trim  
      gr  
    }, trim=trim)
  }
  
  cat("calculating overlaps\n")  
  frag <- rowRanges(object)

  strand(frag) <- "+"
  countsLeftFragmentEnd <- sapply(reads, countOverlaps, query=frag, type=c("start"), maxgap=shift)
  
  strand(frag) <- "-"
  countsRightFragmentEnd <- sapply(reads, countOverlaps, query=frag, type=c("end"), maxgap=shift)
  
  colData(object)$mappedReads = apply(countsLeftFragmentEnd, 2, sum) + apply(countsRightFragmentEnd, 2, sum)
  colData(object)$mappingRatio = colData(object)$mappedReads / (colData(object)$rawReads - colData(object)$lowQualityReads)

  metaDataFrame <- DataFrame(type=rep("countInfo", 5),
                             description=rep("", 5))
  idx <- colnames(colData(object)) %in% c("originalReads", 
                                          "rawReads", 
                                          "lowQualityReads",
                                          "mappedReads",
                                          "mappingRatio")
  mcols(colData(object))[idx,] <- metaDataFrame
  
  assays(object) <- SimpleList(countsLeftFragmentEnd=countsLeftFragmentEnd, 
                               countsRightFragmentEnd=countsRightFragmentEnd)
  
  assays(object) <- SimpleList(countsLeftFragmentEnd=countsLeftFragmentEnd, 
                               countsRightFragmentEnd=countsRightFragmentEnd)
  object
}

##' Count fragment overlaps when sequencing was performed from the second 
##' cutting site
##' 
##' \code{countFragmentOverlapsSecondCutter} counts the number of reads mapping to
##' each cutting site of the second cutter and then summarizes them over the
##' fragment ends of the first cutter stored in \code{rowRanges} of the 
##' \code{FourC} object.
##' 
##' @param object A \code{FourC} object.
##' @param extend Defines whether the read start should be extended by the length 
##' of the second cutter recognition sequence to overlap the cutting sites which is 
##' required for counting. If the cutting site has been trimmed with the primer
##' sequence this has to be set to \code{TRUE}. Default is \code{TRUE}.
##' @param minMapq Minimum mapping quality required for counting the read. 
##' Default is \code{0}. If set to negative values the filtering step is skipped.
##' @param shift Maximum difference in starts or ends between read and fragment
##' positions. Default is \code{0}.
##' @return Updated \code{FourC} object that contains two new \code{assays} 
##' \code{countsLeftFragmentEnd} and \code{countsRightFragmentEnd} with the 
##' count values at the respective fragment end.
##' @author Felix A. Klein, \email{felix.klein@@embl.de}
##' @seealso \code{\link{FourC}}, \code{\link{findViewpointFragments}}, 
##' \code{\link{countFragmentOverlaps}}
##' @examples
##' 
##' ## not run:
##' ## countFragmentOverlapsSecondCutter(fc, extend=TRUE, minMapq=30)
##' 
##' @export
countFragmentOverlapsSecondCutter <- function(object, extend=TRUE, minMapq=0, shift=0){
  stopifnot(class(object)=="FourC")
  if(length(rowRanges(object)) == 0) 
    stop("Add fragments before calling 'findViewpointFragments'")
  
  cat("reading bam files\n")
  bamFiles = file.path(metadata(object)$bamFilePath, colData(object)$bamFile)
  colData(object)$originalReads = sapply(bamFiles, function(bamFile) countBam(bamFile)$records)
  reads = lapply(bamFiles, function(bamfile){  
    what <- c("mapq")
    flag <- scanBamFlag(isUnmappedQuery = FALSE)
    param <- ScanBamParam(what=what,flag = flag)
    
    readGAlignments(bamfile, param=param)
  })  
  colData(object)$rawReads = sapply(reads, length)
  
  ## filter bad quality reads
  if(minMapq >= 0){
    reads = lapply(reads, 
                   function(rga, minMapq){                 
                     if(!any(is.na(mcols(rga)$mapq))){
                       return(rga[mcols(rga)$mapq > minMapq])
                     }
                     warning("mapq quality filter could not be applied due to NA values in the mapq quality scores. All reads are used for counting overlaps.")
                     rga
                   },                     
                   minMapq=minMapq)
    colData(object)$lowQualityReads = colData(object)$rawReads - sapply(reads, length)
  }
  
  ## convert to GRanges to certify strand specific counting
  reads = lapply(reads, granges)
  
  ## trim reads if trim > 0
  if(extend){
    cutterLen = nchar(metadata(object)$reSequence2)
    reads = lapply(reads, function(gr, cutterLen){
      start(gr[strand(gr)=="+"]) = start(gr[strand(gr)=="+"]) - cutterLen
      end(gr[strand(gr)=="-"]) = end(gr[strand(gr)=="-"]) + cutterLen  
      gr  
    }, cutterLen=cutterLen)
  }
  
  cat("calculating overlaps\n")  
  frag <- rowRanges(object)
  
  ref = getReferenceSeq(object)
  site = getSites(metadata(object)$reSequence2, ref)
  
  ## left and right switch because the reads extended from the cutting site
  strand(site) <- "+"
  countsRightSiteEnd <- sapply(reads, countOverlaps, query=site, type=c("start"), maxgap=shift)
  
  strand(site) <- "-"
  countsLeftSiteEnd <- sapply(reads, countOverlaps, query=site, type=c("end"), maxgap=shift)
  
  ## map counts at second cutter sites to first cutter fragments
  ## multiple overlaps are not possible
  fragWithSite = findOverlaps(site, frag, minoverlap = nchar(metadata(object)$reSequence2), select="first")
  
  countsLeftFragmentEnd = matrix(0, nrow=length(frag), ncol=dim(countsRightSiteEnd)[2])
  countsRightFragmentEnd = matrix(0, nrow=length(frag), ncol=dim(countsRightSiteEnd)[2])  
  
  countsLeftFragmentEnd[unique(fragWithSite[!is.na(fragWithSite)]),] = apply(countsLeftSiteEnd, 2, function(siteCounts) tapply(siteCounts, fragWithSite, sum))
  countsRightFragmentEnd[unique(fragWithSite[!is.na(fragWithSite)]),] = apply(countsRightSiteEnd, 2, function(siteCounts) tapply(siteCounts, fragWithSite, sum)) 
    
  colData(object)$mappedReads = apply(countsLeftFragmentEnd, 2, sum) + apply(countsRightFragmentEnd, 2, sum)
  colData(object)$mappingRatio = colData(object)$mappedReads / (colData(object)$rawReads - colData(object)$lowQualityReads)
  
  metaDataFrame <- DataFrame(type=rep("countInfo", 5),
                             description=rep("", 5))
  idx <- colnames(colData(object)) %in% c("originalReads", 
                                          "rawReads", 
                                          "lowQualityReads",
                                          "mappedReads",
                                          "mappingRatio")
  mcols(colData(object))[idx,] <- metaDataFrame
  
  assays(object) <- SimpleList(countsLeftFragmentEnd=countsLeftFragmentEnd, 
                               countsRightFragmentEnd=countsRightFragmentEnd)
  object
}
