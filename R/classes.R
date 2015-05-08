##' FourC-class
##' 
##' @name FourC
##' @aliases FourC-class
##' @docType class
##' 
##' @note The \code{FourC} object extends the \code{DESeqDataSet} class. 
##' 
##' @param colData Column data that contains the required information for each library to set up
##' the \code{FourC} object:
##' \enumerate{
##' \item viewpoint, name of the viewpoint 
##' \item condition, experimental condition
##' \item replicate, replicate number
##' \item bamFile, file name of the bam file
##' \item sequencingPrimer, was the 4C library sequenced from the side of the first restriction enzyme cutting site or second
##' }
##' 
##' @param metadata Experimental data information required for the \code{FourC} object: 
##' \enumerate{
##' \item projectPath, directory where the project will be saved
##' \item fragmentDir, directory in the project directory where the information
##' about restriction fragments will be saved
##' \item referenceGenomeFile, path to the reference genome or a \code{BSgenome} object.
##' \item reSequence1, restriction enzyme recognition pattern of the first 
##' restriction enzyme
##' \item reSequence2, restriction enzyme recognition pattern of the second 
##' restriction enzyme
##' \item primerFile, path to the file containing the primer sequences used for
##' preparing the 4C libraries
##' \item bamFilePath, path to the directory where the bam files are stored
##' }
##' 
##' @keywords package
##' @rdname FourC
##' @export
setClass("FourC",
         contains = "DESeqDataSet")

setValidity( "FourC", function( object ) {
  if (! ("projectPath" %in% names(metadata(object))))
    return( "metadata must contain 'projectPath'" ) 
  if (! ("fragmentDir" %in% names(metadata(object))))
    return( "metadata must contain 'fragmentDir'" ) 
  if (! ("referenceGenomeFile" %in% names(metadata(object))))
    return( "metadata must contain 'referenceGenomeFile'" ) 
  if (! ("reSequence1" %in% names(metadata(object))))
    return( "metadata must contain 'reSequence1'" ) 
  if (! ("reSequence2" %in% names(metadata(object))))
    return( "metadata must contain 'reSequence2'" ) 
  if (! ("primerFile" %in% names(metadata(object))))
    return( "metadata must contain 'primerFile'" ) 
  if (! ("bamFilePath" %in% names(metadata(object))))
    return( "metadata must contain 'bamFilePath'" ) 
  if (! ("bamFile" %in% names(colData(object))))
    return( "colData must contain 'bamFile'" ) 
  if (! ("sequencingPrimer" %in% names(colData(object))))
    return( "colData must contain 'sequencingPrimer'" )
  TRUE
} )

#' @rdname FourC
#' @export 
setMethod("updateObject", "FourC",
    function(object, ..., verbose=FALSE)
    {
        ans <- callNextMethod()
        attr(ans, "vst") <- attr(object, "vst", exact=TRUE)
        attr(ans, "inverse-vst") <- attr(object, "inverse-vst", exact=TRUE)
        ans
    }
)

##' @rdname FourC
##' 
##' @import GenomicRanges SummarizedExperiment splines DESeq2 Biobase
##' @import Biostrings Rsamtools ggbio 
##' @import reshape2 fda GenomicAlignments Matrix rtracklayer LSD
##' 
##' @importFrom gtools combinations
##'  
##' 
##' @examples
##' 
##' 
##' metadata <- list(projectPath=tempdir(),
##'                  fragmentDir="re_fragments",
##'                  referenceGenomeFile=system.file("extdata/dm3_2L_1-6900.fa", 
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
##' @export
FourC <- function(colData, metadata){
  if (! ("viewpoint" %in% names(colData)))
    return( "colData must contain 'viewpoint'" ) 
  if (! ("condition" %in% names(colData)))
    return( "colData must contain 'condition'" ) 
  if (! ("replicate" %in% names(colData)))
    return( "colData must contain 'replicate'" ) 
  row.names(colData) <- with(as.data.frame(colData), paste(viewpoint, condition, replicate, sep="_"))
  colData$dummy = factor(paste0("dummy", rep_len(c(1,2), nrow(colData))))
  tmp <- DESeqDataSet(SummarizedExperiment(colData=colData, 
                                           metadata=metadata,
                                           assays=SimpleList(counts=matrix(seq_len(nrow(colData)),ncol=nrow(colData)))),
                      design=formula(~dummy))
  design(tmp) = formula(~1)
  colData(tmp) = colData(tmp)[,colnames(colData(tmp)) != "dummy"]
  fc <- new('FourC', 
            tmp)
  rm(tmp)
  
  ## remove dummy assay and row data again after validity of DESeqDataSet is passed
  assays(fc) <- SimpleList()
  rowRanges(fc) <- GRanges()
  
  fc
}
