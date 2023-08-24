#' @title Pattern Recognition Analysis (PRA)
#' @description
#' This function carries out the pattern recognition analysis
#' @param dataURA output URA function
#' @param BPname BPname
#' @param thres.role thres.role
#' @importFrom randomForest randomForest
#' @importFrom utils data
#' @return returns list of TSGs and OCGs when biological processes are provided, otherwise a randomForest based classifier that can be used on new data
#' @export
#' @examples
#' data(dataURA)
#' dataPRA <- PRA(dataURA = dataURA,
#' BPname = c("apoptosis","proliferation of cells"),
#' thres.role = 0)

PRA <- function(dataURA, BPname, thres.role = 0){

    # Check user input

    data(DiseaseList)

    if (!is.null(BPname) && all(BPname %in% names(DiseaseList)) == FALSE) {
	stop("BPname should be NULL or a character vector containing one or more BP(s) 
         among possible BPs stored in the DiseaseList object.")
    }
  
    if (is.null(dim(dataURA))) {
	stop("The URA data must be non-empty with genes in rows and BPs in columns")
    }
  
    if (all(colnames(dataURA) %in% names(DiseaseList)) == FALSE) {
	stop("The columns in the URA data must be BPs among possible BPs 
         stored in the DiseaseList")
    }
  
    if (!is.numeric(thres.role)) {
	stop("Thres.role must be numeric")
    }

    data('tabGrowBlock')
    data('knownDriverGenes')

    tabGrowBlock <- get("tabGrowBlock")
    knownDriverGenes <- get("knownDriverGenes")

    names.blocking <- tabGrowBlock[which(tabGrowBlock$Cancer.blocking == "Increasing"), "Disease"]
    names.growing <- tabGrowBlock[which(tabGrowBlock$Cancer.growing == "Increasing"), "Disease"]

    if(is.null(BPname)){

      	# print("random forest")
        common.genes.tsg <- intersect(knownDriverGenes$TSG, rownames(dataURA))
        common.genes.ocg <- intersect(knownDriverGenes$OCG, rownames(dataURA))

        dataTrain <- data.frame(dataURA[c(common.genes.tsg, common.genes.ocg),],
            "labels"=as.factor(c(rep(1,length(common.genes.tsg)),rep(0,length(common.genes.ocg)))))

        fit.rf <- randomForest((labels)~.,data=dataTrain, importance=TRUE)
        return(fit.rf)

    }else{
    	# print("selected processes")
        ind.proc.growing <- which(BPname %in% names.growing)
        ind.proc.blocking <- which(BPname %in% names.blocking)

        names.genes.tsg <- names(which(dataURA[,BPname[ind.proc.growing]] < -thres.role & dataURA[,BPname[ind.proc.blocking]] > thres.role))
        names.genes.ocg <- names(which(dataURA[,BPname[ind.proc.growing]] > thres.role & dataURA[,BPname[ind.proc.blocking]] < -thres.role))

        if(length(names.genes.tsg)>0){
            scores.tsg <- apply(abs(dataURA[names.genes.tsg,,drop=FALSE]), 1, mean)
        }else{
            scores.tsg <- NULL
        }

        if(length(names.genes.ocg)>0){
            scores.ocg <- apply(abs(dataURA[names.genes.ocg,,drop=FALSE]), 1, mean)
        }else{
            scores.ocg <- 0
        }

        return(list("TSG"=scores.tsg, "OCG"=scores.ocg))
    	# dataURA[,BPname]

    }
}

utils::globalVariables(c("DiseaseList"))
