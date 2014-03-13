# depends  
# library(Rsamtools)
# library(plyr)
# library(Rgraphviz)

#' @title Class ContigLinkViewer
#' @description Class for storing descriptive statistics on the pairs of sequence which connect the contigs as well as to estimate assembly quality.
#' @details An assembly experiment from short pairs of sequences rarely gives a unique contig. This class aims (i) at providing
#' some descriptive statistics relating to the alignment of paires on the contigs (ii) at bringing out links between contig. 
#' In this class, the information given from pairs mapping on the same contig or from pairs binding two differents contig is treated differently.
#' Insert size estimated from sequencing pairs mapping on the same contig is used as scaffolding parameter.   
#' @slots ContigLength A named vector giving the length of each contig
#' @slots dfSame A data.frame in sam like format containing the pairs of sequence that have been aligned on the same contig
#' @slots dfDiff A data.frame in sam like format containing the pairs of sequence that have been aligned on different contig
#' @slots ExpInsertSize A numeric vector for the experimental insert size which should have been given with the sequencing protocol
#' @slots InsertSizeSame A numeric vector given descriptive statistics from the inserts Size of the pairs of sequence that have been aligned on the same contig
#' @slots AssQual A list of 4 data.frame. See \code{\link{EstimateAssQual}} method. 
#' @slots ContigLinkCount A list of 2 Objects. See \code{\link{ContigLinkCount}} method. 
#' @author Delphine Charif & Joseph Tran
#' @rdname ContigLinkViewer-class
#' @name ContigLinkViewer-class
#' @exportClass ContigLinkViewer
setClass("ContigLinkViewer",
         representation(
           ContigLength ="vector",
           dfSame="data.frame",
           dfDiff="data.frame",
           ExpInsertSize="numeric",
           InsertSizeSame="numeric",
           InsertSizeDiff="numeric",  
           AssQual="list",
           ContigLinkCount="list"
         )
)

setGeneric(
  name="EstimateInsertSizeFromSame",
  def = function( object, ... ){standardGeneric("EstimateInsertSizeFromSame")}
)

setGeneric(
  name="EstimateInsertSizeFromDiff",
  def = function( object, ... ){standardGeneric("EstimateInsertSizeFromDiff")}
)

setGeneric(
  name="ContigLinkCount",
  def = function( object, MinCount ){standardGeneric("ContigLinkCount")}
)


setGeneric(
  name="EstimateAssQual",
  def = function( object, ... ){standardGeneric("EstimateAssQual")}
)

setGeneric(
  name="SetAssQual<-",
  def = function( object, value ){standardGeneric("SetAssQual<-")}
)

setGeneric(
  name="GetAssQual",
  def = function( object){standardGeneric("GetAssQual")}
)


setGeneric(
  name="SetInsertSizeSame<-",
  def = function( object, value ){standardGeneric("SetInsertSizeSame<-")}
)
setGeneric(
  name="GetInsertSizeSame",
  def = function( object){standardGeneric("GetInsertSizeSame")}
)
setGeneric(
  name="SetInsertSizeDiff<-",
  def = function( object, value ){standardGeneric("SetInsertSizeDiff<-")}
)

setGeneric(
  name="GetInsertSizeDiff",
  def = function( object){standardGeneric("GetInsertSizeDiff")}
)


setGeneric(
  name="SetContigLinkCount<-",
  def = function( object,value ){standardGeneric("SetContigLinkCount<-")}
)

#' @title Setteur method for the slot ContigLinkCount of the class \code{\link{ContigLinkViewer}}.
#' @rdname ContigLinkViewer-class
#' @name SetContigLinkCount
setReplaceMethod(f="SetContigLinkCount", 
                 signature=c("ContigLinkViewer", "list"),
                 definition=function(object,value){ 
                   object@ContigLinkCount <- value
                   return(object) 
                 } ) 

#' @title Setteur method for the slot InsertSizeSame of the class \code{\link{ContigLinkViewer}}.
#' @rdname ContigLinkViewer-class
#' @name SetInsertSizeSame
setReplaceMethod(f="SetInsertSizeSame", 
                 signature=c("ContigLinkViewer","numeric"), 
                 definition=function(object, value){
                   object@dfSame$Isize <- value
                   slot(object,"InsertSizeSame") <- c(summary(value),"sd"=sd(value,na.rm=TRUE))
                   return(object) 
                 } ) 

#' @title Setteur method for the slot InsertSizeDiff of the class \code{\link{ContigLinkViewer}}.
#' @rdname ContigLinkViewer-class
#' @name SetInsertSizeDiff
setReplaceMethod(f="SetInsertSizeDiff", 
                 signature=c("ContigLinkViewer","numeric"),
                 definition=function(object, value){
                   object@dfDiff$Isize <- value
                   slot(object,"InsertSizeDiff") <- c(summary(value),"sd"=sd(value,na.rm=TRUE))
                   return(object) 
                 } )
#' @title Setteur method for the slot AssQual of the class \code{\link{ContigLinkViewer}}.
#' @rdname ContigLinkViewer-class
#' @name SetAssQual
setReplaceMethod(f="SetAssQual", 
                 signature=c("ContigLinkViewer", "list"),
                 definition=function(object,value){ 
                   slot(object,"AssQual") <- value
                   return(object) 
                 }) 
#' @title Getteur method for the slot AssQual of the class \code{\link{ContigLinkViewer}}.
#' @rdname ContigLinkViewer-class
#' @name GetAssQual
setMethod(f="GetAssQual", 
          signature=c("ContigLinkViewer"),
          definition=function(object){ 
            slot(object,"AssQual")
          }) 

#' @title Getteur method for the slot InsertSizeSame of the class \code{\link{ContigLinkViewer}}.
#' @rdname ContigLinkViewer-class
#' @name GetInsertSizeSame
setMethod(f="GetInsertSizeSame",
          signature="ContigLinkViewer",
          definition=function(object){
            slot(object,"InsertSizeSame")
          })

#' @title Getteur method for the slot InsertSizediff of the class \code{\link{ContigLinkViewer}}.
#' @rdname ContigLinkViewer-class
#' @name GetInsertSizeDiff
setMethod(f="GetInsertSizeDiff",
          signature="ContigLinkViewer",
          definition=function(object){
            slot(object,"InsertSizeDiff")
          })

#' @title Constructor method for the class \code{\link{ContigLinkViewer}}.
#' @details This method instanciates the slots dfSame, dfDiff and ExpInsertSize. For pairs mapping on 
#' different contigs, stored in dfDiff, the Ori1 and Ori2 columns are added. These columns indicates the 
#' contigs orientation so as to giving the good orientation (+ - or -> <-) to reads which connect them.    
#' @param BamFile A Bam file resulting from the alignment of paired-end sequences or
#' mated paire sequences on a set of contigs. The Bam file has to be sorted and indexed. 
#' The sequences that do not mapped or which are singletons have to be filtered out.  
#' @param ExperimentalInsertSize The experimental insert size. The mean length of dna fragment that 
#' have been keeped for the sequencing experiment.
#' @author Delphine Charif & Joseph Tran
#' @rdname ContigLinkViewer-class  
#' @name ContigLinkViewer-class
setMethod("initialize",
          signature(.Object="ContigLinkViewer"),
          function(.Object,BamFile, ExperimentalInsertSize){
            
            # Scan Bam File
            param <- ScanBamParam(tag=c("NM", "NH"), what=c("flag","isize"))
            x <- readGAlignmentsListFromBam(BamFile, use.names=TRUE, param=param,asMates=TRUE)
            
            # Get Reference sequence length
            .Object@ContigLength <- seqlengths(x)
            
            # From the List of mates => Create a dataframe with one row per mate
            # From GAlignment to data.frame
            dfx <- as.data.frame(x)
            seqx <- seq(1,dim(dfx)[1],2)
            df <- merge( dfx[seqx,], dfx[seqx+1,],by.x="element",by.y="element",suffixes=c(".1",".2"))
            
            # Split the data.frame into two:
            # -> dfSame: Mate mapped on same Ref
            # -> dfDiff: Mate mapped on diff Ref
            df$isSameRef <- df$seqnames.1 == df$seqnames.2
            df.split <- split(df,df$isSameRef)
            .Object@dfDiff <- df.split[[1]]
            .Object@dfSame <- df.split[[2]]
            
            # Add contig Orientation
            
            .Object@dfDiff$Ori1 <- rep("F",dim(.Object@dfDiff)[1])
            .Object@dfDiff$Ori2 <- rep("F",dim(.Object@dfDiff)[1])
            
            tmpNum <- bamFlagTest( .Object@dfDiff$flag.1,'isFirstMateRead')
            tmp <- paste( .Object@dfDiff$strand.1, .Object@dfDiff$strand.2,sep="")
            .Object@dfDiff$Ori1[which(tmp=="-+" & tmpNum==TRUE)]<-"R"
            .Object@dfDiff$Ori2[which(tmp=="-+" & tmpNum==TRUE)]<-"R"
            .Object@dfDiff$Ori2[which(tmp=="++" & tmpNum==TRUE)]<-"R"
            .Object@dfDiff$Ori1[which(tmp=="++" & tmpNum==FALSE)]<-"R"
            .Object@dfDiff$Ori1[which(tmp=="--" & tmpNum==TRUE)]<-"R"
            .Object@dfDiff$Ori2[which(tmp=="--" & tmpNum==FALSE)]<-"R"
            
            # Set the ExperimentalInsertSize from argument
            .Object@ExpInsertSize = ExperimentalInsertSize
            
            return(.Object)
          }
)

#' @title Estimate the insert size of a library from pairs of sequence which align
#' on the same contig.
#' @details Estimate the insert size of a library from pairs of sequence which align
#' on the same contig and which have the good orientation. The well mapped paired-end sequences
#' (-> <-) belong to this two categories:
#' \itemize{
#' \item{1} The strand of the first sequence in pair is + and its sam flag is 99. 
#' The strand of the second sequence in the pair is - and its flag is 147.
#' The insert size is obtained by substracting the start position of the first sequence on the 
#' contig to the end position of the second sequence on the contig and by adding 1.
#' \item{2} The strand of the first sequence in pair is - and its sam flag is 83. 
#' The strand of the second sequence in the pair is + and its flag is 163.
#' The insert size is obtained by substracting the start position of the second sequence on the 
#' contig to the end position of the first sequence on the contig and by adding 1.
#' }
#' @param object An object of class \code{\link{ContigLinkViewer}}
#' @value An object of class \code{\link{ContigLinkViewer}} with the slot InsertSizeSame.
#' @author Delphine Charif & Joseph Tran
#' @rdname ContigLinkViewer-method
#' @name EstimateInsertSizeFromSame               
setMethod(
  f="EstimateInsertSizeFromSame",
  signature="ContigLinkViewer",
  definition = function( object, ... ){
    
    Isize <- rep(NA,dim(object@dfSame)[1])
    PropPaired1 <- which( object@dfSame$strand.1=="+" & object@dfSame$strand.2=="-" & object@dfSame$flag.1 == 99 & object@dfSame$flag.2 == 147)
    PropPaired2 <- which( object@dfSame$strand.1=="-" & object@dfSame$strand.2=="+" & object@dfSame$flag.1 == 83 & object@dfSame$flag.2 == 163)
    Isize[PropPaired1] <- object@dfSame[PropPaired1,]$end.2 - object@dfSame[PropPaired1,]$start.1 + 1
    Isize[PropPaired2] <- object@dfSame[PropPaired2,]$end.1 - object@dfSame[PropPaired2,]$start.2 + 1
    
    # Echo Warning when the Estimate Isize is +- 1 sd far from the experimental ISize  
    if( object@ExpInsertSize <= (mean(Isize,na.rm=TRUE)-sd(Isize,na.rm=TRUE))  | object@ExpInsertSize >= (mean(Isize,na.rm=TRUE)+sd(Isize,na.rm=TRUE))){
      warning("The expected insert size (experimental) is outside the confidence interval of the insert size estimated from mates mapping on the same reference")
    }
    print(c(summary(Isize),"sd"=sd(Isize,na.rm=TRUE)))
    # Update Object
    SetInsertSizeSame(object) <- Isize
    
    return(object)
  }
)


#' @title  Estimate the insert size of a library from pairs of sequence which align
#' on different contigs.
#' @details  Estimate the insert size of a library from pairs of sequence which align
#' on two contigs. In this case the orientation of the sequences on the contigs does not matter. 
#' There is two ways to obtain the insert size.
#' \itemize{
#'  \item{The sequence is the first in the pair} The insert size is the length of the first contig plus
#'  the end position of the second sequence on the second contig minus the start position of the first 
#'  sequence on the first contig minus 1.
#'  \item{The sequence is the second in the pair} The insert size is the length of the first contig plus
#'  the end position of the first sequence on the first contig minus the start position of the second 
#'  sequence on the second contig minus 1.
#' }
#' @param An object of class \code{\link{ContigLinkViewer}} 
#' @value An object of class \code{\link{ContigLinkViewer}} with the slot InsertSizeDiff.
#' @author Delphine Charif & Joseph Tran
#' @rdname ContigLinkViewer-method
#' @name EstimateInsertSizeFromDiff
setMethod(
  f="EstimateInsertSizeFromDiff",
  signature="ContigLinkViewer",
  definition = function( object, ... ){
    
    tmp <- bamFlagTest(object@dfDiff$flag.1,'isFirstMateRead')
    Isize <- ifelse(tmp,
                    object@ContigLength[object@dfDiff$seqnames.1] + object@dfDiff$end.2 - object@dfDiff$start.1 - 1,
                    object@ContigLength[object@dfDiff$seqnames.1] + object@dfDiff$end.1 - object@dfDiff$start.2 - 1)
    
    print(c(summary(Isize),"sd"=sd(Isize,na.rm=TRUE)))
    # Update Object
    SetInsertSizeDiff(object) <- Isize
    
    return(object)
  }
)

#' @title Estimate the assembly quality
#' @details This method creates the AssQual slot which contain 4 objects:
#' \itemize{
#'  \item{parameters} A vector giving the insert size parameters which will be used to estimate the assembly quality. Min and Max or Mean and sd.
#'  \item{dfSameQualAssByContig} A data.frame summarizing the repartition of the pairs according to the insert size range and the orientation of the pairs.
#'  This data.frame has 6 columns:
#'    \itemize{
#'      \item{seqnames.1} A vector giving the name of the contig
#'      \item{strand.1} A vector giving the strand of the first sequence in pair
#'      \item{strand.2} A vector giving the strand of the second sequence in pair
#'      \item{outof} A vector giving the number of pair that are outside the confidence interval for the insert size
#'      \item{nonav} A vector giving the number of pair for which the insert size can not be computed
#'      \item{within} A vector giving the number of pair with an insert size inside the confidence interval
#'      }
#'  \item{dfSameQualAssTot} A data.frame which sum the dfSameQualAssByContig by orientation only disregarding contig information.
#'  \item{dfDiffQualAssByContig} A data.frame summarizing the repartition of the pairs according to the insert size range and the orientation of the contigs
#'  that are linked. 
#'  This data.frame has 6 columns:
#'  for the insert size
#'      \item{seqnames.1} A vector giving the name of the first contig
#'      \item{seqnames.2} A vector giving the name of the second contig
#'      \item{Ori1}  A vector giving the orientation of the first contig
#'      \item{Ori2}  A vector giving the orientation of the second contig
#'      \item{outof} A vector giving the number of pair that are outside the confidence interval for the insert size
#'      \item{within} A vector giving the number of pair with an insert size inside the confidence interval
#' } 
#' To get the quality, a confident interval has to be computed for the insert size.
#' @param object An object of class \code{\link{ContigLinkViewer}}  
#' @param EstimateCutOff If TRUE, the interval for the insert size is [MeanIsize-(SdNumber*SdIsize), MeanIsize+(SdNumber*SdIsize)]
#' @param Min If EstimateCutOff is FALSE The min value for the insert Size
#' @param Max If EstimateCutOff is FALSE The max value for the insert Size
#' @param SdNumber The number of sd to keep for computing the insert size interval.
#' @value An object of class \code{\link{ContigLinkViewer}} with the AssQual slot
#' @author Delphine Charif & Joseph Tran
#' @rdname ContigLinkViewer-method
#' @name EstimateAssQual
setMethod(
  f="EstimateAssQual",
  signature="ContigLinkViewer",
  definition = function( object,
                         EstimateCutOff = TRUE,
                         Min=NA,
                         Max=NA,
                         SdNumber=NA ){
    
    if(EstimateCutOff == TRUE){
      if(! is.na(SdNumber)){
        MeanIsize=object@InsertSizeSame["Mean"] 
        SdIsize=object@InsertSizeSame["sd"]
        Min <- MeanIsize-(SdNumber*SdIsize)
        Max <- MeanIsize+(SdNumber*SdIsize)
        Parameters <- c("Min"=Min,"Max"=Max,"MeanSize"=MeanIsize,"SdIsize"=SdIsize,"SdNumber"=SdNumber)
      }
      else{
        stop("SdNumber is missing")
      }
    }
    else{
      if(any(Max & Min)){
        Parameters <- c("Min"=Min,"Max"=Max)
      }
      else{
        stop("Min or/and Max is/are missing")
      }
    }
    
    
    dfQualAssByContigSame <- as.data.frame(ddply( object@dfSame, .(seqnames.1,strand.1,strand.2), function(x){
      tmp <- rep(NA,length(x$Isize))
      tmp[which(x$Isize >= Min & x$Isize <= Max)]<-"within"
      tmp[which(x$Isize < Min | x$Isize > Max)]<-"outof"
      tmp[is.na(tmp)]<-"nonav"
      tmp <- factor(tmp,levels=c("outof","nonav","within"))
      table(tmp)
    }))
    
    dfQualAssTotSame <- as.data.frame(ddply(object@dfSame, .(strand.1,strand.2), function(x){
      tmp <- rep(NA,length(x$Isize))
      tmp[which(x$Isize >= Min & x$Isize <= Max)]<-"within"
      tmp[which(x$Isize < Min | x$Isize > Max)]<-"outof"
      tmp[is.na(tmp)]<-"nonav"
      tmp <- factor(tmp,levels=c("outof","nonav","within"))
      table(tmp)
    }))
    
    dfQualAssByContigDiff <- as.data.frame(ddply(object@dfDiff, .(seqnames.1,seqnames.2,Ori1,Ori2), function(x){
      tmp <- rep(NA,length(x$Isize))
      tmp[which(x$Isize >= Min & x$Isize <= Max)]<-"within"
      tmp[which(x$Isize <  Min | x$Isize > Max)]<-"outof"
      tmp <- factor(tmp,levels=c("outof","within"))
      table(tmp)
    }))
    
    SetAssQual(object) <- list("parameters"=Parameters,
                               "dfSameQualAssByContig" = dfQualAssByContigSame,
                               "dfSameQualAssTot"= dfQualAssTotSame,
                               "dfDiffQualAssByContig"=dfQualAssByContigDiff)
    
    return(object)
  }
)

#' @title Count the number of pairs of sequence that connect two contigs.
#' @details This method creates the CLC slot which is a list of 2 objects:
#' \itemize{
#' \item{MinCount} A vector giving the minimum number of pairs which bind two contigs
#' \item{CLC} A data.frame with 5 columns:
#'  \itemize{
#'  \item{seqnames.1} A vector giving the name of the first contig 
#'  \item{seqnames.2} A vector giving the name of the second contig
#'  \item{Ori1}  A vector giving the orientation of the first contig
#'  \item{Ori2}  A vector giving the orientation of the second contig
#'  \item{within} A vector giving the number of pairs above the MinCount parameter that connect two contigs
#' }
#' @param An object of class \code{\link{ContigLinkViewer}}
#' @param MinCount The minimum number of pairs which connect two contigs
#' @value An object of class \code{\link{ContigLinkViewer}} with the ContigLinkCount slot
#' @author Delphine Charif & Joseph Tran
#' @name ContigLinkCount
#' @rdname ContigLinkViewer-method
setMethod(f="ContigLinkCount",
          signature="ContigLinkViewer",
          definition = function(object, MinCount){
            # Test if there is insertSize within
            CLC <- object@AssQual$dfDiffQualAssByContig
            if(! any(which(CLC$within >= MinCount))){
              stop("Zero value upper MinCount")
            }
            CLC <- CLC[which(CLC$within >= MinCount),c(1:4,6)]
            ListCLC <- list("MinCount"=MinCount, "CLC"=CLC)
            SetContigLinkCount(object)<- ListCLC
            return(object)
          })

#' @title Plot method for the class \code{\link{ContigLinkViewer}}
#' @details This method allows to generate a graph to visualize the contigs and
#' their link using the Rgraphviz package. 
#' Contigs are the nodes and the number of pairs supporting the link between two are specified
#' by the legend of the edges.
#' @param object An object of class \code{\link{ContigLinkViewer}}
#' @author Delphine Charif & Joseph Tran
#' @name plot 
#' @rdname ContigLinkViewer-method
setMethod(f="plot",
          signature="ContigLinkViewer",
          definition = function(object, x=NA, ...){
            CLC <- object@ContigLinkCount$CLC
            
            CLC$name1 <- paste(CLC$seqnames.1, CLC$Ori1,sep="_")
            CLC$name2 <- paste(CLC$seqnames.2, CLC$Ori2,sep="_")
            
            graph <- new("graphNEL",nodes=sort(unique(c(CLC$name1,CLC$name2))),edgemode="directed")
            
            
            for(i in 1:dim(CLC)[1]){
              graph <- addEdge(CLC[i,6],CLC[i,7],graph,CLC[i,5])
            }
            
            # Graph Attributes
            
            # Edge Attr
            
            eAttr <- list()
            ew <- as.character(unlist(edgeWeights(graph)))
            ew <- ew[setdiff(seq(along=ew), removedEdges(graph))]
            names(ew) <- edgeNames(graph)
            
            eAttr$label <- ew
            
            # Node Attr
            
            nAttrs <- list()
            col <- rep(object@ContigLength,2)
            names(col) <- c(paste(names(object@ContigLength),"_","F",sep=""),paste(names(object@ContigLength),"_","R",sep=""))
            col <- round(col/10000,0)
            cc <- colorRampPalette(c("blue", "orange", "red"))(max(col)+1)
            nodecol <- cc[col+1]
            names(nodecol)<-names(col)
            
            nAttrs$color <- nodecol
            
            # General Attributes
            attrs=list()
            attrs$edge$fontsize <- "22"
            attrs$node$fontsize <- "18"
            attrs$node$width <- "2"
            attrs$edge$arrowsize <- "0.5"
            
            plot(graph, edgeAttrs =eAttr , nodeAttrs = nAttrs, attrs=attrs)
          })

