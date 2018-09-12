#!/usr/bin/env Rscript
library(optparse) 
sessionInfo()
option_list = list(
    make_option(c("-i", "--input"),
                type="character",
                default=NULL,
                help="input - pattern of score file.", 
                metavar="character"),
    make_option(c("-o", "--out"),
                type="character",
                default=NULL,
                help="name of output plot.", 
                metavar="character")

); 

opt_parser <- OptionParser(option_list=option_list,
                           description = "\nRuns emptydroplets on raw cellranger output.",
                           epilogue = "Example:\n\n  Blah  \n\n");
opt <- parse_args(opt_parser);

if (is.null(opt$input)){
       print_help(opt_parser)
       stop("Missing input pattern!\n", call.=FALSE)
   }
   if (is.null(opt$out)){
          prinpt_help(opt_parser)
          stop("Missing output name!\n", call.=FALSE)
      }

path <- opt$input
name <- opt$name
  
library(ggplot2)
library(reshape2)
      summarize_cell_counts <- function(pattern){
          filenames=list.files(path=".", full.names=TRUE,pattern = pattern)
          message(filenames)
          datalist = lapply(filenames, function(x){read.csv(file=x,header=T)})
          mat <- Reduce(function(x,y) {merge(x,y, all = FALSE,by="Name")}, datalist)
          mat[is.na(mat)] <- 0
          return(mat)
      }

      dat = summarize_cell_counts(pattern = "*scores.csv")


      colnames(dat) = gsub("Score_.home.user.tmp.Standard_Challenge_ATGACCTCTATTACC_mis_0.fatest_","",colnames(dat))
      colnames(dat) = gsub("*.h5","",colnames(dat))



      x = melt(dat)

      p <- ggplot(x,aes(x = variable, y = value)) + geom_boxplot()


      p
