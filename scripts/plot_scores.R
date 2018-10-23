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

pattern <- opt$input
outname <- opt$out

library(ggplot2)
library(reshape2)

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

summarize_cell_counts <- function(pattern){
    basename = basename(pattern);
    dirname  = dirname(pattern);
    filenames=list.files(path=dirname, full.names=TRUE,pattern = basename)

    filenames = grep(".scores.csv",filenames,value=TRUE)
    message(filenames)
    datalist = lapply(filenames, function(x){read.csv(file=x,header=T)})
    mat <- Reduce(function(x,y) {merge(x,y, all = FALSE,by="Name")}, datalist)
    mat[is.na(mat)] <- 0
    return(mat)
}

dat = summarize_cell_counts(pattern = pattern)



colnames(dat) = gsub("*.h5","",colnames(dat))

colnames(dat) = sapply(colnames(dat), function(x) substrRight(x, 5))
colnames(dat)[1] = "0";
colnames(dat) = as.numeric(colnames(dat))


x = melt(dat)
p <- ggplot(x,aes(x = variable, y = value)) + geom_boxplot()
name <- paste("Summary of run:",basename(pattern));
p  <-  p + labs(title = name)
p  <-  p + xlab("Iteration")
p  <-  p + ylab(expression(paste(log[2],frac(P(x/H),P(x/R)), sep = "")))
p <- p + theme(axis.text.x = element_text(angle = 90, size=4,hjust = 1))

ggsave(outname,p, dpi= 300,units="cm", height = 12,width=24,limitsize = TRUE )


