library(xlsx)
    summarize_cell_counts <- function(pattern){
        filenames=list.files(path=".", full.names=TRUE,pattern = pattern)
        datalist = lapply(filenames, function(x){read.csv(file=x,header=T)})
        mat <- Reduce(function(x,y) {merge(x,y, all = TRUE)}, datalist)
        mat[is.na(mat)] <- 0
        return(mat)
    }

ann_names =  c("Blueprint_Encode","HPCA","Mouse-RNAseq","Immgen");

annotations = paste("*_down_cell_ann_",ann_names,".csv",sep = "")

cell_counts  <- lapply(annotations, FUN =  function(x){summarize_cell_counts(pattern = x)} )

