
summarize_cell_counts <- function(pattern){
    filenames=list.files(path=".", full.names=TRUE,pattern = pattern)
    message(filenames)
    datalist = lapply(filenames, function(x){read.csv(file=x,header=T)})
    mat <- Reduce(function(x,y) {merge(x,y, all = FALSE,by="Name")}, datalist)
                  mat[is.na(mat)] <- 0
                  return(mat)
                 }

    dat = summarize_cell_counts(pattern = "*scores.csv")


    colnames(dat) = gsub("Score_model_at","",colnames(dat))
      colnames(dat) = gsub(".h5","",colnames(dat))

      
