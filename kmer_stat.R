library(ggplot2)

kmer <- read.table("/Users/wangpeng/Desktop/Kmer.out", comment.char = "",
                   header = F, sep = "\t")

sort.kmer <- sort(kmer$V1, decreasing = T)

stat <- as.data.frame(table(sort.kmer))

findpeak <- function(x){
  pks <- c()
  signs <- sign(diff(x, na.pad = F))
  #1, -1 or 1, 0 is the peak point
  for(ind in 1:length(signs)){
    if(signs[ind] > 0 && signs[ind+1] < 0)
      pks <- c(pks, ind+1) 
  }
  return(pks)
}


#kmer depth
pp <- findpeak(stat$Freq)[1]
#genome size  kmer_nummber / kmer_depth
gsize <- sum(as.numeric(stat[,1]) * stat[,2]) / pp
gsize/1000/1000


poisdtb <- dpois(1:200, pp) *  gsize 

plot(poisdtb, type='l', lty=1, col="green")

lines(stat[2:200,], type="l", xlim = c(2,200))
