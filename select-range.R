args <- commandArgs(trailingOnly = TRUE)

read.metrics <- function(fname){
  cmd <- paste("sed -n '/^insert_size/,$p' ", fname)
  mf <- pipe(cmd)
  read.table(mf, header=TRUE)
}

select.valid_range <- function(h){
  m<-max(h[,2])
  mode=mean(h[h[,2]==m,1])
  h1<-min(h[h[,2]>=m/2,1])
  h2<-max(h[h[,2]>=m/2,1])
  t1<-min(h[h[,2]>=m/10,1])
  t2<-max(h[h[,2]>=m/10,1])
  c(max(c(t1, mode-(mode-h1)*2)), max(c(t2, mode+(h2-mode)*2)))
}

h<- read.metrics(args[1])
r<-select.valid_range(h)
cat(r[1], r[2],"\n")

