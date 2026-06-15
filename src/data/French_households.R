library(dplyr)
hh = read.csv(file="../data/raw/TD_MEN4_2021.csv",header = TRUE, sep=';')
 
t =table(hh$NPERC)
t=t[-1]
prob.size = t/sum(t)
S=rep(0,100000)
for(n in 1:100000){
  a= which(runif(n=1) <cumsum(prob.size))[1]
  
  S[n] = as.numeric(a)+1
  
}
saveRDS(S, file='../data/processed_data/INSEE_French_HH_size.rds')
