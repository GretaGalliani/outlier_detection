k=3
S=c(1,1,3,2,3,1,2)
empty_list <- vector(mode="list",length=k)

for (j in 1:k){
    empty_list[j] <- vector(mode="list",length=1)
}

for(j in 1:k)
{  index <- which(S==j)
      empty_list[j]<- c(index)
}
