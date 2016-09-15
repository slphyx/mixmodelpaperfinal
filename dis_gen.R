three <- NA
n <- 1000

d1prop <- .45
d2prop <- .35
d3prop <- .2
prop <- c(d1prop,d2prop,d3prop)

d1mu <- 3
d2mu <- 5
d3mu <- 6.5
mus <- c(d1mu,d2mu,d3mu)

d1sd <- 1.2
d2sd <- 1.3
d3sd <- 1.1
sds <- c(d1sd,d2sd,d3sd)

d <- list()
for(i in 1:length(prop)){
  d[[i]] <- rlnorm(n*prop[i],log(mus[i]),log(sds[i]))  
}

three <- unlist(d)


write.table(round(three,2),"three.csv", col.names = FALSE, row.names = FALSE)
