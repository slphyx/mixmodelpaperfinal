three <- NA
n <- 2000

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


four <- NA


d1prop <- .3
d2prop <- .3
d3prop <- .2
prop <- c(d1prop,d2prop,d3prop,.1)

d1mu <- 3
d2mu <- 5
d3mu <- 6.5
mus <- c(d1mu,d2mu,d3mu, 8.5)

d1sd <- 1.2
d2sd <- 1.3
d3sd <- 1.1
sds <- c(d1sd,d2sd,d3sd, 1.5)

d <- list()
for(i in 1:length(prop)){
  d[[i]] <- rlnorm(n*prop[i],log(mus[i]),log(sds[i]))  
}

four <- unlist(d)


write.table(round(four,2),"four.csv", col.names = FALSE, row.names = FALSE)


five <- NA


d1prop <- .2
d2prop <- .2
d3prop <- .2
prop <- c(d1prop,d2prop,d3prop, .2,.2)

d1mu <- 2
d2mu <- 4.5
d3mu <- 7
mus <- c(d1mu,d2mu,d3mu, 9.5, 15)

d1sd <- 1.2
d2sd <- 1.1
d3sd <- 1.1
sds <- c(d1sd,d2sd,d3sd, 1.1, 1.1)

d <- list()
for(i in 1:length(prop)){
  d[[i]] <- rlnorm(n*prop[i],log(mus[i]),log(sds[i]))  
}

five <- unlist(d)


write.table(round(five,2),"five.csv", col.names = FALSE, row.names = FALSE)