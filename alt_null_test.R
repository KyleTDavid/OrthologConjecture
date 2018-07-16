args=commandArgs(trailingOnly = TRUE)

#read in and subset data
df <- read.table(file =args[1],header=F,sep='\t',row.names = NULL)
alt <- data.frame(df)
null <- data.frame(df)
S <- subset(x = df,subset = df$V3=='S')
D <- subset(x = df,subset = df$V3=='D')

#integrate area under both density plots for overlap coefficient
a <- S$V16
b <- D$V16

# define limits of a common grid, adding a buffer so that tails aren't cut off
lower <- min(c(a, b))-1
upper <- max(c(a, b))+1

# generate kernel densities
da <- density(a, from=lower, to=upper)
db <- density(b, from=lower, to=upper)
d <- data.frame(x=da$x, a=da$y, b=db$y)

# calculate intersection densities
d$w <- pmin(d$a, d$b)

# integrate areas under curves
library(sfsmisc)
total <- integrate.xy(d$x, d$a) + integrate.xy(d$x, d$b)
intersection <- integrate.xy(d$x, d$w)

# compute overlap coefficient
OVL <- intersection / min(integrate.xy(d$x, d$a),integrate.xy(d$x, d$b))

#mean divergence (all nodes)
DW <- mean(df$V16)

#mean divergence (speciation)
DWS <- mean(S$V16)

#get Hedges' D
HD <- (mean(D$V16)-mean(S$V16))/(sqrt(((length(S$V16)*var(S$V16))+((length(D$V16)*var(D$V16))))/(length(D$V16)+length(S$V16)-2)))

# count Ds and Ss in empirical data
empD<-sum(df$V3=="D")
empS<-sum(df$V3=="S")

# get diff between D and S
EC <- abs(mean(D$V16)-mean(S$V16))

# generate null values
nullist <- list()
for (i in 1:1000){
  null$V3 <- "S"
  nullDs<-sample((1:length(df$V1)), size = empD,replace = FALSE)
  null$V3<-replace(x=null$V3,list = nullDs,values="D")
  nullS <- subset(x = null,subset = null$V3=='S')
  nullD <- subset(x = null,subset = null$V3=='D')
  nullist[i] <- as.numeric(abs(mean(nullS$V16)-mean(nullD$V16)))
  nullist<-unlist(nullist)
}

# Empirical P
Ep <- (sum(EC < nullist))/1000

altlist <- list()
for (i in 1:1000){
  alt$V3 <- "S"
  dm<-exp(df$V16-mean(df$V16))
  altD<-sample((1:length(df$V1)), size = empD,replace = FALSE,prob = dm)
  alt$V3<-replace(x=alt$V3,list = altD,values="D")
  altS <- subset(x = alt,subset = alt$V3=='S')
  altD <- subset(x = alt,subset = alt$V3=='D')
  altlist[i] <- as.numeric(abs(mean(altS$V16)-mean(altD$V16)))
  altlist<-unlist(altlist)
}

# Alt P
Ap <- (sum(mean(altlist) < nullist))/1000

DS = (mean(D$V16)-mean(S$V16))

cat(Ep,Ap,OVL,HD,empD,empS,"\n")
