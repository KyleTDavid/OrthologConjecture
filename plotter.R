#first arg is masterdata file, second arg is output filename

#read in and subset data for graphing
args=commandArgs(trailingOnly = TRUE)
tmp <- read.table(file="masterdata.1.8.18_realnodes.txt",header=F,sep='\t',row.names = NULL,stringsAsFactors = F)
df <- tmp[tmp$V17=="AV",]

S <- subset(x = df,subset = df$V3=='S')
D <- subset(x = df,subset = df$V3=='D')

DS <- density(log(D$V16))
SS <- density(log(S$V16))

alt <- data.frame(df)
null <- data.frame(df)

empD<-sum(df$V3=="D")
empS<-sum(df$V3=="S")

# get diff between D and S
EC <- (mean(D$V16)-mean(S$V16))

# generate null values
nullist <- list()
for (i in 1:1000){
  null$V3 <- "S"
  nullDs<-sample((1:length(df$V1)), size = empD,replace = FALSE)
  null$V3<-replace(x=null$V3,list = nullDs,values="D")
  nullS <- subset(x = null,subset = null$V3=='S')
  nullD <- subset(x = null,subset = null$V3=='D')
  nullist[i] <- as.numeric(mean(nullD$V16)-mean(nullS$V16))
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
  altlist[i] <- as.numeric(mean(altD$V16)-mean(altS$V16))
  altlist<-unlist(altlist)
}

# Alt P
Ap <- (sum(mean(altlist) < nullist))/1000


png(filename=args[2],width = 300,height=250,units='mm',res=300)

#set layout
layout(matrix(c(1,1,1,1,1,1,1,1,2,2,2,2), 3, 4.9, byrow = TRUE))

#plot speciation curve and fill in
plot(SS,cex.lab=1.2,ylim=c(0,0.5), xlim=c(-10,4),xlab='',bty='n',main='',ylab='',labels=F)
polygon(SS$x,SS$y,col=rgb(1,0,0,0.5),lty = 0)

#plot duplication curve and fill in
par(new=T)
plot(cex.lab=1.2,DS,ylim=c(0,0.5), xlim=c(-10,4),xlab='',ylab='',xaxt='n',yaxt='n',bty='n',main='',labels=F)
polygon(DS$x,DS$y,col=rgb(0,0,1,0.5),lty = 0)

#legend
#legend("topright",c("Speciation","Duplication"),col=c(rgb(1,0,0,0.5),rgb(0,0,1,0.5)),pch=20,bty="n",cex=1.5)

nulldensity<-density(nullist)
altdensity<-density(altlist)
plot(0,xlim=c((min(EC,(altdensity$x),nulldensity$x)),(max(EC,(altdensity$x),nulldensity$x))),ylim=c(0,max(altdensity$y,nulldensity$y)),xlab='',ylab='')
polygon(nulldensity$x,nulldensity$y,col="gray",lty = 0)
polygon(altdensity$x,altdensity$y,col="gray17",lty = 0)


abline(v=EC,lwd=3,col="indianred")


dev.off()








