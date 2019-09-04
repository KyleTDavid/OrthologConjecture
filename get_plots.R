library(ggplot2)
library(gridExtra)
library(sfsmisc)
library(ape)

setwd("~/Desktop/ortholog/")

#branch_model.py output
df_raw <- read.table("masterdata_unfiltered.txt",header=F,sep='\t',row.names = NULL,stringsAsFactors = F)
names(df_raw) <- c("tree_num","tree_likelihood","node_type","MRCA","age","bl_1","dN_1","dS_1","w_1","leaves_1",
               "bl_2","dN_2","dS_2","w_2","leaves_2","delta_w")
outgroups <- c("Bilateria", "Caenorhabditis elegans", "Chordata", "Ciona", "Ciona intestinalis B CG-2006", "Ciona savignyi", "Drosophila melanogaster", "Ecdysozoa", "Opisthokonta",
               "Saccharomyces cerevisiae", "Vertebrata", "Petromyzon marinus")

#filter trees without duplications
df <- df_raw[df_raw$tree_num %in% unique(df_raw[df_raw$node_type=="D",]$tree_num),]
#filter trees without speciations
df <- df[df$tree_num %in% unique(df[df$node_type=="S",]$tree_num),]

df <- df[!(df$MRCA %in% outgroups),]

#filter synonymous subsitutions > 2 to avoid saturation issues
df <- df[df$dS_1<=2.0,]
df <- df[df$dS_2<=2.0,]

#filter synonymous substitutions <0.01, unreliable for dN/dS estimation
df <- df[df$dS_1>=0.01,]
df <- df[df$dS_2>=0.01,]

#filter dN/dS estimates > 10, too big to be believable/workable
df <- df[df$w_1<=10,]
df <- df[df$w_2<=10,]

#how many trees in filtered dataset?
length(unique(df$tree_num))
#how many species?
length(unique(df[grepl(" ", df$MRCA),]$MRCA))

dfS <- df[df$node_type=="S",]
#how many speciations?
length(dfS$tree_num)

dfD <- df[df$node_type=="D",]
#how many duplications?
length(dfD$tree_num)

#get null dataframe df with binomial variable bin
make_null <- function(df, bin) {
  nulldf <- df
  var1 <- sort(unique(nulldf[[bin]]))[1]
  var2 <- sort(unique(nulldf[[bin]]))[2]
  null_list <- sample((1:length(nulldf[,1])), size=length(nulldf[nulldf[[bin]]==var1,][,1]), replace = F)
  nulldf[[bin]] <- var2
  nulldf[[bin]] <- replace(x=nulldf[[bin]],list = null_list ,values=var1)
  return(nulldf)
}

run <- function(df, bin, var) {
  nullval <- list()
  bin1 <- sort(unique(df[[bin]]))[1]
  bin2 <- sort(unique(df[[bin]]))[2]
  for (i in 1:1000) {
    null <- lapply(split(df, f = df$tree_num), make_null, bin=bin)
    null <- do.call("rbind", null)
    nullval[i] <- abs(mean(null[null[[bin]]==bin1,][[var]])-
                        mean(null[null[[bin]]==bin2,][[var]]))
  }
  emp <- abs(mean(df[df[[bin]]==bin1,][[var]])-
        mean(df[df[[bin]]==bin2,][[var]]))
  return((sum(emp<unlist(nullval)))/1000)
}

run(df, "node_type", "delta_w")

abs(mean(df[df[["node_type"]]=="D",][["delta_w"]])-
      mean(df[df[["node_type"]]=="S",][["delta_w"]]))
sum((0.4123797<testnulls)/1000)

# 
# #simulation tests
# permt <- function(dfA,dfB,var) {
#   df1 <- data.frame(dfA)
#   df2 <- data.frame(dfB)
#   df1$count = 1
#   df2$count = 2
#   df_null <- data.frame(rbind(df1,df2))
#   df_null_vals <- list()
#   for (i in 1:1000){
#     df_null$count <- 1
#     null2s<-sample((1:length(df_null[,1])),
#                    size = length(df2[,1]),replace = FALSE)
#     df_null$count<-replace(x=df_null$count,list = null2s,values=2)
#     null1 <- df_null[df_null$count==1,]
#     null2 <- df_null[df_null$count==2,]
#     df_null_vals[i] <- as.numeric(abs(mean(null1[[var]])-mean(null2[[var]])))
#     df_null_vals<-unlist(df_null_vals)
#   }
#   df_p <- sum(abs(mean(dfA[[var]]) - mean(dfB[[var]]))
#               < df_null_vals)/1000
#   return(df_p)
# }

#Hedge's Gs
hg <- function(array1, array2) { return((mean(array2)-mean(array1))/(sqrt(((length(array1)*var(array1))+((length(array2)*var(array2))))
                                                                          /(length(array2)+length(array1)-2))))
}

#OVL
#integrate area under both density plots for overlap coefficient
OVL <- function(df) {
  a <- df[df$node_type=="S",]$delta_w
  b <- df[df$node_type=="D",]$delta_w
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
  total <- integrate.xy(d$x, d$a) + integrate.xy(d$x, d$b)
  intersection <- integrate.xy(d$x, d$w)
  # compute overlap coefficient
  OVL <- intersection / min(integrate.xy(d$x, d$a),integrate.xy(d$x, d$b))
  return(OVL)
}

#the good stuff
#permt(dfS,dfD,"delta_w")
hg(dfS$delta_w, dfD$delta_w)
OVL(df)
ks.test(dfS$delta_w, dfD$delta_w)
mean(dfS$delta_w)
sd(dfS$delta_w)
mean(dfD$delta_w)
sd(dfD$delta_w)

#aggregate by clade
for (clade in unique(df$MRCA)) {
  cat(clade,"\t",run(df[df$MRCA==clade,],"node_type", "delta_w"),"\n")
}

#same species paralogs
dfD$cat <- ifelse(grepl(" ", dfD$MRCA), "IN", "OUT")

dfin <- dfD[grepl(" ", dfD$MRCA),]
mean(dfin$delta_w)
sd(dfin$delta_w)
#different species paralogs
dfout <- dfD[!(grepl(" ", dfD$MRCA)),]
mean(dfout$delta_w)
sd(dfout$delta_w)

run(dfD,"cat","delta_w")
abs(mean(dfD[dfD[["cat"]]=="IN",][["delta_w"]])-
      mean(dfD[dfD[["cat"]]=="OUT",][["delta_w"]]))
sum((0.1410249<paralogcatnull)/1000)

hg(dfD[dfD$cat=="IN",]$delta_w, dfD[dfD$cat=="OUT",]$delta_w)


#top part of time plot with inparalogs removed
dfout <-df[!(df$node_name=="IN"),]
dfout$category <- factor(dfout$node_type, levels=c("S","D"))
ggplot(dfout,aes(-age,delta_w)) + scale_fill_manual(values = c(rgb(1,0,0,0.5),rgb(0,0,1,0.5))) + scale_color_manual(values = c('#7f84ff','#ff876e')) + geom_smooth(aes(group=category,colour=category,fill=category),linetype=0) + coord_cartesian(xlim=c(0,-4.5),ylim=c(0,0.8)) +
  theme(legend.position="none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.title = element_blank(),  axis.line = element_line(colour = "black")) 

dfin <-df[df$node_name=="IN",]
dfinag <- aggregate(delta_w ~ MRCA, dfin, FUN = function(x) c(mn = mean(x), n = length(x) ) )


dfw1 <- df[c("node_type","w_1")]
dfw2 <- df[c("node_type","w_2")]
names(dfw2) <- c("node_type","w_1")
dfw <- rbind(dfw1, dfw2)
names(dfw) <- c("node_type","w")
permt(dfw[dfw$node_type=="S",], dfw[dfw$node_type=="D",], "dwp")
hg(dfw[dfw$node_type=="S",]$w, dfw[dfw$node_type=="D",]$w)
mean(dfw[dfw$node_type=="S",]$w)
sd(dfw[dfw$node_type=="S",]$w)
mean(dfw[dfw$node_type=="D",]$w)
sd(dfw[dfw$node_type=="D",]$w)

#delta parent stuff
dfp = df
dfp$upleaves = paste(df$leaves_1, df$leaves_2, sep=", ")
dfp$parent1 = dfp$w_1[match(dfp$upleaves, dfp$leaves_1)]
dfp$parent2 = dfp$w_2[match(dfp$upleaves, dfp$leaves_2)]
dfp$parent = ifelse(is.na(dfp$parent1), dfp$parent2, dfp$parent1)
dfp = dfp[,!(names(dfp) %in% c("parent1", "parent2"))]
dfp = na.omit(dfp)
dfDp = dfp[dfp$node_type=="D",]
dfSp = dfp[dfp$node_type=="S",]
mean(c(dfSp$w_1-dfSp$parent,dfSp$w_2-dfSp$parent))
sd(c(dfSp$w_1-dfSp$parent,dfSp$w_2-dfSp$parent))
mean(c(dfDp$w_1-dfDp$parent,dfDp$w_2-dfDp$parent))
sd(c(dfDp$w_1-dfDp$parent,dfDp$w_2-dfDp$parent))
dfp$dwp1 <- dfp$w_1 - dfp$parent
dfp$dwp2 <- dfp$w_2 - dfp$parent

dfdwp1 <- data.frame("node_type"=dfp$node_type,"dwp"=apply(dfp[c("dwp1","dwp2")],1,min),"cat"="min","tree_num"=dfp$tree_num, stringsAsFactors = F)
dfdwp2 <- data.frame("node_type"=dfp$node_type,"dwp"=apply(dfp[c("dwp1","dwp2")],1,max),"cat"="max","tree_num"=dfp$tree_num, stringsAsFactors = F)
dfdwp <- rbind(dfdwp1, dfdwp2)

mean(dfdwp[dfdwp$node_type=="D",]$dwp)
sd(dfdwp[dfdwp$node_type=="D",]$dwp)


#permt(dfdwp[dfdwp$node_type=="S",], dfdwp[dfdwp$node_type=="D",], "dwp")
hg(dfdwp[dfdwp$node_type=="S",]$dwp, dfdwp[dfdwp$node_type=="D",]$dwp)
#permt(dfdwp[dfdwp$cat=="min" & dfdwp$node_type=="D",], dfdwp[dfdwp$cat=="min" & dfdwp$node_type=="S",], "dwp")
run(dfdwp1, "node_type", "dwp")
run(dfdwp2, "node_type", "dwp")



hg(dfdwp[dfdwp$cat=="min" & dfdwp$node_type=="D",]$dwp, dfdwp[dfdwp$cat=="min" & dfdwp$node_type=="S",]$dwp)
#permt(dfdwp[dfdwp$cat=="max" & dfdwp$node_type=="D",], dfdwp[dfdwp$cat=="max" & dfdwp$node_type=="S",], "dwp")
hg(dfdwp[dfdwp$cat=="max" & dfdwp$node_type=="D",]$dwp, dfdwp[dfdwp$cat=="max" & dfdwp$node_type=="S",]$dwp)


mean(dfD[dfD$age>=0.1 & dfD$age<=0.4,]$delta_w)
sd(dfD[dfD$age>=0.1 & dfD$age<=0.4,]$delta_w)

mean(dfout[dfout$age>=0.1 & dfout$age<=0.4,]$delta_w)
sd(dfout[dfout$age>=0.1 & dfout$age<=0.4,]$delta_w)

mean(dfdwp[dfdwp$cat=="max" & dfdwp$node_type=="S",]$dwp)/mean(dfdwp[dfdwp$cat=="min" & dfdwp$node_type=="S",]$dwp)
mean(dfdwp[dfdwp$cat=="max" & dfdwp$node_type=="D",]$dwp)/mean(dfdwp[dfdwp$cat=="min" & dfdwp$node_type=="D",]$dwp)

#fig 2 density plot 
SS <- density(log(dfS$delta_w))
DS <- density(log(dfD$delta_w))
#plot speciation curve and fill in
plot(SS,cex.lab=1.2,ylim=c(0,0.5), xlim=c(-10,4),xlab='',bty='n',main='',ylab='',xaxt='n',yaxt='n')
polygon(SS$x,SS$y,col=rgb(1,0,0,0.5),lty = 0)
#plot duplication curve and fill in
par(new=T)
plot(cex.lab=1.2,DS,ylim=c(0,0.5), xlim=c(-10,4),xlab='',ylab='',xaxt='n',yaxt='n',bty='n',main='',labels=F)
polygon(DS$x,DS$y,col=rgb(0,0,1,0.5),lty = 0)
#x-axis
axis(1,at=c(log(10),log(1),log(0.1),log(0.01),log(0.001),log(0.0001)),labels=F)
#y-axis
axis(2,at=seq(0,0.5,by=0.1),labels=F)

#fig2 box
p1 <- ggplot(aes(node_type,dwp, fill=node_type), data=dfdwp[dfdwp$cat=="min",]) + geom_boxplot(outlier.shape = NA) + coord_cartesian(ylim = c(-1.5,2.0)) + 
  theme(legend.position="none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.title = element_blank(),  axis.line = element_line(colour = "black")) 

p2 <- ggplot(aes(node_type,dwp, fill=node_type), data=dfdwp[dfdwp$cat=="max",]) + geom_boxplot(outlier.shape = NA) + coord_cartesian(ylim = c(-1.5,2.0)) + 
  theme(legend.position="none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.title = element_blank(),  axis.line = element_line(colour = "black")) 

grid.arrange(p1,p2,nrow=1)

#fig3 gam
df$category <- factor(df$node_type, levels=c("S","D"))
ggplot(df,aes(-age,delta_w)) + scale_fill_manual(values = c(rgb(1,0,0,0.5),rgb(0,0,1,0.5))) + scale_color_manual(values = c('#7f84ff','#ff876e')) + geom_smooth(aes(group=category,colour=category,fill=category),linetype=0) + coord_cartesian(xlim=c(0,-4.5),ylim=c(0,0.8)) +
  theme(legend.position="none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.title = element_blank(),  axis.line = element_line(colour = "black")) 
#fig3 tree
tree <- read.tree(file="4lab/Supplemental/Supplemental_V.nh.txt")
plot(ladderize(tree),show.tip.label = F)

#fig S2
dfS2 <- df[!(grepl(" ", df$MRCA)),]
dfS2$category <- factor(dfS2$node_type, levels=c("S","D"))
ggplot(dfS2,aes(-age,delta_w)) + scale_fill_manual(values = c(rgb(1,0,0,0.5),rgb(0,0,1,0.5))) + scale_color_manual(values = c('#7f84ff','#ff876e')) + geom_smooth(aes(group=category,colour=category,fill=category),linetype=0) + coord_cartesian(xlim=c(0,-4.5),ylim=c(0,0.8)) +
  theme(legend.position="none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.title = element_blank(),  axis.line = element_line(colour = "black")) 

#fig S3
ONE <- df[c("node_type","bl_1","dN_1","dS_1","w_1")]
names(ONE) <- c("node_type","bl","dN","dS","w")
TWO <- df[c("node_type","bl_2","dN_2","dS_2","w_2")]
names(TWO) <- c("node_type","bl","dN","dS","w")
dNdS <- rbind(ONE,TWO)
dNdS$category <- factor(dNdS$node_type, levels=c("S","D"))
p1 <- ggplot(dNdS,aes(-bl,dN)) + scale_fill_manual(values = c(rgb(1,0,0,0.5),rgb(0,0,1,0.5))) + scale_color_manual(values = c('#7f84ff','#ff876e')) + geom_smooth(aes(group=category,colour=category,fill=category),linetype=0) + 
  theme(legend.position="none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.title = element_blank(),  axis.line = element_line(colour = "black")) 
p2 <- ggplot(dNdS,aes(-bl,dS)) + scale_fill_manual(values = c(rgb(1,0,0,0.5),rgb(0,0,1,0.5))) + scale_color_manual(values = c('#7f84ff','#ff876e')) + geom_smooth(aes(group=category,colour=category,fill=category),linetype=0) + 
  theme(legend.position="none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.title = element_blank(),  axis.line = element_line(colour = "black")) 
grid.arrange(p1,p2,nrow=2)






#do trees with more duplication nodes have higher delta omega?
fraction_duplications <- function(df) {
  return(length(df[df$node_type=="D",]$node_type)/length(df$node_type))
}

test1 <- lapply(split(df, f = df$tree_num), fraction_duplications)
test1 <- do.call("rbind",testout)
test1 <- data.frame(tree_num=rownames(test1), percentage=test1)

test2 <- aggregate(delta_w ~ tree_num, df[df$node_type=="S",], mean)

merge <- merge(test1, test2)

plot(merge$percentage, merge$delta_w)

summary(lm(merge$percentage ~ merge$delta_w))






