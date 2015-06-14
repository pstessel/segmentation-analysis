setwd("/Volumes/HD2/Users/pstessel/Documents/Github/R/450_R_Code")
load("apphappyData.RData")
q24 <- read.csv(file="q24.csv", header = FALSE, sep = ",", quote = "\"",
         dec = ".", fill = TRUE)
ls()
require(cluster)
require(useful)
require(Hmisc)
numdata <- apphappy.3.num.frame
labdata <- apphappy.3.labs.frame
backup <- numdata

### Imputing Data

random.imp <- function (a){
  missing <- is.na(a)
  n.missing <- sum(missing)
  a.obs <- a[!missing]
  imputed <- a
  imputed[missing] <- sample (a.obs, n.missing, replace=TRUE)
  return (imputed)
}

q5r1.imp <- random.imp(numdata$q5r1)
q12.imp <- random.imp(numdata$q12)
q57.imp <- random.imp(numdata$q57)

### numdata.fixed <- data.frame(numdata, q5r1.imp, q12.imp, q57.imp)

### numdata.fixed$q5r1 <- NULL
### numdata.fixed$q12 <- NULL
### numdata.fixed$q57 <- NULL

### Following are EDA

str(numdata)
summary(numdata)
a=table(numdata$q1)
a
barplot(a)
b=table(numdata$q1,numdata$q2r1)
b
barplot(b)
hist(numdata$q1)

a <- table(numdata$q26r17)
a
barplot(a)


numsub <- subset(numdata, select=c("q24r1","q24r2","q24r3","q24r2r3", "q25r1","q25r2","q25r3",
                                   "q26r3","q26r4","q26r5"))
rcorr(as.matrix(numsub), type="pearson")
                 
rcorr(as.matrix(numsub), type="pearson")

###
numdata$avg_q24r2r3 <- (numdata$q24r2 + numdata$q24r3)/2
###

### numsub <- subset(numdata, select=c("q24r1","avg_q24r2r3", "q24r4", "q24r5", "q24r6", "q24r7", "q24r8",
                                   "q24r9", "q24r10", "q24r11", "q24r12"))
### rcorr(as.matrix(numsub), type="pearson")

numsub <- subset(numdata, select=c("q24r1",  "q24r2",	"q24r3",	"q24r4",	"q24r5",	"q24r6",	"q24r7",
                                    "q24r8",	"q24r9",	"q24r10",	"q24r11",	"q24r12"))
rcorr(as.matrix(numsub), type="pearson")


###
numdata$avg_q25r1r2r3 <- (numdata$q25r1 + numdata$q25r2 + numdata$q25r3)/3
###



str(numsub)
summary(numsub)

require(corrplot)
numsubcorrelation <- cor(numsub)
corrplot(numsubcorrelation)

mcor <- cor(numsub)
corrplot(mcor, method="shade", shade.col=NA, tl.col="black")

#corrplot(numsubcorrelation, method="shade", addCoef.col="black", #addCoefasPercent=TRUE ,type="lower", shade.col=NA, tl.col="black", #tl.srt=45, addcolorlabel="no", order="AOE",insig = "p-value")





### Create a 'scree' plot to determine the num of clusters

wssplot <- function(numsub, nc=15, seed=1234) {
  wss <- (nrow(numsub)-1)*sum(apply(numsub,2,var))
  for (i in 2:nc) {
    set.seed(seed)
    wss[i] <- sum(kmeans(numsub, centers=i)$withinss)}
  plot(1:nc, wss, type="b", xlab="Number of Clusters",
       ylab="Within groups sum of squares")} 

wssplot(numsub)

### Create a Kmeans with 5 clusters

clusterresults <- kmeans(numsub,3)
clusterresults
clusterresults$withinss
clusterresults$tot.withinss
clusterresults$totss
clusterresults$betweenss
clusterresults$size

### Create a PC plot

plot(clusterresults, data=numsub)

clusterresults$centers

names(clusterresults)
head(clusterresults$cluster)

newdf <- as.data.frame(clusterresults$cluster)

write.csv(newdf, file = "clusterresults.csv")

write.csv(numsub, file = "numsub.csv")

### PCA Analysis to 'see' the 'clusters'

fit <- princomp(numsub, cor=TRUE)
summary(fit)
loadings(fit)
plot(fit,type="lines") 
fit$scores
biplot(fit)
str(fit$scores)
score <- data.frame(fit$scores)
names(score)
write.csv(score, file = "score.csv")
newdf <- read.csv("clusterresults.csv")
newdf$clus <- newdf$clusterresults.cluster
PC_Clusterdata <- cbind(score,newdf)
head(score)
head(newdf)
head(PC_Clusterdata)
PC_Clusterdata$clus <- factor(PC_Clusterdata$clus)
head(PC_Clusterdata)

g1<-ggplot(PC_Clusterdata, aes(x = Comp.1, y = Comp.2)) + geom_point(aes(color=factor(clus)))

g1

g2<-ggplot(PC_Clusterdata, aes(x = Comp.3, y = Comp.4)) + geom_point(aes(color=factor(clus)))

g2

###################
my_hist3d <- function(x, y, freq=FALSE, nclass="auto") {
  n<-length(x)
  if (nclass == "auto") { nclass<-ceiling(sqrt(nclass.Sturges(x))) }
  breaks.x <- seq(min(x),max(x),length=(nclass+1))
  breaks.y <- seq(min(y),max(y),length=(nclass+1))
  h <- NULL
  for (i in 1:nclass) 
    for (j in 1:nclass) 
      h <- c(h, sum(x <= breaks.x[j+1] & x >= breaks.x[j] & y <= breaks.y[i+1] & y >= breaks.y[i] ) )
  if (freq) h <- h / n
  xx <- as.factor(round(mean(breaks.x[1:2])+(0:(nclass-1))*diff(breaks.x[1:2]), 1))
  yy <- as.factor(round(mean(breaks.y[1:2])+(0:(nclass-1))*diff(breaks.y[1:2]), 1))
  res <- cbind(expand.grid(xx,yy), h)
  colnames(res) <- c(deparse(substitute(x)),deparse(substitute(y)),'Frequency')
  formu <- as.formula(paste("Frequency ~ ", paste(colnames(res)[1:2], collapse= "+")))
  cloud(formu, res, panel.3d.cloud=panel.3dbars, col.facet='lightblue', 
        xbase=1, ybase=1, scales=list(arrows=FALSE, col=1), 
        par.settings = list(axis.line = list(col = "transparent")))
}

library(latticeExtra)

############################

s <- read.csv("score.csv")
x <- s$Comp.1
y <- s$Comp.2
my_hist3d(x, y, nclass=10)




### Create a MDS to 'see' the 'clusters'

d <- dist(numsub)
fit <- cmdscale(d,eig=TRUE, k=2)
names(fit)
fit

x <- fit$points[,1]
y <- fit$points[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", 
     main="Metric MDS", type="n")
text(x, y, labels = row.names(numsub), cex=.7)

mdsscore <- data.frame(fit$points)
write.csv(mdsscore, file = "mdsscore.csv")
mds <- read.csv("mdsscore.csv")
MDS_Clusterdata <- cbind(mds,newdf)
head(MDS_Clusterdata)
MDS_Clusterdata$clus <- factor(MDS_Clusterdata$clus)

g1<-ggplot(MDS_Clusterdata, aes(x = X1, y = X2)) + geom_point(aes(color=factor(clus)))
g1


x <- mds$X1
y <- mds$X2
my_hist3d(x, y, nclass=10)



### Create a dataset with the original data with the cluster info
### This will be useful for creating profiles for the clusters



newdf <- read.csv("clusterresults.csv")
combdata <- cbind(numsub,newdf,numdata$q1)
head(combdata)

require(reshape)
combdata <- rename(combdata, c(clusterresults.cluster="cluster"))
head(combdata)


aggregate(combdata,by=list(byvar=combdata$cluster), mean)

### My Code ###################################################

is.na(x) # returns TRUE of x is missing
y <- c(1,2,3,NA)
is.na(y) # returns a vector (F F F T)

### Imputing Data

random.imp <- function (a){
  missing <- is.na(a)
  n.missing <- sum(missing)
  a.obs <- a[!missing]
  imputed <- a
  imputed[missing] <- sample (a.obs, n.missing, replace=TRUE)
  return (imputed)
}

q5r1.imp <- random.imp(numdata$q5r1)
q12.imp <- random.imp(numdata$q12)
q57.imp <- random.imp(numdata$q57)

numdata.fixed <- data.frame(numdata, q5r1.imp, q12.imp, q57.imp)

missing <- is.na(numdata$q5r1)

numdata.fixed <- numdata.fixed[ -c(q12, q57) ]

numdata.fixed$q5r1

ndf_labels <- names(numdata.fixed)

### Output to an Excel Spreadsheet
library(xlsx)
write.xlsx(ndf_labels, "/label_names.csv")

library(xlsx)
write.xlsx(numdata.fixed, file="numdata", sheetName="Sheet1")

# Write to a file, suppress row names
write.csv(numdata.fixed, "data.csv", row.names=FALSE)

# Same, except that instead of "NA", output blank cells
write.csv(combdata, "cluster.csv", row.names=FALSE, na="")

# Use tabs, suppress row names and column names
write.table(data, "data.csv", sep="\t", row.names=FALSE, col.names=FALSE) 

hcresults <- hclust (d=dist(numsub))

plot(hcresults)
