## Ovarian cancer datas
install.packages("readxl")

install.packages("ggplot2")
library(ggplot2)
library(readxl)
rna <- read_excel("C:/Users/Sindhuja Sridharan/Downloads/rna.xlsx.xlsx")
pnnlprot <- read_excel("C:/Users/Sindhuja Sridharan/Downloads/pnnlproteome.xlsx.xlsx")
jhuprot <- read_excel("C:/Users/Sindhuja Sridharan/Downloads/jhuproteome.xlsx.xlsx")

# PROTEINS transpose ---------------------------------------------------------------------------
tmp      <- as.matrix(jhuprot[,])
labels  <- tmp[,1]
jhuprot  <- as.matrix(jhuprot[,2:length(jhuprot)])
jhuprot_ <- t(jhuprot)
colnames(jhuprot_) <- labels

tmp      <- as.matrix(pnnlprot[,])
labels  <- tmp[,1]
pnnlprot  <- as.matrix(pnnlprot[,2:length(pnnlprot)])
pnnlprot_ <- t(pnnlprot)
colnames(pnnlprot_) <- labels

# RNA transpose
tmp      <- as.matrix(rna[,])
labels  <- tmp[,1]
rna  <- as.matrix(rna[,2:length(rna)])
rna_ <- t(rna)
colnames(rna_) <- labels


# JOIN proteome, use jhu for intersection

idx <- which(rownames(pnnlprot_) %in% rownames(jhuprot_))
pnnlprot_ <- pnnlprot_[-idx,]
protein <- rbind(jhuprot_,pnnlprot_)


# Find all intersections

idx <- which(rownames(protein) %in% rownames(rna_))
length(idx)  # only 105
protein<-protein[idx,]

idx <- which(rownames(rna_) %in% rownames(protein))
rna <- rna_[idx,]

##############################################################################################
#   REMOVE severely missing data
##############################################################################################

# Remove severely missing data >50% ?
x <- 52    #   can change the number (we have 105 observations)


# PROTEINS------------------------------------------------------

columns <- NULL
howmany <- NULL
for (j in 1:ncol(protein)){
  for (i in 1:nrow(protein)){
    if (is.na(protein[i,j])==TRUE){
      columns <- append(columns,i)
    }
  }
  howmany[j] <- length(columns)
  columns<-NULL
}

idx<- which((howmany>x)==TRUE)
protein<-protein[,-idx]

# RNA------------------------------------------------------

columns <- NULL
howmany <- NULL
for (j in 1:ncol(rna)){
  for (i in 1:nrow(rna)){
    if (is.na(rna[i,j])==TRUE){
      columns <- append(columns,i)
    }
  }
  howmany[j] <- length(columns)
  columns<-NULL
}

idx<- which((howmany>x)==TRUE)
if (length(idx)>0){
  rna<-rna[,-idx]
}

###########################################################################################

#### Mean imputation
##############################
### Protein
#############################

for (j in 1:ncol(protein)){
  mean_ <- mean(protein[,j], na.rm=TRUE)
  for (i in 1:nrow(protein)){
    if (is.na(protein[i,j])==TRUE){
      protein[i,j] <- mean_
    }
  }
}

#imputed_Data_prot <-ifelse(is.na(protein), mean(protein, na.rm=TRUE),protein)
#sum(is.na(imputed_Data_prot))

#############################
### RNA
#############################
#imputed_Data_rna <-ifelse(is.na(rna), mean(rna, na.rm=TRUE),rna)
#sum(is.na(imputed_Data_rna))

for (j in 1:ncol(rna)){
  mean_ <- mean(rna[,j], na.rm=TRUE)
  for (i in 1:nrow(rna)){
    if (is.na(rna[i,j])==TRUE){
      rna[i,j] <- mean_
    }
  }
}
#############################
###################
###Q-Q plots
###################

par(mfrow=c(1,2))
qqnorm(protein, main = "Normal Q-Q Plot",plot.it = TRUE)
qqnorm(rna, main = "Normal Q-Q Plot",plot.it = TRUE)

##############################
## Vizualize the imputed data
#############################
par(mfrow=c(1,2))
hist(protein)
hist(rna)


#################################
####Normalize data#######
##############################

#source("https://bioconductor.org/biocLite.R")
#biocLite("preprocessCore")
library(preprocessCore)
norm_prot <- normalize.quantiles(protein)
norm_rna <- normalize.quantiles(rna)

#####################################
# PUT COLNAMES AND ROWNAMES BACK ON

colnames(norm_prot)<-colnames(protein)
rownames(norm_prot)<-rownames(protein)

colnames(norm_rna)<-colnames(rna)
rownames(norm_rna)<-rownames(rna)

##############################
## Vizualize the normalized data
#############################

par(mfrow=c(1,2))
hist(norm_prot)
hist(norm_rna)

###################################
###PCR
###################################
require(pls)
set.seed (1000)
x_<- data.matrix(norm_rna)
y_ <- data.matrix(norm_prot[,1])

pcr_model <- pcr(y_~ x_, validation = "CV")
#summary(pcr_model)
#validationplot(pcr_model, val.type='MSEP')

mse_v<-numeric()
for(i in 1:10) {
  pcr.pred <-as.vector(predict(pcr_model, x_, ncomp = i))
  mse_v[i] <- mean((pcr.pred - y_[,1])^2)
}
mse_v

#visualize MSE
y_m <- mse_v
x_m <- 1:10
mse_p <- data.frame(x_m, y_m)

ggplot()+
  geom_point(data=mse_p, aes(x=x_m, y=y_m), size=2)+
  geom_line(data=mse_p, aes(x=x_m, y=y_m), size=1)

##############################################
############PLS
###############################################
set.seed (1500)
x_<- data.matrix(norm_rna)
y_ <- data.matrix(norm_prot[,1])

pls_model <- plsr(y_~ x_, validation = "CV")
#summary(pls_model)
#validationplot(pls_model, val.type='MSEP')


msepls_v<-numeric()
for(i in 1:10) {
  pls.pred <-as.vector(predict(pls_model, x_, ncomp = i))
  msepls_v[i] <- mean((pls.pred - y_[,1])^2)
}
msepls_v

#visualize MSE
y_m <- msepls_v
x_m <- 1:10
msepls_p <- data.frame(x_m, y_m)

ggplot()+
  geom_point(data=msepls_p, aes(x=x_m, y=y_m), size=2)+
  geom_line(data=msepls_p, aes(x=x_m, y=y_m), size=1)

#############################################
#shuffle independently each X matrix
###########################################
#########PLS
set.seed(3000)
sam <-data.frame(norm_rna)
samshuff <-sample(sam)
y_ <- data.matrix(norm_prot[,1])
z<-data.matrix(samshuff)
pls_modelshu <- plsr(y_~ z, validation = "CV")
#summary(pls_modelshu)

msepls1_v<-numeric()
for(i in 1:10) {
  pls.pred1 <-predict(pls_modelshu, z, ncomp = i)
  msepls1_v[i] <- mean((pls.pred1 - y_[,1])^2)
}
msepls1_v

#visualize MSE
y_m <- msepls1_v
x_m <- 1:10
msepls1_p <- data.frame(x_m, y_m)

ggplot()+
  geom_point(data=msepls1_p, aes(x=x_m, y=y_m), size=2)+
  geom_line(data=msepls1_p, aes(x=x_m, y=y_m), size=1)
##################PCR
set.seed(3500)

pcr_modelshuf <- pcr(y_~ z, validation = "CV")
#summary(pcr_model1)

mse1_v<-numeric()
for(i in 1:10) {
  pcr.pred1 <-as.vector(predict(pcr_modelshuf, z, ncomp = i))
  mse1_v[i] <- mean((pcr.pred1 - y_[,1])^2)
}
mse1_v

#visualize MSE
y_m <- mse1_v
x_m <- 1:10
mse1_p <- data.frame(x_m, y_m)

ggplot()+
  geom_point(data=mse1_p, aes(x=x_m, y=y_m), size=2)+
  geom_line(data=mse1_p, aes(x=x_m, y=y_m), size=1)


matplot(x_m, cbind(mse_v,msepls_v,mse1_v,msepls1_v))

######################################
###Analysing raw data
#####################################
# Histograms and density lines
###   RNA
par(mfrow=c(3, 3))
rnanames <- dimnames(rna)[[2]]
for (i in 2:10) {
  hist(rna[,i], main=rnanames[i], probability=TRUE, col="gray", border="white")
  d <- density(rna[,i])
  lines(d, col="red")
}

#####PROTEIN
par(mfrow=c(3, 3))
protnames <- dimnames(protein)[[2]]
for (i in 2:10) {
  hist(protein[,i], main=protnames[i], probability=TRUE, col="gray", border="white")
  d <- density(protein[,i], na.rm=TRUE)
  lines(d, col="red")
}

#####Mean and Variance

meanrna <- apply(rna, 2, mean)
varrna <- apply(rna, 2, var)
sdrna <- apply(rna, 2, sd)
range(varrna)
range(meanrna)
hist(meanrna, xlab="Histogram of Mean RNA") 
hist(varrna)
plot(sort(meanrna))
plot(sort(varrna))

d1 <- density(meanrna)
plot(d1)
q1<- qqnorm(meanrna, main = "meanrna")


meanprotein <- apply(protein, 2, mean,na.rm=TRUE)
varprotein <- apply(protein, 2, var,na.rm=TRUE)
sdprotein <- apply(protein, 2, sd,na.rm=TRUE)
range(varprotein)
range(meanprotein)
hist(meanprotein, xlab="Histogram of Mean RNA") 
hist(varprotein)
plot(sort(meanprotein))
plot(sort(varprotein))

d2 <- density(meanprotein,na.rm=TRUE)
plot(d2)
q2<- qqnorm(meanprotein, main = "meanprotein")


################## Low varince filter
newrna <- subset(varrna, varrna > .5)
plot(sort(newrna))
variancex <- data.frame(newrna)
rownames(variancex) == colnames(rna)
write.csv(variancex, file = "var1.csv")
write.csv(t(rna), file = "var2.csv")
nrnaa <- read_excel("C:/Users/Sindhuja Sridharan/Desktop/Projectdiff/nrnaa.xlsx",col_names = TRUE)
tmp1     <- as.matrix(nrnaa[,])
labels1  <- tmp1[,1]
nrnaa  <- as.matrix(nrnaa[,2:length(nrnaa)])
nrnaa <- t(nrnaa)
colnames(nrnaa) <- labels1

######

PCA1 <- prcomp(nrnaa, center = TRUE, scale. = TRUE)

summary(PCA1)

plot(PCA1, type = "b") ##Scree plot

PCAscores1 <- PCA1$x
PCAloadings1 <- PCA1$rotation

plot(PCAscores1[,1:5],  # x and y data
     pch=21,           # point shape
     cex=1.5, 
     col= c("red","blue","green","yellow","black"),
     main="Scores"     # title of plot
)
