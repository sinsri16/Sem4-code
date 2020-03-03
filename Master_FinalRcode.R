## Ovarian cancer datas
install.packages("readxl")
library(readxl)
install.packages("ggplot2")
library(ggplot2)
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
summary(pcr_model)
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

##Fold change
error_mse1<-numeric()
for(i in 1:10) {
  datamse1 <- data.matrix(mse_v)
  error_mse1[i] <- (datamse1[i]/datamse1[7])
}
error_mse1

##############################################
############PLS
###############################################
set.seed (1500)
x_<- data.matrix(norm_rna)
y_ <- data.matrix(norm_prot[,1])

pls_model <- plsr(y_~ x_, validation = "CV")
summary(pls_model)
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

##Fold change
error_msepls1<-numeric()
for(i in 1:10) {
  datamse2 <- data.matrix(msepls_v)
  error_msepls1[i] <- (datamse2[i]/datamse2[7])
}
error_msepls1

###############################################
####PCR AFTER SHUFFLE
################################################
require(pls)
set.seed(2000)
xshuffle <- x_[sample(nrow(x_)),]
y_ <- data.matrix(norm_prot[,1])

pcr_model1 <- pcr(y_~ xshuffle, validation = "CV")
summary(pcr_model1)
#validationplot(pcr_model1, val.type='MSEP')

mse1_v<-numeric()
for(i in 1:10) {
  pcr.pred1 <-as.vector(predict(pcr_model1, xshuffle, ncomp = i))
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


##Fold change
error_mse3<-numeric()
for(i in 1:10) {
  datamse3 <- data.matrix(mse1_v)
  error_mse3[i] <- datamse3[i]/datamse3[7]
}
error_mse3

#####################################################
#####PLS AFTER SHUFFLE
#####################################################
set.seed (2500)
xshuffle <- x_[sample(nrow(x_)),]
y_ <- data.matrix(norm_prot[,1])

pls_model1 <- plsr(y_~ xshuffle, validation = "CV")
summary(pls_model1)
#validationplot(pls_model1, val.type='MSEP')


msepls1_v<-numeric()
for(i in 1:10) {
  pls.pred1 <-predict(pls_model1, xshuffle, ncomp = i)
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

##Fold change
error_mse4<-numeric()
for(i in 1:10) {
  datamse4 <- data.matrix(msepls1_v)
  error_mse4[i] <- datamse4[i]/datamse4[7]
}
error_mse4

########################
#####Error rates
#######################
error_mse1
error_msepls1
error_mse3
error_mse4
################################################################
########Visualizing the mses together
#################################################################

#matplot(x_m, cbind(mse_v,msepls_v,mse1_v,msepls1_v))
matplot(x_m,mse_v,ylim = range(c(mse_v-error_mse1, mse_v+error_mse1)))
arrows(x_m,mse_v-error_mse1,x_m, mse_v+error_mse1,length=0.05, angle=90, code=3)
matpoints(x_m, mse_v, type = "p")
matlines (x_m, mse_v, type = "l")


matplot(x_m,msepls_v)
arrows(msepls_v-error_msepls1, x_m, msepls_v+error_msepls1, x_m, length=0.05, angle=90, code=3)
matplot(x_m, mse1_v)
arrows(mse1_v-error_mse3,x_m, mse1_v+error_mse3,x_m, length=0.05, angle=90, code=3)
matplot(x_m,msepls1_v)
arrows(msepls1_v-error_mse4,x_m, msepls1_v+error_mse4,x_m, length=0.05, angle=90, code=3)

legend("right", inset=.05,col = 1:4, legend=c("mse_v","msepls_v","mse1_v","msepls1_v"), pch=1, horiz=FALSE)
matpoints(x_m, cbind(mse_v,msepls_v,mse1_v,msepls1_v),col = 1:4, type = "p")
matlines (x_m, cbind(mse_v,msepls_v,mse1_v,msepls1_v),col = 1:4, type = "l")



