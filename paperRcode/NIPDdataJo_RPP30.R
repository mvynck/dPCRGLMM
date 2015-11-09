require(ggplot2)
require(lme4)
data<-read.csv("~/Dropbox/Work/Matthijs-Olivier/dPCR/CNVReplicates/Data/NIPD/NIPDmerged.csv",sep=";")
sample.ID<-data[seq(1,42,3),1]
data<-data[,-1]

genenames<-unlist(strsplit(colnames(data),split="\\."))[seq(1,56,4)]
CNV<-array(0,length(genenames))
CNV.low<-array(0,length(genenames))
CNV.high<-array(0,length(genenames))
for(k in 1:14){
CNV.calc<-data.frame(genenames,CNV,CNV.low,CNV.high)
for(j in 1:14){
print(paste("Patient ",k,", gene ",j,sep=""))
flush.console()
target <- data[(c(k-1)*3+1):(c(k-1)*3+3),c((j-1)*2+1,(j-1)*2+2)]
reference <- data[(c(k-1)*3+1):(c(k-1)*3+3),27:28]
target <- target[complete.cases(target),]
reference <- reference[complete.cases(reference),]

samples<-rep(letters[1],sum(target[1,]))
for(i in 2:nrow(target)){
	samples<-c(samples,rep(letters[i],sum(target[i,])))
}

Y<-rep(1,as.vector(t(as.matrix(target))[1]))
for(i in 2:length(as.vector(t(as.matrix(target))))){
	Y<-c(Y,rep(i%%2,as.vector(t(as.matrix(target)))[i]))
}

targets<-rep(1,length(Y))

db <- data.frame(samples,Y,targets)

samples<-rep(letters[1],sum(reference[1,]))
for(i in 2:nrow(reference)){
	samples<-c(samples,rep(letters[i],sum(reference[i,])))
}

Y<-rep(1,as.vector(t(as.matrix(reference))[1]))
for(i in 2:length(as.vector(t(as.matrix(reference))))){
	Y<-c(Y,rep(i%%2,as.vector(t(as.matrix(reference)))[i]))
}

targets<-rep(0,length(Y))

db2 <- data.frame(samples,Y,targets)
db.final <- rbind(db,db2)
if(sum(target[,1]>50)){
model.fit <- glmer(Y~targets+(1|samples),data=db.final,family=binomial(cloglog),nAGQ=1)
(CNV.calc[j,2]<-exp(model.fit@beta[2])*2)
(CNV.calc[j,3]<-(exp(model.fit@beta[2])-1.96*sqrt(vcov(model.fit)[2,2]))*2)
(CNV.calc[j,4]<-(exp(model.fit@beta[2])+1.96*sqrt(vcov(model.fit)[2,2]))*2)
}
}
setwd("~/Dropbox/Work/Matthijs-Olivier/dPCR/CNVReplicates/Data/NIPD/results")
write.csv(CNV.calc,paste(sample.ID[k],".csv",sep=""))
}

setwd("~/Dropbox/Work/Matthijs-Olivier/dPCR/CNVReplicates/Data/NIPD/results/")
files<-list.files()
for(k in 1:14){
setwd("~/Dropbox/Work/Matthijs-Olivier/dPCR/CNVReplicates/Data/NIPD/results/")

CNV.calc<-read.csv(files[k])[-14,]
CNV.calc<-CNV.calc[,-1]

setwd("~/Dropbox/Work/Matthijs-Olivier/dPCR/CNVReplicates/Data/NIPD/figures")
png(paste(unlist(strsplit(files[k],"\\."))[1],".png",sep=""),width=1000,height=600)
q<-ggplot(CNV.calc, aes(x = CNV.calc[1:13,1], y = CNV.calc[1:13,2])) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymax = CNV.calc[1:13,4], ymin = CNV.calc[1:13,3])) +
  theme(axis.title.x=element_blank()) + 
  ylab("Copy number") +
  geom_hline(yintercept=c(0,1,2,3),linetype="dashed")
print(q)
dev.off()
}