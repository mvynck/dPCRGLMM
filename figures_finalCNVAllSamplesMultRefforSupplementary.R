load('~/Dropbox/Work/Matthijs-Olivier/dPCR/CNVReplicates/Rcode/NIPDdataJo/antonovScripts/decreasingStability/refGeneStability05_sample.RData')
datamat<-datamat[c(1:24,27:32),]
a<-fit.set(datamat)
CNV.calc<-get.estimates(datamat,a)
ggplot(CNV.calc, aes(x = CNV.calc[1:8,1], y = CNV.calc[1:8,2])) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymax = CNV.calc[1:8,4], ymin = CNV.calc[1:8,3])) +
  theme_bw()+theme(axis.title.x=element_blank()) + 
  ylab("Copy number") +
  geom_hline(yintercept=c(0,1,2,3),linetype="dashed")

load('~/Dropbox/Work/Matthijs-Olivier/dPCR/CNVReplicates/Rcode/NIPDdataJo/antonovScripts/decreasingStability/refGeneStability06_sample.RData')
datamat<- datamat[c(1:11,21,22,26:34),]
a<-fit.set(datamat)
CNV.calc<-get.estimates(datamat,a)
ggplot(CNV.calc, aes(x = CNV.calc[1:3,1], y = CNV.calc[1:3,2])) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymax = CNV.calc[1:3,4], ymin = CNV.calc[1:3,3])) +
  theme_bw()+theme(axis.title.x=element_blank()) + 
  ylab("Copy number") +
  geom_hline(yintercept=c(0,1,2,3),linetype="dashed")


load('~/Dropbox/Work/Matthijs-Olivier/dPCR/CNVReplicates/Rcode/NIPDdataJo/antonovScripts/decreasingStability/refGeneStability09_sample.RData')
datamat<-datamat[c(3:15,24:28),]
a<-fit.set(datamat)
CNV.calc<-get.estimates(datamat,a)
ggplot(CNV.calc, aes(x = CNV.calc[1:3,1], y = CNV.calc[1:3,2])) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymax = CNV.calc[1:3,4], ymin = CNV.calc[1:3,3])) +
  theme_bw()+theme(axis.title.x=element_blank()) + 
  ylab("Copy number") +
  geom_hline(yintercept=c(0,1,2,3),linetype="dashed")

load('~/Dropbox/Work/Matthijs-Olivier/dPCR/CNVReplicates/Rcode/NIPDdataJo/antonovScripts/decreasingStability/refGeneStability13_sample.RData')
datamat<-datamat[1:39,]
a<-fit.set(datamat)
CNV.calc<-get.estimates(datamat,a)
ggplot(CNV.calc, aes(x = CNV.calc[1:8,1], y = CNV.calc[1:8,2])) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymax = CNV.calc[1:8,4], ymin = CNV.calc[1:8,3])) +
  theme_bw()+theme(axis.title.x=element_blank()) + 
  ylab("Copy number") +
  geom_hline(yintercept=c(0,1,2,3),linetype="dashed")
  
  
load('~/Dropbox/Work/Matthijs-Olivier/dPCR/CNVReplicates/Rcode/NIPDdataJo/antonovScripts/decreasingStability/refGeneStability15_sample.RData')
datamat<-datamat[c(1:9,16:33),]
a<-fit.set(datamat)
CNV.calc<-get.estimates(datamat,a)
ggplot(CNV.calc, aes(x = CNV.calc[1:3,1], y = CNV.calc[1:3,2])) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymax = CNV.calc[1:3,4], ymin = CNV.calc[1:3,3])) +
  theme_bw()+theme(axis.title.x=element_blank()) + 
  ylab("Copy number") +
  geom_hline(yintercept=c(0,1,2,3),linetype="dashed")
  
load('~/Dropbox/Work/Matthijs-Olivier/dPCR/CNVReplicates/Rcode/NIPDdataJo/antonovScripts/decreasingStability/refGeneStability16_sample.RData')
datamat<-datamat[c(1:27,31:33,40:42),]
a<-fit.set(datamat)
CNV.calc<-get.estimates(datamat,a)
ggplot(CNV.calc, aes(x = CNV.calc[1:8,1], y = CNV.calc[1:8,2])) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymax = CNV.calc[1:8,4], ymin = CNV.calc[1:8,3])) +
  theme_bw()+theme(axis.title.x=element_blank()) + 
  ylab("Copy number") +
  geom_hline(yintercept=c(0,1,2,3),linetype="dashed")
  
load('~/Dropbox/Work/Matthijs-Olivier/dPCR/CNVReplicates/Rcode/NIPDdataJo/antonovScripts/decreasingStability/refGeneStability18_sample.RData')
a<-fit.set(datamat)
CNV.calc<-get.estimates(datamat,a)
ggplot(CNV.calc, aes(x = CNV.calc[1:8,1], y = CNV.calc[1:8,2])) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymax = CNV.calc[1:8,4], ymin = CNV.calc[1:8,3])) +
  theme_bw()+theme(axis.title.x=element_blank()) + 
  ylab("Copy number") +
  geom_hline(yintercept=c(0,1,2,3),linetype="dashed")

load('~/Dropbox/Work/Matthijs-Olivier/dPCR/CNVReplicates/Rcode/NIPDdataJo/antonovScripts/decreasingStability/refGeneStability19_sample.RData')
a<-fit.set(datamat)
CNV.calc<-get.estimates(datamat,a)
ggplot(CNV.calc, aes(x = CNV.calc[1:3,1], y = CNV.calc[1:3,2])) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymax = CNV.calc[1:3,4], ymin = CNV.calc[1:3,3])) +
  theme_bw()+theme(axis.title.x=element_blank()) + 
  ylab("Copy number") +
  geom_hline(yintercept=c(0,1,2,3),linetype="dashed")


load('~/Dropbox/Work/Matthijs-Olivier/dPCR/CNVReplicates/Rcode/NIPDdataJo/antonovScripts/decreasingStability/refGeneStability20_sample.RData')
a<-fit.set(datamat)
CNV.calc<-get.estimates(datamat,a)
ggplot(CNV.calc, aes(x = CNV.calc[1:8,1], y = CNV.calc[1:8,2])) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymax = CNV.calc[1:8,4], ymin = CNV.calc[1:8,3])) +
  theme_bw()+theme(axis.title.x=element_blank()) + 
  ylab("Copy number") +
  geom_hline(yintercept=c(0,1,2,3),linetype="dashed")