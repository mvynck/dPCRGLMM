load('~/Dropbox/Work/Matthijs-Olivier/dPCR/CNVReplicates/Rcode/NIPDdataJo/antonovScripts/decreasingStability/refGeneStability05_sample.RData')
datamat$target[!datamat$gene=="RPP30"]<-TRUE
datamat$target[datamat$gene=="RPP30"]<-FALSE
a<-fit.set(datamat)
CNV.calc<-get.estimates(datamat,a)
ggplot(CNV.calc, aes(x = CNV.calc[1:13,1], y = CNV.calc[1:13,2])) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymax = CNV.calc[1:13,4], ymin = CNV.calc[1:13,3])) +
  theme_bw()+theme(axis.title.x=element_blank()) + 
  ylab("Copy number") +
  geom_hline(yintercept=c(0,1,2,3),linetype="dashed")

load('~/Dropbox/Work/Matthijs-Olivier/dPCR/CNVReplicates/Rcode/NIPDdataJo/antonovScripts/decreasingStability/refGeneStability06_sample.RData')
datamat$target[!datamat$gene=="RPP30"]<-TRUE
datamat$target[datamat$gene=="RPP30"]<-FALSE
a<-fit.set(datamat)
CNV.calc<-get.estimates(datamat,a)
ggplot(CNV.calc, aes(x = CNV.calc[1:11,1], y = CNV.calc[1:11,2])) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymax = CNV.calc[1:11,4], ymin = CNV.calc[1:11,3])) +
  theme_bw()+theme(axis.title.x=element_blank()) + 
  ylab("Copy number") +
  geom_hline(yintercept=c(0,1,2,3),linetype="dashed")


load('~/Dropbox/Work/Matthijs-Olivier/dPCR/CNVReplicates/Rcode/NIPDdataJo/antonovScripts/decreasingStability/refGeneStability09_sample.RData')
datamat$target[!datamat$gene=="RPP30"]<-TRUE
datamat$target[datamat$gene=="RPP30"]<-FALSE
a<-fit.set(datamat)
CNV.calc<-get.estimates(datamat,a)
ggplot(CNV.calc, aes(x = CNV.calc[1:11,1], y = CNV.calc[1:11,2])) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymax = CNV.calc[1:11,4], ymin = CNV.calc[1:11,3])) +
  theme_bw()+theme(axis.title.x=element_blank()) + 
  ylab("Copy number") +
  geom_hline(yintercept=c(0,1,2,3),linetype="dashed")
  
  
  
  
  
load('~/Dropbox/Work/Matthijs-Olivier/dPCR/CNVReplicates/Rcode/NIPDdataJo/antonovScripts/decreasingStability/refGeneStability13_sample.RData')
datamat$target[!datamat$gene=="RPP30"]<-TRUE
datamat$target[datamat$gene=="RPP30"]<-FALSE
a<-fit.set(datamat)
CNV.calc<-get.estimates(datamat,a)
ggplot(CNV.calc, aes(x = CNV.calc[1:13,1], y = CNV.calc[1:13,2])) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymax = CNV.calc[1:13,4], ymin = CNV.calc[1:13,3])) +
  theme_bw()+theme(axis.title.x=element_blank()) + 
  ylab("Copy number") +
  geom_hline(yintercept=c(0,1,2,3),linetype="dashed")
  
  
  
load('~/Dropbox/Work/Matthijs-Olivier/dPCR/CNVReplicates/Rcode/NIPDdataJo/antonovScripts/decreasingStability/refGeneStability15_sample.RData')
datamat$target[!datamat$gene=="RPP30"]<-TRUE
datamat$target[datamat$gene=="RPP30"]<-FALSE
a<-fit.set(datamat)
CNV.calc<-get.estimates(datamat,a)
ggplot(CNV.calc, aes(x = CNV.calc[1:11,1], y = CNV.calc[1:11,2])) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymax = CNV.calc[1:11,4], ymin = CNV.calc[1:11,3])) +
  theme_bw()+theme(axis.title.x=element_blank()) + 
  ylab("Copy number") +
  geom_hline(yintercept=c(0,1,2,3),linetype="dashed")
  
  
load('~/Dropbox/Work/Matthijs-Olivier/dPCR/CNVReplicates/Rcode/NIPDdataJo/antonovScripts/decreasingStability/refGeneStability16_sample.RData')
datamat$target[!datamat$gene=="RPP30"]<-TRUE
datamat$target[datamat$gene=="RPP30"]<-FALSE
a<-fit.set(datamat)
CNV.calc<-get.estimates(datamat,a)
ggplot(CNV.calc, aes(x = CNV.calc[1:13,1], y = CNV.calc[1:13,2])) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymax = CNV.calc[1:13,4], ymin = CNV.calc[1:13,3])) +
  theme_bw()+theme(axis.title.x=element_blank()) + 
  ylab("Copy number") +
  geom_hline(yintercept=c(0,1,2,3),linetype="dashed")
  
  
load('~/Dropbox/Work/Matthijs-Olivier/dPCR/CNVReplicates/Rcode/NIPDdataJo/antonovScripts/decreasingStability/refGeneStability18_sample.RData')
datamat$target[!datamat$gene=="RPP30"]<-TRUE
datamat$target[datamat$gene=="RPP30"]<-FALSE
a<-fit.set(datamat)
CNV.calc<-get.estimates(datamat,a)
ggplot(CNV.calc, aes(x = CNV.calc[1:13,1], y = CNV.calc[1:13,2])) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymax = CNV.calc[1:13,4], ymin = CNV.calc[1:13,3])) +
  theme_bw()+theme(axis.title.x=element_blank()) + 
  ylab("Copy number") +
  geom_hline(yintercept=c(0,1,2,3),linetype="dashed")
  

load('~/Dropbox/Work/Matthijs-Olivier/dPCR/CNVReplicates/Rcode/NIPDdataJo/antonovScripts/decreasingStability/refGeneStability19_sample.RData')
datamat$target[!datamat$gene=="RPP30"]<-TRUE
datamat$target[datamat$gene=="RPP30"]<-FALSE
a<-fit.set(datamat)
CNV.calc<-get.estimates(datamat,a)
ggplot(CNV.calc, aes(x = CNV.calc[1:11,1], y = CNV.calc[1:11,2])) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymax = CNV.calc[1:11,4], ymin = CNV.calc[1:11,3])) +
  theme_bw()+theme(axis.title.x=element_blank()) + 
  ylab("Copy number") +
  geom_hline(yintercept=c(0,1,2,3),linetype="dashed")
  


load('~/Dropbox/Work/Matthijs-Olivier/dPCR/CNVReplicates/Rcode/NIPDdataJo/antonovScripts/decreasingStability/refGeneStability20_sample.RData')
datamat$target[!datamat$gene=="RPP30"]<-TRUE
datamat$target[datamat$gene=="RPP30"]<-FALSE
a<-fit.set(datamat)
CNV.calc<-get.estimates(datamat,a)
ggplot(CNV.calc, aes(x = CNV.calc[1:13,1], y = CNV.calc[1:13,2])) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymax = CNV.calc[1:13,4], ymin = CNV.calc[1:13,3])) +
  theme_bw()+theme(axis.title.x=element_blank()) + 
  ylab("Copy number") +
  geom_hline(yintercept=c(0,1,2,3),linetype="dashed")
  