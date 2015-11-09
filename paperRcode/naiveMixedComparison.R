setwd('~/Dropbox/work/Matthijs-Olivier/dPCR/CNVReplicates/Rcode/NIPDdataJo/')
source("cleanCode_v0_1Antonov.R")
setwd('~/Dropbox/work/Matthijs-Olivier/dPCR/CNVReplicates/Data/NIPD')


#gene name;	pos;	neg;	ref;	target;	repl;
results<-matrix(NA,196,6)
results<-data.frame(results)

data.full<-read.csv("NIPDmerged.csv",sep=";")
for(i in 1:length(unique(data.full[,1]))){
	data.temp<-data.full[which(data.full[,1]==unique(data.full[,1])[i]),]
	for(j in 1:14){
		data.temp.sub<-data.temp[,((j*2):(j*2+1))]
		data.temp.sub<-data.frame(gene=paste0("Gene",j),data.temp.sub,target=TRUE)
		data.temp.sub<-data.temp.sub[!is.na(data.temp.sub[,2]),]
		data.temp.sub<-data.temp.sub[!(data.temp.sub[,2])<10,]
		if(nrow(data.temp.sub)>1){
		colnames(data.temp.sub)[2:3]<-c("pos","neg")
		db<-generate.df(data.temp.sub)
		model.rand<-glmer(Y~1+(1|repl),data=db,family=binomial(cloglog),nAGQ=1,verbose=F)
		model.naive<-glm(Y~1,data=db,family=binomial(cloglog))
		results[(j+(i-1)*14),]<-c(i,j,model.rand@beta[1],as.numeric(coef(model.naive)),as.numeric(vcov(model.rand)),vcov(model.naive))
		} else {
			results[(j+(i-1)*14),1:2]<-c(i,j)
		}
	}
}
colnames(results)<-c("sample","gene","mixed.est","naive.est","mixed.var","naive.var")
results2<-results[!is.na(results[,6]),]
hist(results2[,6]/results2[,5],breaks=seq(0,1.05,0.05),xlim=c(0,1.2),ylim=c(0,120),xlab="Variance ratio",main="")
hist(results2[,4]/results2[,3],breaks=20,xlim=c(0.992,1.002),xlab="Estimate ratio",main="")
save(results,file="naive_absquant.RData")