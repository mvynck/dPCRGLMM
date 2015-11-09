source("cleanCode_v0_1Antonov.R")

#gene name;	pos;	neg;	ref;	target;	repl;

data.full<-read.csv("NIPDmerged.csv",sep=";")
subjects<-data.full[13:15,1]
data.full<-data.full[,-1]
subject<-5
all.CNV<-list()
#for(subject in 1:length(unique(subjects))){
	cat("Subject",subject,"\n\n")
	data<-data.full[(subject*3-2):(subject*3),]
	
	#targets
	AMELX <- cbind(c(rep("AMELX",3)),data[,1],data[,2],c(rep(TRUE,3)))
	ATP11C <- cbind(c(rep("ATP11C",3)),data[,11],data[,12],c(rep(TRUE,3)))
	IL13RA1 <- cbind(c(rep("IL13RA1",3)),data[,5],data[,6],c(rep(TRUE,3)))
	AMELY <- cbind(c(rep("AMELY",3)),data[,17],data[,18],c(rep(TRUE,3)))
	SRY <- cbind(c(rep("SRY",3)),data[,23],data[,24],c(rep(TRUE,3)))
	C18ORF62 <- cbind(c(rep("C18ORF62",3)),data[,19],data[,20],c(rep(TRUE,3)))
	LAMA3 <- cbind(c(rep("LAMA3",3)),data[,3],data[,4],c(rep(TRUE,3)))
	SMCHD1 <- cbind(c(rep("SMCHD1",3)),data[,9],data[,10],c(rep(TRUE,3)))

	#references
	CLIC6 <- cbind(c(rep("CLIC6",3)),data[,21],data[,22],c(rep(FALSE,3)))
	DSCR3 <- cbind(c(rep("DSCR3",3)),data[,15],data[,16],c(rep(FALSE,3)))
	SYNJ1 <- cbind(c(rep("SYNJ1",3)),data[,25],data[,26],c(rep(FALSE,3)))
	ITGBL1 <- cbind(c(rep("ITGBL1",3)),data[,13],data[,14],c(rep(FALSE,3)))
	LNX2 <- cbind(c(rep("LNX2",3)),data[,7],data[,8],c(rep(FALSE,3)))
	RPP30 <- cbind(c(rep("RPP30",3)),data[,27],data[,28],c(rep(FALSE,3)))

	#dataframe
	datamat<-data.frame(rbind(CLIC6,DSCR3,SYNJ1,AMELX,AMELY,SRY,ATP11C,IL13RA1,C18ORF62,ITGBL1,LAMA3,LNX2,SMCHD1,RPP30))
	colnames(datamat)<-c("gene","pos","neg","target")
	datamat$gene<-as.character(datamat$gene)
	datamat$pos<-as.numeric(as.character(datamat$pos))
	datamat$neg<-as.numeric(as.character(datamat$neg))
	datamat$target<-as.logical(datamat$target)
	datamat<-datamat[!is.na(datamat$pos),]
	
	
	ref.data<-datamat[datamat$target==F,]
	ref.genes<-unique(ref.data$gene)
	stability<-determine.stability(datamat)
	CNV.est<-list()
	for(i in 1:length(ref.genes)){
		use.ref <- order(abs(stability[[2]][,3]),decreasing=T)[1:i]
		use.ref.data <- ref.data[!is.na(match(ref.data$gene,stability[[2]][use.ref,1])),]
		use.data <- rbind(datamat[datamat$target==T,],use.ref.data)
		CNV.est[[i]]<-get.estimates(use.data,fit.set(use.data))
	}
	pdf(paste0("refGeneStability",unique(subjects)[1],".pdf"))
	par(mfrow=c(1,3))
	MAD<-array(0,length(ref.genes))
	for(i in 1:length(ref.genes)){
		MAD[i]<-mean(abs(CNV.est[[i]][,2]-round(CNV.est[[i]][,2])))
	}
	plot(MAD,type="l",ylim=c(0,max(MAD)*1.1),xlab="Number of reference genes",ylab="Mean absolute deviation")
	uncertainty<-array(0,length(ref.genes))
	for(i in 1:length(ref.genes)){
		uncertainty[i]<-mean(CNV.est[[i]][,4]-CNV.est[[i]][,3])
	}
	plot(uncertainty,type="l",ylim=c(0,max(uncertainty)*1.1),xlab="Number of reference genes",ylab="Confidence interval width")
	
	integer<-array(0,length(ref.genes))
	for(i in 1:length(ref.genes)){
		integer[i]<-sum(floor(CNV.est[[i]][,4])-floor(CNV.est[[i]][,3]))
	}
	barplot(integer,xlab="Number of reference genes",ylab="Number of integer copy numbers",names.arg=1:6)
	dev.off()
	all.CNV[[subject]]<-CNV.est
	pdf(paste0("refGeneStability",unique(subjects)[1],"StabPlot.pdf"))
		plot.stability(stability)
	dev.off()
#}
save.image(paste0("refGeneStability",subjects[1],".RData"))