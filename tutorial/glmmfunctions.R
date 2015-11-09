library(lme4)

generate.df <- function(datamat){
	cat("Converting data ...\n")
	#target
	target.df <- datamat[datamat$target==T,]
	if(nrow(target.df)!=0){
	ID<-letters[1]
	Y <- c(rep(1,sum(target.df$pos)),rep(0,sum(target.df$neg)))
	target <- rep("yes",sum(c(target.df$pos,target.df$neg)))
	sample <- rep(ID,sum(c(target.df$pos,target.df$neg)))
	repl<-NULL
	for(i in 1:nrow(target.df)){
		repl<-c(repl,rep(paste(ID,i,sep=""),target.df$pos[i]))
	}
	for(i in 1:nrow(target.df)){
		repl<-c(repl,rep(paste(ID,i,sep=""),target.df$neg[i]))
	}
	db<-data.frame(Y,target,sample,repl)
	}
	
	#reference
	ref.df<-datamat[datamat$target==F,]
	if(nrow(ref.df)!=0){
		for(i in 2:length(unique(datamat$gene))){
			ID<-letters[i]
			ref.df.subset <- ref.df[which(ref.df$gene==unique(ref.df$gene)[i-1]),]
			Y <- c(rep(1,sum(ref.df.subset$pos)),rep(0,sum(ref.df.subset$neg)))
			target <- rep("no",sum(c(ref.df.subset$pos,ref.df.subset$neg)))
			sample <- rep(ID,sum(c(ref.df.subset$pos,ref.df.subset$neg)))
			
			repl<-NULL
			for(j in 1:nrow(ref.df.subset)){
				repl<-c(repl,rep(paste(ID,j,sep=""),ref.df.subset$pos[j]))
			}
			
			for(j in 1:nrow(ref.df.subset)){
				repl<-c(repl,rep(paste(ID,j,sep=""),ref.df.subset$neg[j]))
			}	
			
			db2<-data.frame(Y,target,sample,repl)
			db<-rbind(db,db2)
		}
	}
	return(db)
}



fit.GLMM <- function(datamat){
	#target
	db <- generate.df(datamat)
	model.rand<-NULL
	type<-length(unique(datamat$gene))>1
	repl<-max(table(datamat$gene))>1
	cat("Fitting model: \n")
	n.ref<-length(unique(datamat$gene[datamat$target==F]))
	mult.ref<-n.ref>1
	if(mult.ref&&type){cat("\tmultiple (n=", n.ref ,") reference genes\n",sep="")}else if(type&&!mult.ref){cat("\tsingle reference gene\n")}
	if(type==T&&repl==T){
		cat("\tcopy number variation\n\twith replicates ... \n")
		model.rand<-glmer(Y~target+(1|repl)+(1|sample),data=db,family=binomial(cloglog),nAGQ=1,verbose=F)
	} else if(type==T&&repl==F) {
		cat("\tcopy number variation\n\twithout replicates ... \n")
		model.rand<-glmer(Y~target+(1|sample),data=db,family=binomial(cloglog),nAGQ=1,verbose=F)	
	} else if(type==F&&repl==T){
		cat("\tabsolute quantification\n\twith replicates ... \n")
		model.rand<-glmer(Y~1+(1|repl),data=db,family=binomial(cloglog),nAGQ=1,verbose=F)
	} else if(type==F&&repl==F){
		cat("\tabsolute quantification\n\twithout replicates ... \n")
		model.rand<-glm(Y~1,data=db,family=binomial(cloglog))		
	}
	return(model.rand)
}

fit.set <- function(datamat){
	n.target<-length(unique(datamat$gene[datamat$target==T]))
	if(n.target>1){
		#multiple targets
		cat("Analyzing multiple targets ... \n")
		model.fit<-vector("list",n.target)
		for(i in 1:n.target){
			target.name<-as.character(unique(datamat$gene[datamat$target==T])[i])
			cat("Analyzing target ", i, " of ", n.target," (",target.name,") ...\n",sep="")
			datamat.temp <- rbind(datamat[datamat$gene==target.name,],datamat[datamat$target==F,])
			model.fit[[i]]<-fit.GLMM(datamat.temp)
			cat("\n\n")
		}
	} else if (n.target==1){
		#single target
		cat("Analyzing single target ... \n")
		model.fit<-fit.GLMM(datamat)
	} else {
		cat("No target\n")
	}
	return(model.fit)
}

get.estimates <- function(datamat,model.fit){
	genenames<-as.character(unique(datamat$gene[datamat$target==T]))
	est<-array(0,length(model.fit))
	est.low<-array(0,length(model.fit))
	est.high<-array(0,length(model.fit))

	est.list<-data.frame(genenames,est,est.low,est.high)
	
	#multiple models
	if(is(model.fit)[1]!="glmerMod"&&is(model.fit)[1]!="glm"){
		for(i in 1:length(model.fit)){
			model.rand <- model.fit[[i]]
			if(is(model.fit[[i]])[1]=="glmerMod"){
				if(dim(vcov(model.rand))[1]==1){
					#absolute quantification
					est.list[i,2]<-(exp(-(model.rand@beta[1]+0.5*summary(model.rand)$varcor$repl[1])))*2
					est.list[i,3]<-(exp(-(model.rand@beta[1]+0.5*summary(model.rand)$varcor$repl[1])-1.96*sqrt(vcov(model.rand)[1,1])))*2
					est.list[i,4]<-(exp(-(model.rand@beta[1]+0.5*summary(model.rand)$varcor$repl[1])+1.96*sqrt(vcov(model.rand)[1,1])))*2					
				} else {
					#copy number variation
					est.list[i,2]<-(exp(-(model.rand@beta[2]+summary(model.rand)$varcor$sample[1]/(length(levels(model.rand@frame$sample))-1))))*2
					est.list[i,3]<-(exp(-(model.rand@beta[2]+summary(model.rand)$varcor$sample[1]/(length(levels(model.rand@frame$sample))-1))-1.96*sqrt(vcov(model.rand)[2,2])))*2
					est.list[i,4]<-(exp(-(model.rand@beta[2]+summary(model.rand)$varcor$sample[1]/(length(levels(model.rand@frame$sample))-1))+1.96*sqrt(vcov(model.rand)[2,2])))*2					
				}
			} else {
				est.list[i,2]<-(exp(-coef(model.rand)))*2
				est.list[i,3]<-(exp(-coef(model.rand)-1.96*sqrt(vcov(model.rand))))*2
				est.list[i,4]<-(exp(-coef(model.rand)+1.96*sqrt(vcov(model.rand))))*2				
			}
		}
		
	#single model
	} else {
			if(is(model.fit)[1]=="glmerMod"){
				if(dim(vcov(model.fit))[1]==1){
					#absolute quantification
					est.list[1,2]<-(exp(-(model.fit@beta[1]+0.5*summary(model.fit)$varcor$repl[1])))*2
					est.list[1,3]<-(exp(-(model.fit@beta[1]+0.5*summary(model.fit)$varcor$repl[1])-1.96*sqrt(vcov(model.fit)[1,1])))*2
					est.list[1,4]<-(exp(-(model.fit@beta[1]+0.5*summary(model.fit)$varcor$repl[1])+1.96*sqrt(vcov(model.fit)[1,1])))*2	
				} else {
					#copy number variation
					est.list[1,2]<-(exp(-(model.fit@beta[2]+summary(model.fit)$varcor$sample[1]/(length(levels(model.fit@frame$sample))-1))))*2
					est.list[1,3]<-(exp(-(model.fit@beta[2]+summary(model.fit)$varcor$sample[1]/(length(levels(model.fit@frame$sample))-1))-1.96*sqrt(vcov(model.fit)[2,2])))*2
					est.list[1,4]<-(exp(-(model.fit@beta[2]+summary(model.fit)$varcor$sample[1]/(length(levels(model.fit@frame$sample))-1))+1.96*sqrt(vcov(model.fit)[2,2])))*2					
				}
			} else {
				#GLM only for absolute quantification with no replicates
				est.list[1,2]<-(exp(-coef(model.fit)))*2
				est.list[1,3]<-(exp(-coef(model.fit)-1.96*sqrt(vcov(model.fit))))*2
				est.list[1,4]<-(exp(-coef(model.fit)+1.96*sqrt(vcov(model.fit))))*2				
			}
			est.list<-est.list[1,]
	}
	return(est.list)
}


determine.stability <- function(datamat){

	ref.data<-datamat[datamat$target==F,]
	ref.genes<-unique(ref.data$gene)
	est<-array(NA,(length(ref.genes)-1))
	est.low<-array(NA,(length(ref.genes)-1))
	est.high<-array(NA,(length(ref.genes)-1))
	gene<-array(NA,(length(ref.genes)-1))
	stability <- data.frame(gene,est,est.low,est.high)
	stab.genes<-list()
	stab.summary<-matrix(NA,length(ref.genes),2)
	for(i in 1:length(ref.genes)){
		cat("Determining stability...\n\tReference gene ",i," of ",length(ref.genes),"\n\n")
		target.ref.ID <- which(ref.data$gene==ref.genes[i])
		target.ref.data <- ref.data[target.ref.ID,]
		target.ref.data$target <- TRUE
		ref.ref.data <- ref.data[-target.ref.ID,]
		ref.ref.genes<-unique(ref.ref.data$gene)			
		for(j in 1:length(ref.ref.genes)){
			ref.ref.ID <- which(ref.ref.data$gene==ref.ref.genes[j])
			data.subset <- rbind(target.ref.data,ref.ref.data[ref.ref.ID,])
			stability[j,]<-get.estimates(data.subset,fit.set(data.subset))
		}
		stab.genes[[i]]<-stability
		stab.summary[i,]<-c(mean(stab.genes[[i]]$est),sum((stab.genes[[i]]$est-2)^2))
	}
	stab.summary<-data.frame(as.character(ref.genes),stab.summary)
	colnames(stab.summary)<-c("Gene","Mean","Deviation")
	return(list(stab.genes,stab.summary))
}

plot.stability <- function(stability){
	barplot(sort(abs(stability[[2]][,3])),names.arg=stability[[2]][,1][order(abs(stability[[2]][,3]))])
}