data.tutorial <- read.csv("http://users.ugent.be/~mvynck/glmm/tutorial.csv")
head(data.tutorial)

#install.packages("lme4")
source("http://users.ugent.be/~mvynck/glmm/glmmfunctions.R")

colnames(data.tutorial)<-c("gene","pos","neg","target")
data.tutorial$gene<-as.character(data.tutorial$gene)
data.tutorial$pos<-as.numeric(as.character(data.tutorial$pos))
data.tutorial$neg<-as.numeric(as.character(data.tutorial$neg))
data.tutorial$target<-as.logical(data.tutorial$target)

fit.example1 <- fit.set(data.tutorial[1,])
(example1 <- get.estimates(data.tutorial[1,],fit.example1))

fit.example2 <- fit.set(data.tutorial[1:3,])
(example2 <- get.estimates(data.tutorial[1:3,],fit.example2))

fit.example3 <- fit.set(data.tutorial[1:4,])
(example3 <- get.estimates(data.tutorial[1:4,],fit.example3))

fit.example4 <- fit.set(data.tutorial[c(1:4,23:25),])
(example4 <- get.estimates(data.tutorial[1:4,],fit.example4))