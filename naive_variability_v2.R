setwd('~/Dropbox/work/Matthijs-Olivier/dPCR/CNVReplicates/Rcode/NIPDdataJo/')
source("cleanCode_v0_1Antonov.R")
setwd('~/Dropbox/work/Matthijs-Olivier/dPCR/CNVReplicates/Data/NIPD')
load('~/Dropbox/Work/Matthijs-Olivier/dPCR/CNVReplicates/Data/NIPD/naive_variability.RData')
genes<-c("AMELX","LAMA3","IL13RA1","LNX2","SMCHD1","ATP11C","ITGBL1","DSCR3","AMELY","C18ORF62","CLIC6","SRY","SYNJ","RPP30")[order(c(10,5,11,2,4,12,3,9,13,6,8,14,7,1))]

par(las=2)
par(mar=c(1,3,1,1))
par(mfrow=c(14,1))


IDs<-141:154
plot(0, 0, xlim = c(1, 14), ylim = c(-2.3,-2.1), type = "n",xlab="",ylab="Estimate",main="",xaxt="n",yaxt="n")
axis(side=2,at=c(-2.30, -2.20, -2.10), labels=c("-2.30", "-2.20", "-2.10"))

points(1:14,results$repl1[IDs][order(c(10,5,11,2,4,12,3,9,13,6,8,14,7,1))],pch=16,col=rgb(0,0,0,0.3))
points(1:14,results$repl2[IDs][order(c(10,5,11,2,4,12,3,9,13,6,8,14,7,1))],pch=16,col=rgb(0,0,0,0.3))
points(1:14,results$repl3[IDs][order(c(10,5,11,2,4,12,3,9,13,6,8,14,7,1))],pch=16,col=rgb(0,0,0,0.3))


IDs<-155:168
plot(0, 0, xlim = c(1, 14), ylim = c(-2.35,-2.15), type = "n",xlab="",ylab="Estimate",main="",xaxt="n",yaxt="n")
axis(side=2,at=c(-2.35, -2.25, -2.15))

points(1:14,results$repl1[IDs][order(c(10,5,11,2,4,12,3,9,13,6,8,14,7,1))],pch=16,col=rgb(0,0,0,0.3))
points(1:14,results$repl2[IDs][order(c(10,5,11,2,4,12,3,9,13,6,8,14,7,1))],pch=16,col=rgb(0,0,0,0.3))
points(1:14,results$repl3[IDs][order(c(10,5,11,2,4,12,3,9,13,6,8,14,7,1))],pch=16,col=rgb(0,0,0,0.3))

IDs<-169:182
plot(0, 0, xlim = c(1, 14), ylim = c(-2.15,-1.95), type = "n",xlab="",ylab="Estimate",main="",xaxt="n",yaxt="n")
axis(side=2,at=c(-2.15, -2.05, -1.95))

points(1:14,results$repl1[IDs][order(c(10,5,11,2,4,12,3,9,13,6,8,14,7,1))],pch=16,col=rgb(0,0,0,0.3))
points(1:14,results$repl2[IDs][order(c(10,5,11,2,4,12,3,9,13,6,8,14,7,1))],pch=16,col=rgb(0,0,0,0.3))
points(1:14,results$repl3[IDs][order(c(10,5,11,2,4,12,3,9,13,6,8,14,7,1))],pch=16,col=rgb(0,0,0,0.3))

IDs<-183:196
plot(0, 0, xlim = c(1, 14), ylim = c(-2.05,-1.85), type = "n",xlab="",ylab="Estimate",main="",xaxt="n",yaxt="n")
axis(side=2,at=c(-2.05, -1.95, -1.85))

points(1:14,results$repl1[IDs][order(c(10,5,11,2,4,12,3,9,13,6,8,14,7,1))],pch=16,col=rgb(0,0,0,0.3))
points(1:14,results$repl2[IDs][order(c(10,5,11,2,4,12,3,9,13,6,8,14,7,1))],pch=16,col=rgb(0,0,0,0.3))
points(1:14,results$repl3[IDs][order(c(10,5,11,2,4,12,3,9,13,6,8,14,7,1))],pch=16,col=rgb(0,0,0,0.3))

plot(0, 0, xlim = c(1, 14), ylim = c(-2.6,-2.1), type = "n",xlab="",ylab="Estimate",main="",xaxt="n",yaxt="n")
axis(side=2,at=c(-2.60, -2.35, -2.10))

#axis(1,at=1:14,labels=genes)
points(1:14,results$repl1[1:14][order(c(10,5,11,2,4,12,3,9,13,6,8,14,7,1))],pch=16,col=rgb(0,0,0,0.3))
points(1:14,results$repl2[1:14][order(c(10,5,11,2,4,12,3,9,13,6,8,14,7,1))],pch=16,col=rgb(0,0,0,0.3))
points(1:14,results$repl3[1:14][order(c(10,5,11,2,4,12,3,9,13,6,8,14,7,1))],pch=16,col=rgb(0,0,0,0.3))

IDs<-15:28
plot(0, 0, xlim = c(1, 14), ylim = c(-2.3,-2.0), type = "n",xlab="",ylab="Estimate",main="",xaxt="n",yaxt="n")
axis(side=2,at=c(-2.30, -2.15, -2.00))

#axis(1,at=1:14,labels=genes)
points(1:14,results$repl1[IDs][order(c(10,5,11,2,4,12,3,9,13,6,8,14,7,1))],pch=16,col=rgb(0,0,0,0.3))
points(1:14,results$repl2[IDs][order(c(10,5,11,2,4,12,3,9,13,6,8,14,7,1))],pch=16,col=rgb(0,0,0,0.3))
points(1:14,results$repl3[IDs][order(c(10,5,11,2,4,12,3,9,13,6,8,14,7,1))],pch=16,col=rgb(0,0,0,0.3))


IDs<-29:42
plot(0, 0, xlim = c(1, 14), ylim = c(-2.4,-2.0), type = "n",xlab="",ylab="Estimate",main="",xaxt="n",yaxt="n")
axis(side=2,at=c(-2.40, -2.20, -2.00),labels=c("-2.40", "-2.20", "-2.00"))

#axis(1,at=1:14,labels=genes)
points(1:14,results$repl1[IDs][order(c(10,5,11,2,4,12,3,9,13,6,8,14,7,1))],pch=16,col=rgb(0,0,0,0.3))
points(1:14,results$repl2[IDs][order(c(10,5,11,2,4,12,3,9,13,6,8,14,7,1))],pch=16,col=rgb(0,0,0,0.3))
points(1:14,results$repl3[IDs][order(c(10,5,11,2,4,12,3,9,13,6,8,14,7,1))],pch=16,col=rgb(0,0,0,0.3))

#sample 10 excluded from all analyses
#IDs<-43:56
#plot(0, 0, xlim = c(1, 14), ylim = range(results[IDs,3:5],na.rm=T), type = "n",xlab="",ylab="Estimate",main="",xaxt="n",yaxt="n")
#axis(1,at=1:14,labels=genes)
#points(1:14,results$repl1[IDs][order(c(10,5,11,2,4,12,3,9,13,6,8,14,7,1))],pch=16,col=rgb(0,0,0,0.3))
#points(1:14,results$repl2[IDs][order(c(10,5,11,2,4,12,3,9,13,6,8,14,7,1))],pch=16,col=rgb(0,0,0,0.3))
#points(1:14,results$repl3[IDs][order(c(10,5,11,2,4,12,3,9,13,6,8,14,7,1))],pch=16,col=rgb(0,0,0,0.3))

IDs<-57:70
plot(0, 0, xlim = c(1, 14), ylim = c(-2.3,-2.0), type = "n",xlab="",ylab="Estimate",main="",xaxt="n",yaxt="n")
axis(side=2,at=c(-2.30, -2.15, -2.00),labels=c("-2.30", "-2.15", "-2.00"))

#axis(1,at=1:14,labels=genes)
points(1:14,results$repl1[IDs][order(c(10,5,11,2,4,12,3,9,13,6,8,14,7,1))],pch=16,col=rgb(0,0,0,0.3))
points(1:14,results$repl2[IDs][order(c(10,5,11,2,4,12,3,9,13,6,8,14,7,1))],pch=16,col=rgb(0,0,0,0.3))
points(1:14,results$repl3[IDs][order(c(10,5,11,2,4,12,3,9,13,6,8,14,7,1))],pch=16,col=rgb(0,0,0,0.3))

IDs<-71:84
plot(0, 0, xlim = c(1, 14), ylim = c(-2.4,-2.0), type = "n",xlab="",ylab="Estimate",main="",xaxt="n",yaxt="n")
axis(side=2,at=c(-2.40, -2.20, -2.00),labels=c("-2.40", "-2.20", "-2.00"))

#axis(1,at=1:14,labels=genes)
points(1:14,results$repl1[IDs][order(c(10,5,11,2,4,12,3,9,13,6,8,14,7,1))],pch=16,col=rgb(0,0,0,0.3))
points(1:14,results$repl2[IDs][order(c(10,5,11,2,4,12,3,9,13,6,8,14,7,1))],pch=16,col=rgb(0,0,0,0.3))
points(1:14,results$repl3[IDs][order(c(10,5,11,2,4,12,3,9,13,6,8,14,7,1))],pch=16,col=rgb(0,0,0,0.3))

IDs<-85:98
plot(0, 0, xlim = c(1, 14), ylim = c(-2.5,-2.2), type = "n",xlab="",ylab="Estimate",main="",xaxt="n",yaxt="n")
axis(side=2,at=c(-2.50, -2.35, -2.20))

#axis(1,at=1:14,labels=genes)
points(1:14,results$repl1[IDs][order(c(10,5,11,2,4,12,3,9,13,6,8,14,7,1))],pch=16,col=rgb(0,0,0,0.3))
points(1:14,results$repl2[IDs][order(c(10,5,11,2,4,12,3,9,13,6,8,14,7,1))],pch=16,col=rgb(0,0,0,0.3))
points(1:14,results$repl3[IDs][order(c(10,5,11,2,4,12,3,9,13,6,8,14,7,1))],pch=16,col=rgb(0,0,0,0.3))

IDs<-99:112
plot(0, 0, xlim = c(1, 14), ylim = c(-2.5,-2.3), type = "n",xlab="",ylab="Estimate",main="",xaxt="n",yaxt="n")
axis(side=2,at=c(-2.50, -2.40, -2.30),labels=c("-2.50", "-2.40", "-2.30"))

#axis(1,at=1:14,labels=genes)
points(1:14,results$repl1[IDs][order(c(10,5,11,2,4,12,3,9,13,6,8,14,7,1))],pch=16,col=rgb(0,0,0,0.3))
points(1:14,results$repl2[IDs][order(c(10,5,11,2,4,12,3,9,13,6,8,14,7,1))],pch=16,col=rgb(0,0,0,0.3))
points(1:14,results$repl3[IDs][order(c(10,5,11,2,4,12,3,9,13,6,8,14,7,1))],pch=16,col=rgb(0,0,0,0.3))

IDs<-113:126
plot(0, 0, xlim = c(1, 14), ylim = c(-3.9,-3.7), type = "n",xlab="",ylab="Estimate",main="",xaxt="n",yaxt="n")
axis(side=2,at=c(-3.90, -3.80, -3.70),labels=c("-3.90", "-3.80", "-3.70"))

#axis(1,at=1:14,labels=genes)
points(1:14,results$repl1[IDs][order(c(10,5,11,2,4,12,3,9,13,6,8,14,7,1))],pch=16,col=rgb(0,0,0,0.3))
points(1:14,results$repl2[IDs][order(c(10,5,11,2,4,12,3,9,13,6,8,14,7,1))],pch=16,col=rgb(0,0,0,0.3))
points(1:14,results$repl3[IDs][order(c(10,5,11,2,4,12,3,9,13,6,8,14,7,1))],pch=16,col=rgb(0,0,0,0.3))


IDs<-127:140
plot(0, 0, xlim = c(1, 14), ylim = c(-2.3,-2.1), type = "n",xlab="",ylab="Estimate",main="",xaxt="n",yaxt="n")
axis(side=2,at=c(-2.30, -2.20, -2.10),labels=c("-2.30", "-2.20", "-2.10"))

axis(1,at=1:14,labels=genes)
points(1:14,results$repl1[IDs][order(c(10,5,11,2,4,12,3,9,13,6,8,14,7,1))],pch=16,col=rgb(0,0,0,0.3))
points(1:14,results$repl2[IDs][order(c(10,5,11,2,4,12,3,9,13,6,8,14,7,1))],pch=16,col=rgb(0,0,0,0.3))
points(1:14,results$repl3[IDs][order(c(10,5,11,2,4,12,3,9,13,6,8,14,7,1))],pch=16,col=rgb(0,0,0,0.3))



