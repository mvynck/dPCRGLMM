############################################################################################
# Copyright 2017 Matthijs Vynck
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Topic: GLMMs for relative expression reference gene selection
# Author: Matthijs.Vynck@UGent.be
# Modified by: -
# Date created: March 2, 2017
# Last updated: March 2, 2017
# Please consider citing Vynck et al. (2016) Biomol Detect Quantif, 9, 1-13. when using.
# More details can be found in section 2.3 (methods) and 3.2 (rationale) of that paper.
############################################################################################

###########################################
# you can run the code until line 65
# line 65 and further is a small tutorial
###########################################

# we will use the lme4 library for GLMM estimation
library(lme4)

# this function calculates the stability
calc.stab.relexp <- function(glmm.df){
	# calculate the variability for all pairwise comparisons
	n <- length(unique(glmm.df$gene))
	var.est <- matrix(NA,n,n)
	for(i in 1:n){
		for(j in 1:n){
			if(i!=j){
			glmm.df.subset <- glmm.df[glmm.df$gene%in%unique(glmm.df$gene)[c(i,j)],]
			fit <- glmer(matrix(c(pos,neg),ncol=2)~1+				#intercept
					dummy(gene,unique(gene)[1])+					#relative expression
					(1|well)+										#replicate effect
					(1|sample),										#sample effect
					data=glmm.df.subset,family=binomial(cloglog), nAGQ=1, verbose=F)
			var.est[i,j]<-vcov(fit)[2,2]
			} else {}
		}
	}

	# calculate reference gene stability
	ref.stab <- rowSums(var.est,na.rm=TRUE)/(n-1)
	ref.df <- data.frame(gene=unique(glmm.df$gene),stability=ref.stab)
	
	# final stability: sorted from most stable to least stable
	# note that a stability value close to zero means that the gene is more stable!
	(ref.final <- ref.df[order(ref.df$stability),])
	return(ref.final)
}




###########################################
# EXAMPLE
###########################################

# construct a matrix containing all data for the candidate reference gene pool
# pos.droplets: the number of positive droplets
# neg.droplets: the number of negative droplets
# gene: a gene identifier, unique for each gene
# sample: sample identiifier, unique for each sample
# well: the well in which the sample was amplified, should e.g. have the same identifier for two genes amplified in duplex

pos.droplets <- c(1000,1100,1500,1600,1700,1760,2020,2080,1200,1300,1600,1700,1800,1860,2160,2260)
neg.droplets <- c(14000,14000,14000,14000,14000,14000,14000,14000,14000,14000,14000,14000,14000,14000,14000,14000)

# suppose we have 4 genes, 2 samples, all performed in duplicate, singleplex
gene <- c(1,1,2,2,3,3,4,4,1,1,2,2,3,3,4,4)
sample <- c(1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2)
well <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)

# put all data into a dataframe
glmm.df <- data.frame(pos=pos.droplets,neg=neg.droplets,gene=factor(gene),sample=factor(sample),well=factor(well))

# you should now have something that looks like this:
#    pos   neg gene sample well
#1  1000 14000    1      1    1
#2  1100 14000    1      1    2
#3  1500 14000    2      1    3
#4  1600 14000    2      1    4
#5  1700 14000    3      1    5
#6  1760 14000    3      1    6
#7  2020 14000    4      1    7
#8  2080 14000    4      1    8
#9  1200 14000    1      2    9
#10 1300 14000    1      2   10
#11 1600 14000    2      2   11
#12 1700 14000    2      2   12
#13 1800 14000    3      2   13
#14 1860 14000    3      2   14
#15 2100 14000    4      2   15
#16 2160 14000    4      2   16

# ready to calculate stability!
(stab <- calc.stab.relexp(glmm.df))

# ordered from most to least stable: 4 - 3 - 2 - 1

# what happens if e.g. gene 3 in sample 2 has a highly variable positive copy number count?
# let's change the third last number (that is gene 3, sample 2, replicate 2) to 4000 positive copies
# we would expect that gene 3 is now classified as being least stable...
pos.droplets <- c(1000,1100,1500,1600,1700,1760,2020,2080,1200,1300,1600,1700,1800,4000,2160,2260)

# adjust dataframe
glmm.df <- data.frame(pos=pos.droplets,neg=neg.droplets,gene=factor(gene),sample=factor(sample),well=factor(well))

# calculate stability
(stab2 <- calc.stab.relexp(glmm.df))

# note that the order of the other genes has changed as well, this is not unexpected:
# the large variability of gene 3 has a large impact on all stability estimates (note that they are all much larger than before)
# the pairwise effects may have been altered dramatically
# it is best to remove gene 3 and recalculate stability

glmm.df <- glmm.df[glmm.df$gene!=3,]
(stab3 <- calc.stab.relexp(glmm.df))

# the original order (4 - 2 - 1) is retrieved
# this may not always be the case because the effect of gene 3 has been removed altogether