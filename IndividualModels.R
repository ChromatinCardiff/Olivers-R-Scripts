library(plotrix)
library(boot)
library(RColorBrewer)
library(plotrix)
library(fields)

put.inNorm <- function(i, j, k) {
	pv <- pnorm(i, mean = j, sd = k)
	#print(paste("Val:", i, "mean:", j, "sd:", k, sep=" "))
	if(pv > 0.5)
	{
		pv <- 1 - pv
	}
	-log10(pv)
}

logPFreq <- function(x) {
	pv <- pnorm(x, mean = FreqSingleMean, sd = FreqSingleSd)

	if(pv > 0.5)
	{
		pv <- 1 - pv
	}
	-log10(pv)
}

BuildDistMatrix <- function(test)
{
	distMat <- matrix(0, nrow(test), nrow(test))
	rSum <- rowSums(test)
	for(i in 1:nrow(test))
	{
		#print(i)
		#region <- test[i]
		comp <- abs(rSum - rSum[i])
		#comp <- apply(test, 1, FUN=function(x) abs(sum(x) - sum(region)))
		distMat[,i] <- comp#abs(rowSums(comp))
	}
	distMat
}

MakeIndividualModels <- function(binSize=10, ...)
{
	inList <- list(...)
	InLen <- length(inList)
	print(paste("Using: ", (InLen-1), " replicates, 1 mutant", sep=""))
		
	for(i in 1:InLen)
	{
		inList[[i]] <- data.frame(inList[[i]])
	}
	
	mix <- matrix(0, nrow(inList[[1]]), ncol(inList[[1]]))
	
	for(i in 1:(InLen-1))
	{
		mix <- mix + inList[[i]]
	}
	
	mix <- mix / (InLen - 1)
	mix[is.na(mix)] <- 0
	
	dT <- matrix(0, nrow(inList[[1]]), ncol(inList[[1]]))
	
	for(i in 1:(InLen-1))
	{
		dT <- dT + abs(mix - as.matrix(inList[[i]]))
	}

	pq <- rowSums(dT)
	
	#for(i in 1:InLen)
	#{
	#	inList[[i]]$Rank <- pq
	#	inList[[i]] <- inList[[i]][order(inList[[i]]$Rank),]
	#	inList[[i]]$Rank <- NULL
	#}

	mix <- matrix(0, nrow(inList[[1]]), ncol(inList[[1]]))
	
	for(i in 1:(InLen-1))
	{
		mix <- mix + inList[[i]]
	}
	mix <- mix / (InLen - 1)
	mix[is.na(mix)] <- 0
	
	means <- matrix(0, nrow(mix), ncol(mix))
	sds <- matrix(0, nrow(mix), ncol(mix))
	
	relDiffs <- list()
	mix <- as.matrix(mix)
	
	for(i in 1:InLen)
	{
		relDiffs[[i]] <- mix - as.matrix(inList[[i]])
		relDiffs[[i]] <- cbind(relDiffs[[i]], pq)
	}

	for(j in 1:ncol(inList[[1]]))
	{
		print(paste("Processing Column ", j, sep=""))
		
		entries <- vector()
		if(j < 5)
		{
			entries <- 1:9
		}
		else if(j > ncol(inList[[1]]) - 4)
		{
			entries <- (ncol(inList[[1]])-8):ncol(inList[[1]])
		}
		else
		{
			entries <- (j-4):(j+4)
		}
		test <- mix[, entries]
		
		allTests <- as.matrix(dist(test, method="manhattan"))
		#allTests <- BuildDistMatrix(test)
		print(paste("Built Dist Matrix ", j, sep=""))
		for(i in 1:nrow(inList[[1]]))
		{
			useDist <- allTests[,i]
			useDist <- cbind(useDist, 1:nrow(test))
			useDist <- useDist[order(useDist[,1]),]
			useDist <- useDist[1:(min(c(floor(nrow(useDist)/10)), 500)),]
			#region <- mix[i, entries]
			
			#Collect the most similar sub-sequences
			#comp <- t(apply(test, 1, FUN=function(x) x - region))
			#comp <- cbind(comp, abs(rowSums(comp)))
			#comp <- cbind(comp, 1:nrow(test))
			#comp <- comp[order(comp[,10]),]
			#comp <- comp[1:floor(nrow(comp)/10),]
			#goGet <- comp[,11]
			goGet <- useDist[,2]
		
			#filter out genes most dissimilar based on replicate variance
			varComp <- relDiffs[[1]][i,ncol(inList[[1]])+1]
			diffs <- abs(c(relDiffs[[1]][goGet, ncol(inList[[1]])+1]) - varComp)
			#comp <- comp[rev(order(diffs)),]
			#comp <- comp[1:floor(nrow(comp)/2),]
			#goGet <- comp[,11]
			
			useDist <- useDist[order(diffs),]
			useDist <- useDist[1:floor(nrow(useDist)/2),]
			goGet <- useDist[,2]
			
			results <- vector()
			
			for(k in 1:(InLen-1))
			{
				results <- c(results, relDiffs[[k]][goGet, j])
			}
			means[i,j] <- mean(results)
			sds[i,j] <- sd(results)
		}
		print(paste("Processed rows, for column ", j, sep=""))
	}
	results <- list()
	for(k in 1:InLen)
	{
		print("Testing all samples under the model...")
		relDiffs[[k]] <- relDiffs[[k]][,-(ncol(inList[[1]]) + 1)]
	
		results[[k]] <- matrix(0, nrow(relDiffs[[k]]), ncol(relDiffs[[k]]))

		for(i in 1:nrow(inList[[1]]))
		{
			results[[k]][i,] <- mapply(put.inNorm, relDiffs[[k]][i,], means[i,], sds[i,])
		}

		results[[k]][is.infinite(results[[k]])] <- 35
	}
	
	print("Scaled Multiple Test Correction...")
	
	maxImg <- matrix(0, nrow(results[[1]]), ncol(results[[1]]))

	for(i in 1:nrow(results[[1]]))
	{
		vecImg <- vector()
		#print(i)
		for(j in 1:(InLen - 1))
		{
			if(j == 1)
			{
				vecImg <- results[[j]][i,]
			}
			else
			{
				vecImg <- mapply(max, vecImg, results[[j]][i,])
			}
		}
		
		maxImg[i,] <- vecImg
	}

	subImg <- results[[InLen]] - maxImg
	subImg[subImg < 0] <- 0.5
	
	retVal <- list()
	retVal[[1]] <- means
	retVal[[2]] <- sds
	retVal[[3]] <- results
	retVal[[4]] <- binSize
	retVal[[5]] <- subImg
	retVal
}

meanBoot <- function(data, indices) {
		d <- data[indices]
		mean(d)
	}

PlotAllRankGraphs <- function(indexes, division, FilePrefix, results, ...)
{
	inL <- list(...)
	RankGenesCore(indexes, division, save=TRUE, justone=FALSE, FilePrefix, results, inL)
}

PlotTopRankMutant <- function(indexes, division, results, ...)
{
	inL <- list(...)
	RankGenesCore(indexes, division, save=FALSE, justone=TRUE, prefix="", results, inL)
}

RankGenesCore <- function(indexes, division, save, justone, prefix, results, inList)
{
	InLen <- length(inList)
	len <- length(inList[[1]][1,])
	muPSums <- rowSums(results[[5]][,indexes[1]:indexes[2]])

	for(i in 1:InLen)
	{
		inList[[i]] <- as.matrix(inList[[i]])
		inList[[i]] <- inList[[i]][order(muPSums),]
	}	
	
	xNam <- seq(-floor(len/2)*results[[4]], floor(len/2)*results[[4]], 50)
	xNam <- xNam[1:(ceiling((len*results[[4]])/50))]
	
	GraphCount <- division
	interGraphs <- 1:GraphCount
	if(justone)
	{
		interGraphs <- GraphCount
	}
	
	colors <- c("cadetblue1", "aquamarine1", "darkolivegreen2", "chocolate1")

	allwt <- list()
	for(i in 1:(InLen-1))
	{
		allwt[[i]] <- list()
	}
	
	for(i in interGraphs)
	{
		members <- floor(length(muPSums)/GraphCount)
		muMeans <- colMeans(inList[[InLen]][(((i-1) * members) + 1) : (i * members),])
	
		for(j in 1:len)
		{
			for(k in 1:(InLen-1))
			{
				allwt[[k]][[j]] <- boot(inList[[k]][(((i-1) * members) + 1) : (i * members),j], meanBoot, 200)
			}
		}
	
		allresults <- matrix(0, 200 * (InLen-1), len)

		for(p in 1:len)
		{
			nextVec <- vector()
			for(k in 1:(InLen-1))
			{
				nextVec <- c(nextVec, allwt[[k]][[p]]$t)
			}
			allresults[,p] <- nextVec
		}

		AllMeans <- apply(allresults, 2, function(x) mean(x))
		AllSds <- apply(allresults, 2, function(x) sd(x))
		AllUpper <- AllMeans + (AllSds * 2)
		AllLower <- AllMeans - (AllSds * 2)
	
		muPvs <- sapply(1:length(muMeans), function(i) put.inNorm(muMeans[i], AllMeans[i], AllSds[i]))

		muPvs[is.infinite(muPvs)] <- 35
		topRes <- max(AllUpper, muMeans)
		lowres <- min(AllLower, muMeans)
		
		if(save)
		{
			pdf(paste(prefix, "_", i, ".pdf", sep=""), width=15, height=10)
		}
		
		par(mfrow=c(2,1), mar=c(0,4,2,2))
		par(fig=c(0,1,0.2,1))
		locations <- c(seq(1,len, len/48), len)
		plot(1:len, AllMeans, ylim=c(lowres,topRes), type="l", lty=rep(1,4), lwd=rep(2,4), ylab="Read Mean Frequencey", main=paste("Chromatin Structure for ", ((i-1)*(100/GraphCount)), "-", (i*(100/GraphCount)), "% of Ranked Mutant Regions", sep=""), xaxt="n")
		polygon(c(1:len, len:1), c(AllUpper, rev(AllLower)), density=NA, col="darkolivegreen2")
		points(1:len, AllMeans, type="l", lwd=1)
		points(1:len, muMeans, type="l", lwd=1.5, col="darkslateblue")
		abline(v=indexes[1])
		abline(v=indexes[2])
		
		par(fig=c(0,1, 0, 0.2), new=TRUE, mar=c(3,4,0,2))
		plot(1:len, muPvs, col="midnightblue", type="h", ylim=c(0,30), lwd="3", ylab="Signif. -log(p)", xaxt = "n", xlab="Base Pairs Distance around Feature")
		abline(h=-log10(5e-02))
		abline(h=-log10((5e-02)/len))
		abline(v=indexes[1])
		abline(v=indexes[2])
		staxlab(1, at=locations, labels=xNam, srt=45,  cex.axis=0.1, line.spacing=0.8, nlines=2, cex=0.7)
		
		if(save)
		{
			dev.off()
		}
	}
}

SaveIndividualModel <- function(results, prefix="model")
{
	InLen <- length(results[[3]])
	len <- length(results[[1]][,1])
	vec1 <- rep(InLen, len)
	vec2 <- rep(results[[4]], len)
	
	finalMatrix <- cbind(vec1, vec2, results[[1]], results[[2]])
	
	for(i in 1:InLen)
	{
		finalMatrix <- cbind(finalMatrix, results[[3]][[i]])
	}
	
	finalMatrix <- cbind(finalMatrix, results[[5]])
	
	colnames(finalMatrix) <- 1:ncol(finalMatrix)
	
	write.table(finalMatrix, paste(prefix, ".indiv", sep=""), sep="\t", quote=FALSE, row.names=FALSE)
}

LoadIndividualModel <- function(fileName)
{
	loaded <- read.table(fileName,header=TRUE)
	
	InLen <- loaded[1,1]
	binSize <- loaded[1,2]
	print(paste("Loaded ", (InLen - 1), ", replicate tests, 1 mutant test. Binsize: ", binSize, sep=""))
	size <- (ncol(loaded) - 2) / (InLen + 2)
	
	retVal <- list()
	retVal[[1]] <- loaded[,3:(size+2)]
	retVal[[2]] <- loaded[,(size+3):((size*2) + 2)]
	
	results <- list()
	
	for(i in 1:InLen)
	{
		results[[i]] <- loaded[,(size * (2+(i-1)) + 3):(size * (2+i) + 2)]
	}
	
	retVal[[3]] <- results
	retVal[[4]] <- binSize
	retVal[[5]] <- loaded[,(size * (2+InLen) + 3):(size * (3+InLen) + 2)]
	retVal
}

PlotIndividualModelTest <- function(results, save=FALSE, prefix="")
{
	rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))(100)
	InLen <- length(results[[3]])
	len <- length(results[[3]][[1]][1,])
	reps <- nrow(results[[3]][[1]])
	muReps <- nrow(results[[5]])

	colR <- ncol(results[[3]])
	
	GImg <- list()
	
	for(i in 1:InLen)
	{
		GImg[[i]] <- t(apply(results[[3]][[i]], 2, function(x) sort(x, decreasing=TRUE)))
	}

	xNam <- seq(-floor(len/2)*results[[4]], floor(len/2)*results[[4]], 50)
	xNam <- xNam[1:(ceiling(len*results[[4]]/50))]
	
	if(save)
	{
		pdf(file=paste(prefix, "_all.pdf",sep=""), height=10, width=15)
	}
	par(mfrow=c(InLen,1), mar=c(3,4,2,2))

	locations <- seq(0,1, 1/((ceiling(len*results[[4]]/50))-1))
	
	for(i in 1:(InLen - 1))
	{
		image.plot(GImg[[i]], col=rf, zlim=c(-log10(5e-02), 5), ylim=c(0,0.55), main =paste("Rep ", i, sep=""), ylab="Proportion of Total Region Set", xaxt = "n")
		abline(h=0.05)
		staxlab(1, at=locations, labels=xNam, srt=45,  cex.axis=0.1, line.spacing=0.8, nlines=2, cex=0.7)
	
	}

	image.plot(GImg[[InLen]], col=rf, zlim=c(-log10(5e-02), 5), ylim=c(0,0.55), main = "Mutant", ylab="Proportion of Total Region Set", xaxt = "n")
	abline(h=0.05)
	staxlab(1, at=locations, labels=xNam, srt=45,  cex.axis=0.1, line.spacing=0.8, nlines=2, cex=0.7)
	
	if(save)
	{
		dev.off()
	}
}

PlotMultiTestCorrection <- function(results, save=FALSE, prefix="mtcorrect")
{
	GImg <- list()
	len <- length(results[[3]][[1]][1,])
	InLen <- length(results[[3]])
	locations <- seq(0,1, 1/(ceiling(len*results[[4]]/50)-1))
	xNam <- seq(-floor(len/2)*results[[4]], floor(len/2)*results[[4]], 50)
	xNam <- xNam[1:(ceiling(len*results[[4]]/50))]
	
	for(i in 1:InLen)
	{
		GImg[[i]] <- results[[3]][[i]]
	}
	
	maxImg <- matrix(0, nrow(GImg[[1]]), ncol(GImg[[1]]))

	for(i in 1:nrow(GImg[[1]]))
	{
		vecImg <- vector()
		#print(i)
		for(j in 1:(InLen - 1))
		{
			if(j == 1)
			{
				vecImg <- GImg[[j]][i,]
			}
			else
			{
				vecImg <- mapply(max, vecImg, GImg[[j]][i,])
			}
		}
		
		maxImg[i,] <- vecImg
	}

	subImg <- GImg[[InLen]] - maxImg
	
	subImg <- apply(subImg, 2, function(x) sort(x, decreasing=TRUE))
	maxImg <- apply(maxImg, 2, function(x) sort(x, decreasing=TRUE))
	GImg[[InLen]] <- apply(GImg[[InLen]], 2, function(x) sort(x, decreasing=TRUE))

	if(save)
	{
		pdf(file=paste(prefix, ".pdf",sep=""), height=10, width=15)
	}

	par(mfrow=c(3,1), mar=c(3,4,2,2))

	image.plot(t(maxImg),zlim=c(-log10(5e-02), 5), ylim=c(0,0.55), main = "Replicate Maximum", ylab="Proportion of Total Region Set", xaxt = "n")
	abline(h=0.05)
	staxlab(1, at=locations, labels=xNam, srt=45,  cex.axis=0.1, line.spacing=0.8, nlines=2, cex=0.7)

	image.plot(t(GImg[[InLen]]), zlim=c(-log10(5e-02), 5), ylim=c(0,0.55), main = "Mutant", ylab="Proportion of Total Region Set", xaxt = "n")
	abline(h=0.05)
	staxlab(1, at=locations, labels=xNam, srt=45,  cex.axis=0.1, line.spacing=0.8, nlines=2, cex=0.7)

	image.plot(t(subImg), zlim=c(-log10(5e-02), 16), ylim=c(0,0.55), main = "Rep.Max corrected Mutant significance levels", ylab="Proportion of Total Region Set", xaxt = "n")
	staxlab(1, at=locations, labels=xNam, srt=45,  cex.axis=0.1, line.spacing=0.8, nlines=2, cex=0.7)

	if(save)
	{
		dev.off()
	}
}

PlotAllRegionGraphs <- function(prefix, results, ...)
{
	inList <- list(...)
	len <- length(inList[[1]][1,])
	InLen <- length(inList)
	rnames <- row.names(inList[[1]])
	XForHeat <- seq(-floor(len/2) * results[[4]], floor(len/2) * results[[4]], results[[4]])
	XForHeat <- XForHeat[1:len]
	count <- 1
	for(genePlot in 1:nrow(inList[[1]]))
	{
		print(paste("Printed Graph: ", rnames[genePlot], ", Number: ", count, sep=""))
		count <- count  + 1
		pdf(paste(prefix, rnames[genePlot], ".pdf", sep=""), width=15,height=10)

		par(mfrow=c(2,1), mar=c(4,4,4,4))

		plot(XForHeat, inList[[1]][genePlot,], type="l", col="cadetblue1", ylab="Read Counts", xlab="Base Pair Distance from Feature", main=paste("Chromatin Structure for ", rnames[genePlot], sep=""))
		
		if(InLen > 2)
		{
			for(i in 2:(InLen -1))
			{
				points(XForHeat, inList[[i]][genePlot,], type="l", col="aquamarine1")
			}
		}
		points(XForHeat, inList[[InLen]][genePlot,], type="l", col="chocolate1")

		plot(XForHeat, results[[5]][genePlot,], type="h", col="darkmagenta", ylab="Negative Log P of Mutant divergence", xlab="Base Pair Distance from Feature", ylim=c(0,8))
		abline(h=-log10(5e-02))

		dev.off()

	}
}

PlotSumSquares <- function(results, MaxNumber=15)
{
	SumSq <- vector()
	SumSq[1] <- (nrow(results[[5]])-1) * sum(apply(results[[5]], 2, var))
	for (i in 2:MaxNumber) 
	{
		SumSq[i] <- sum(kmeans(results[[5]], centers=i)$withinss)
	}
	plot(1:MaxNumber, SumSq, type="b", xlab="Clusters Count", ylab="Sum of Squares, within clusters.", main="Kmeans: Sum of squares reduction by cluster count. Look for the knee!")
}

ClusterRegions <- function(results, clusters)
{
	kmeans(results[[5]], clusters)
}

PlotOneCluster <- function(kclusters, ClusterNumber, ...)
{
	inList <- list(...)
	PlotClustersCore(kclusters, "" ,FALSE, TRUE, ClusterNumber, inList)
} 

PlotAllClusters <- function(kclusters, prefix="", ...)
{
	inList <- list(...)
	PlotClustersCore(kclusters, prefix, TRUE, FALSE, 1, inList)
}

#TODO this here function!
PlotClustersCore <- function(kclusters, prefix="", save=FALSE, JustOne=TRUE, clusterNo=1, inList) 
{
	tot <- nrow(inList[[1]])
	len <- length(inList[[1]][1,])
	numClus <- length(levels(as.factor(kclusters$cluster)))
	xNam <- seq(-floor(len/2)*results[[4]], floor(len/2)*results[[4]], 50)
	xNam <- xNam[1:(ceiling((len*results[[4]])/50))]
	locations <- c(seq(1,len, len/48), len)
	toDraw <- 1:numClus
	InLen <- length(inList)
	
	if(JustOne)
	{
		toDraw <- clusterNo
	}
	
	allwt <- list()
	for(i in 1:(InLen-1))
	{
		allwt[[i]] <- list()
	}

	for(i in toDraw)
	{
		print(paste("Running Bootstrap Model and Drawing Cluster: ", i, sep=""))
		reMuCounts <- inList[[InLen]][kclusters$cluster == i,]
		muMeans <- colMeans(reMuCounts)
		HowMany <- signif((nrow(reMuCounts)/tot) * 100, 4)
		
		for(j in 1:len)
		{
			for(k in 1:(InLen-1))
			{
				allwt[[k]][[j]] <- boot(inList[[k]][kclusters$cluster == i,j], meanBoot, 200)
			}
		}
	
		allresults <- matrix(0, 200 * (InLen - 1), len)

		for(p in 1:len)
		{
			nextVec <- vector()
			for(k in 1:(InLen-1))
			{
				nextVec <- c(nextVec, allwt[[k]][[p]]$t)
			}
			allresults[,p] <- nextVec
		}

		AllMeans <- apply(allresults, 2, function(x) mean(x))
		AllSds <- apply(allresults, 2, function(x) sd(x))
		AllUpper <- AllMeans + (AllSds * 2)
		AllLower <- AllMeans - (AllSds * 2)
	
		muPvs <- sapply(1:length(muMeans), function(i) put.inNorm(muMeans[i], AllMeans[i], AllSds[i]))

		muPvs[is.infinite(muPvs)] <- 35
		
		if(save)
		{
			pdf(paste(prefix, "_", i, ".pdf", sep=""), width=15, height=10)
		}
		
		topres <- max(AllUpper, muMeans)
		lowres <- min(AllLower, muMeans)
		
		par(mfrow=c(2,1), mar=c(0,4,2,2))
		par(fig=c(0,1,0.2,1))
		plot(1:len, AllMeans, ylim=c(lowres,topres), type="l", ylab="Read Mean Frequencey", main=paste("Chromatin Structure for cluster ", i , "/", numClus ,". ", HowMany, "% of set", sep=""), xaxt="n")
		polygon(c(1:len, len:1), c(AllUpper, rev(AllLower)), density=NA, col="bisque2")
		points(1:len, AllMeans, type="l", lwd=1)
		points(1:len, muMeans, type="l", lwd=2, col="blueviolet")
	
		par(fig=c(0,1, 0, 0.2), new=TRUE, mar=c(3,4,0,2))
		plot(1:len, muPvs, col="forestgreen", type="h", ylim=c(0,30), lwd="3", ylab="Signif. -log(p)", xaxt = "n", xlab="Base Pairs Distance around Feature")
		abline(h=-log10(5e-02))
		abline(h=-log10((5e-02)/len))
		staxlab(1, at=locations, labels=xNam, srt=45,  cex.axis=0.1, line.spacing=0.8, nlines=2, cex=0.7)
	
		if(save)
		{
			dev.off()
		}
	}
}
