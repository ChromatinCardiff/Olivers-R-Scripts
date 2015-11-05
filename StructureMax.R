library(Hmisc)
library(nnet)

StructureMax <- function(MeanSample, PeakL, PeakR, prealign=TRUE, binSize=10, stepSize=5, splineSize=10000)
{
	l1Len <- length(MeanSample[1,])
	l1Entries <- length(MeanSample[,1])
	
	officialTestRange <- 100
	basesFromPeakToSkip <- 10
	x <- 1:l1Len
	scopeWidth <- officialTestRange/binSize
	lb2 <- PeakR - scopeWidth
	hb2 <- PeakR + scopeWidth

	lb <- PeakL - scopeWidth
	hb <- PeakL + scopeWidth
	
	baseToSpline <- splineSize / l1Len
	nucSize <- 160/binSize
	geneVecAdd <- floor(scopeWidth * (baseToSpline))
	rangeSize <- floor(((baseToSpline) / binSize) * (officialTestRange/2))
	sideAlignIndexes <- 1:floor((basesFromPeakToSkip * (baseToSpline/binSize)) / stepSize) 
	iniAlignIndexes <- 1:floor((basesFromPeakToSkip * (baseToSpline/binSize)) / 2)
	
	low <- floor((lb + (scopeWidth / 2)) * (baseToSpline))
	top <- floor((hb - (scopeWidth / 2)) * (baseToSpline))

	low2 <- floor((lb2 + (scopeWidth / 2)) * (baseToSpline))
	top2 <- floor((hb2 - (scopeWidth / 2)) * (baseToSpline))

	XForHeat <- seq(-(floor(l1Len/2)*binSize), floor(l1Len/2)*binSize, binSize)
	xNam <- seq(-(floor(l1Len/2)*binSize), floor(l1Len/2)*binSize, 50)

	results <- sapply(1:nrow(MeanSample), function(p) spline(x, MeanSample[p,], n=splineSize, ties="ordered"))
	results2 <- sapply(1:nrow(MeanSample), function(p) spline(x, MeanSample[p,], n=splineSize, ties="ordered"))

	Peaks <- colMeans(MeanSample)
	TotalMean <- mean(Peaks)
	cnt = 0
	PeakIndex <- vector()

	PeakIndA <- seq(PeakL %% nucSize, PeakL, nucSize)
	PeakIndA <- PeakIndA[PeakIndA > scopeWidth]
	PeakIndB <- seq(PeakR, l1Len, nucSize)
	PeakIndB <- PeakIndB[PeakIndB < (l1Len - scopeWidth)]
	PeakIndex <- c(PeakIndA, PeakIndB)

	mOnePeak <- length(which(PeakIndex < (l1Len/2)))
	plusOnePeak <- mOnePeak + 1
	lastPeak <- length(PeakIndex)

	##############################################################################################################################################################
	#START FIRST ALIGNMENT#######################################################################################################################################
	##############################################################################################################################################################
	if(prealign)
	{
		mids <- vector()

		for(i in 1:length(results[1,]))
		{
			region <- cbind(results[,i][[2]][low:top], low:top)
			maxInd <- which.is.max(region[,1])
			
			while(maxInd == 1 || maxInd == length(region[,1]))
			{
				if(maxInd == 1)
				{
					region[iniAlignIndexes,1] <- 0
				}
				else
				{
					tinds <- (length(region[,1]) - iniAlignIndexes) + 1
					region[tinds, 1] <- 0
				}
				maxInd <- which.is.max(region[,1])
			}
			
			orreg <- region[rev(order(region[,1])),]
			mids[i] <- orreg[1,2]
		}

		meMids <- median(mids)

		shifts1 <- mids - meMids

		SumMix <- rep(0, splineSize)

		for(i in 1:length(shifts1))
		{
			results[,i][[2]] <- Lag(results[,i][[2]], -shifts1[i])
			results[,i][[2]][is.na(results[,i][[2]])] <- 0
			SumMix <- SumMix + results[,i][[2]]
		}

		fixedMix1 <- SumMix / (length(shifts1))
	}
	else
	{
		fixedMix1 <- rep(0, splineSize)
		for(i in 1:l1Entries)
		{
			fixedMix1 <- fixedMix1 + results[,i][[2]]
		}
		fixedMix1 <- fixedMix1 / l1Entries
		shifts1 <- rep(0, l1Entries)
	}
	
	AllPeak <- vector()
	for(l in 1:mOnePeak) # First peaks from PeakIndex
	{
		lbA <- PeakIndex[l] - scopeWidth
		hbA <- PeakIndex[l] + scopeWidth
		lowA <- floor((lbA) * (baseToSpline))
		topA <- floor((hbA) * (baseToSpline))
	
		peakSamp <- fixedMix1[lowA:topA]
		sortSamp <- peakSamp[rev(order(peakSamp))]
		AllPeak[l] <- sortSamp[1]
	}

	ItMat1 <- matrix(0, l1Entries, length(fixedMix1))

	for(i in 1:l1Entries)
	{
		ItMat1[i,] <- results[,i][[2]]
	}

	currSum <- sum(AllPeak)
	RemedialShifts1 <- vector()
	rowLen <- nrow(MeanSample)

	BaseTestr <- colMeans(ItMat1)

	#START FIRST ITERATIVE ADJUSTMENT
	for(k in 1:rowLen)
	{
		if(k %% 500 == 0 || k == rowLen)
		{
			print(paste("Done ", k, " left alignments", sep=""))
		}
		
		PeakSumList <- vector()
		cnt <- 0
		loShiftRange <- -rangeSize - shifts1[k]
		hiShiftRange <- rangeSize - shifts1[k]

		ToShift <- seq(loShiftRange, hiShiftRange, stepSize)
	
		GeneVec <- (c(rep(0,geneVecAdd), results[,k][[2]], rep(0, geneVecAdd)) / rowLen)
		ItTestr <- BaseTestr - (results[,k][[2]]/rowLen)
	
		for(shif in ToShift)
		{
			cnt <- cnt + 1
			AllPeak <- vector()
		
			for(l in 1:mOnePeak) # First peaks from PeakIndex
			{
				lbB <- PeakIndex[l] - scopeWidth
				hbB <- PeakIndex[l] + scopeWidth
				lowB <- floor((lbB) * (baseToSpline))
				topB <- floor((hbB) * (baseToSpline))
	
				tempMini <- ItTestr[lowB:topB] + GeneVec[(geneVecAdd + lowB - shif):(geneVecAdd + topB-shif)]
				AllPeak[l] <- max(tempMini)
			}
		
			PeakSumList[cnt] <- sum(AllPeak)
		}
		
		slidingPenalty <- 0
		maxInd <- which.is.max(PeakSumList)
		
		while((maxInd <= 1 + slidingPenalty || maxInd >= length(ToShift) - slidingPenalty) && slidingPenalty < (length(ToShift)/2))
		{
			indsToUse <- sideAlignIndexes
			if(slidingPenalty > 0)
			{
				inMAx <- max(indsToUse)
				indsToUse <- c(sideAlignIndexes, (inMAx + 1):(slidingPenalty))
			}
			if(maxInd <= 1 + slidingPenalty)
			{
				PeakSumList[indsToUse] <- 0
			}
			else
			{
				tinds <- (length(ToShift) - indsToUse) + 1
				PeakSumList[tinds] <- 0
			}
			slidingPenalty <- max(indsToUse) + 1
			maxInd <- which.is.max(PeakSumList)
		}
		SortShifts <- ToShift[maxInd]
		SortResults <- max(PeakSumList)
	
		if(SortResults > currSum)
		{
			currSum <- SortResults
			ItMat1[k,] <- Lag(results[,k][[2]], SortShifts)
			ItMat1[k,][is.na(ItMat1[k,])] <- 0
		
			RemedialShifts1[k] <- SortShifts
			BaseTestr <- ItTestr + (ItMat1[k,] / rowLen)
		}
		else
		{
			ItMat1[k,] <- results[,k][[2]]
			RemedialShifts1[k] <- 0
		}
	}

	##############################################################################################################################################################
	#START SECOND ALIGNMENT#######################################################################################################################################
	##############################################################################################################################################################
	if(prealign)
	{
		mids <- vector()

		for(i in 1:length(results2[1,]))
		{
			region <- cbind(results2[,i][[2]][low2:top2], low2:top2)
			
			maxInd <- which.is.max(region[,1])
			while(maxInd == 1 || maxInd == length(region[,1]))
			{
				if(maxInd == 1)
				{
					region[iniAlignIndexes,1] <- 0
				}
				else
				{
					tinds <- (length(region[,1]) - iniAlignIndexes) + 1
					region[tinds, 1] <- 0
				}
				maxInd <- which.is.max(region[,1])
			}
			
			orreg <- region[rev(order(region[,1])),]
			mids[i] <- orreg[1,2]
		}

		#mixByCentre <- mix[order(mids),]

		meMids <- median(mids)

		shifts2 <- mids - meMids

		SumMix <- rep(0, splineSize)

		for(i in 1:length(shifts2))
		{
			results2[,i][[2]] <- Lag(results2[,i][[2]], -shifts2[i])
			results2[,i][[2]][is.na(results2[,i][[2]])] <- 0
			SumMix <- SumMix + results2[,i][[2]]
		}

		fixedMix2 <- SumMix / (length(shifts2))
	}
	else
	{
		fixedMix2 <- rep(0, splineSize)
		for(i in 1:l1Entries)
		{
			fixedMix2 <- fixedMix2 + results[,i][[2]]
		}
		fixedMix2 <- fixedMix2 / l1Entries
		shifts2 <- rep(0, l1Entries)
	}
	
	AllPeak <- vector()
	pcnt <- 0
	for(l in plusOnePeak:lastPeak) # Last seven peaks from PeakIndex
	{
		pcnt <- pcnt + 1
		lbA <- PeakIndex[l] - scopeWidth
		hbA <- PeakIndex[l] + scopeWidth
		lowA <- floor((lbA) * (baseToSpline))
		topA <- floor((hbA) * (baseToSpline))
	
		peakSamp <- fixedMix2[lowA:topA]
		sortSamp <- peakSamp[rev(order(peakSamp))]
		AllPeak[pcnt] <- sortSamp[1]
	}

	ItMat2 <- matrix(0, l1Entries, length(fixedMix2))
	for(i in 1:length(shifts2))
	{
		ItMat2[i,] <- results2[,i][[2]]
	}
	RemedialShifts2 <- vector()
	currSum <- sum(AllPeak)

	BaseTestr <- colMeans(ItMat2)

	#START SECOND ITERATIVE ADJUSTMENT
	for(k in 1:rowLen)
	{
		if(k %% 500 == 0 || k == rowLen)
		{
			print(paste("Done ", k, " right alignments", sep=""))
		}
		PeakSumList <- vector()
		cnt <- 0
		loShiftRange <- -rangeSize - shifts2[k]
		hiShiftRange <- rangeSize - shifts2[k]
		
		ToShift <- seq(loShiftRange, hiShiftRange, stepSize)
	
		ToGetFromIt <- c(1:(k-1), (k+1):rowLen)
	
		if(k == 1)
		{
			ToGetFromIt <- 2:rowLen
		}
		else if(k == rowLen)
		{
			ToGetFromIt <- 1:(rowLen-1)
		}
	
		GeneVec <- (c(rep(0,geneVecAdd), results2[,k][[2]], rep(0, geneVecAdd)) / rowLen)
		ItTestr <- BaseTestr - (results2[,k][[2]]/rowLen)
	
		for(shif in ToShift)
		{
			cnt <- cnt + 1
			AllPeak <- vector()
			pcnt <- 0
		
			for(l in plusOnePeak:lastPeak) # last peaks from PeakIndex
			{
				pcnt <- pcnt + 1
				lbB <- PeakIndex[l] - scopeWidth
				hbB <- PeakIndex[l] + scopeWidth
				lowB <- floor((lbB) * (baseToSpline))
				topB <- floor((hbB) * (baseToSpline))
	
				tempMini <- ItTestr[lowB:topB] + GeneVec[(geneVecAdd + lowB - shif):(geneVecAdd + topB-shif)]
				AllPeak[pcnt] <- max(tempMini)
			}
		
			PeakSumList[cnt] <- sum(AllPeak)
		}
		
		slidingPenalty <- 0
		maxInd <- which.is.max(PeakSumList)
		
		while((maxInd <= 1 + slidingPenalty || maxInd >= length(ToShift) - slidingPenalty) && slidingPenalty < (length(ToShift)/2))
		{
			indsToUse <- sideAlignIndexes
			if(slidingPenalty > 0)
			{
				inMAx <- max(indsToUse)
				indsToUse <- c(sideAlignIndexes, (inMAx + 1):(slidingPenalty))
			}
			if(maxInd <= 1 + slidingPenalty)
			{
				PeakSumList[indsToUse] <- 0
			}
			else
			{
				tinds <- (length(ToShift) - indsToUse) + 1
				PeakSumList[tinds] <- 0
			}
			slidingPenalty <- max(indsToUse) + 1
			maxInd <- which.is.max(PeakSumList)
		}
		SortShifts <- ToShift[maxInd]
		SortResults <- max(PeakSumList)
	
		if(SortResults > currSum)
		{
			currSum <- SortResults
			ItMat2[k,] <- Lag(results2[,k][[2]], SortShifts)
			ItMat2[k,][is.na(ItMat2[k,])] <- 0
		
			RemedialShifts2[k] <- SortShifts
			BaseTestr <- ItTestr + (ItMat2[k,] / rowLen)
		}
		else
		{
			ItMat2[k,] <- results2[,k][[2]]
			RemedialShifts2[k] <- 0
		}

	}

	savS1 <- (RemedialShifts1 + (shifts1))
	savS2 <- (RemedialShifts2 + (shifts2))
	
	retVal <- list()
	retVal[[1]] <- ItMat1
	retVal[[2]] <- ItMat2
	retVal[[3]] <- savS1
	retVal[[4]] <- savS2
	retVal[[5]] <- l1Len
	retVal[[6]] <- binSize
	retVal[[7]] <- splineSize
	
	retVal
}

SaveShifts <- function(Shifts, prefix)
{
	mat <- cbind(Shifts[[3]], Shifts[[4]], rep(Shifts[[7]], length(Shifts[[3]])))
	
	write.table(mat, paste(prefix, ".shifts", sep=""), sep="\t", quote=FALSE, row.names=FALSE)
}

LoadShifts <- function(filename)
{
	AllShift <- as.matrix(read.table(filename, header=TRUE))
	retVal <- list()
	retVal[[1]] <- 0
	retVal[[2]] <- 0
	retVal[[3]] <- as.vector(AllShift[,1])
	retVal[[4]] <- as.vector(AllShift[,2])
	retVal[[5]] <- 0
	retVal[[6]] <- 0
	retVal[[7]] <- AllShift[1,3]
	
	retVal
}

PlotStructMax <- function(MeanSample, Shifts, save=FALSE, prefix="")
{
	splineSize <- length(Shifts[[1]][1,])
	xValsSmall <- seq(-floor(Shifts[[5]]/2)*Shifts[[6]], (floor(Shifts[[5]]/2)*Shifts[[6]]), Shifts[[6]])
	xValsSmall <- xValsSmall[1:length(MeanSample[1,])]

	xValsBigR <- seq(0, (Shifts[[5]]/2)*Shifts[[6]], (Shifts[[5]]*Shifts[[6]])/splineSize)[1:(splineSize/2)]
	xValsBigL <- seq(-(Shifts[[5]]/2)*Shifts[[6]], 0, Shifts[[5]]*Shifts[[6]]/splineSize)[1:(splineSize/2)]

	it1 <- colMeans(Shifts[[1]])
	it2 <- colMeans(Shifts[[2]])
	me1 <- colMeans(MeanSample)

	combiShifts1 <- (Shifts[[3]] / (splineSize/Shifts[[5]])) * Shifts[[6]]
	combiShifts2 <- (Shifts[[4]] / (splineSize/Shifts[[5]])) * Shifts[[6]]

	p1 <- density(combiShifts1)
	p2 <- density(combiShifts2)

	maxRes <- max(it1, it2, me1)
	minRes <- min(it1, it2, me1)
	DmaxY <- max(c(p1$y, p2$y))
	DmaxX <- max(c(p1$x, p2$x))
	DminX <- min(c(p1$x, p2$x))
	
	if(save)
	{
		pdf(paste(prefix, "LRMax.pdf", sep=""), height=10, width=15)
	}
	plot(xValsSmall, me1, type="l", ylim=c(minRes+1,maxRes*1.15), lwd=3, xlab="Base Pair Distance around mean", ylab="Mean Read Counts", main="Chromatin Structural Maximum Alignment, L & R")
	points(xValsBigL, it1[1:(splineSize/2)], type="l", col=rgb(1,0,0,0.6), lwd=4)
	points(xValsBigR, it1[((splineSize/2) + 1):splineSize], type="l", col=rgb(1,0,0,0.6), lwd=2)
	points(xValsBigL, it2[1:(splineSize/2)], type="l", col=rgb(0,1,0,0.6), lwd=2)
	points(xValsBigR, it2[((splineSize/2) + 1):splineSize], type="l", col=rgb(0,1,0,0.6), lwd=4)
	abline(v=0, col="blue", lwd=2)
	legend("topright", c("Aligned on TSS", "Max 5' of TSS", "Max 3' of TSS"), lty=c(1,1,1), lwd=c(3, 3,3), col=c("black", rgb(1,0,0,0.5),rgb(0,1,0,0.5)))

	par(fig=c(0.07,0.37,0.52,0.92), new=TRUE)

	plot(p1, xlim=c(DminX,DmaxX), ylim=c(0, DmaxY), xlab="", ylab="Density", main="TSS position vs Structural Max")
	polygon(p1, col=rgb(1,0,0,0.5), border="blue")
	polygon(p2, col=rgb(0,1,0,0.5), border="blue")
	abline(v=0, col="blue")
	
	if(save)
	{
		dev.off()
	}
}

SplineChange <- function(wt, scount)
{
	x <- 1:ncol(wt)
	wres <- sapply(1:nrow(wt), function(p) spline(x, wt[p,], n=scount, ties="ordered"))
	rmat <- matrix(0, nrow(wt), scount)
	for(i in 1:nrow(wt))
	{
		rmat[i,] <- wres[,i][[2]]
	}
	rmat
}

DoSplineShifts <- function(wt, shift, upscale, txt)
{
	print(paste("Shifting ", txt, sep=""))
	oriLen <- ncol(wt)
	wtTemp <- SplineChange(wt, upscale)
	wtRet <- wtTemp
	for(i in 1:nrow(wt))
	{
		wtRet[i,] <- Lag(wtTemp[i,], shift[i])
		wtRet[i,][is.na(wtRet[i,])] <- 0
	}
	SplineChange(wtRet, oriLen)
}

ApplyStructMax <- function(Shifts, Sample, Right=FALSE)
{
	if(Right)
	{
		shifted <- DoSplineShifts(Sample, Shifts[[4]], Shifts[[7]], "Right..")
	}
	else
	{
		shifted <- DoSplineShifts(Sample, Shifts[[3]], Shifts[[7]], "Left..")
	}
	
	shifted
}

ApplyDualStructMax <- function(Shifts, Sample)
{
	shifted1 <- DoSplineShifts(Sample, Shifts[[4]], Shifts[[7]], "Right..")
	shifted2 <- DoSplineShifts(Sample, Shifts[[3]], Shifts[[7]], "Left..")
	
	len <- length(Sample[1,])
	mid <- floor(len/2)
	print("Merging Left and Right..")
	shifted <- cbind(shifted2[,1:mid], shifted1[,(mid+1):len])
	
	shifted
}

FindPeaks <- function(MeanSamples)
{
	Peaks <- colMeans(MeanSamples)
	cnt = 0
	PeakIndex <- vector()

	for(i in 5:(length(Peaks)-5))
	{
		if(Peaks[i-4] < Peaks[i] && Peaks[i-3] < Peaks[i] && Peaks[i-2] < Peaks[i] && Peaks[i-1] < Peaks[i])
		{
			if(Peaks[i+4] < Peaks[i] && Peaks[i+3] < Peaks[i] && Peaks[i+2] < Peaks[i] && Peaks[i+1] < Peaks[i])
			{
				cnt <- cnt + 1
				PeakIndex[cnt] <- i
			}
		}	
	}
	
	mid <- floor(length(Peaks)/2)
	minusMid <- PeakIndex - mid
	ChosenPeaks <- vector()
	cnt <- 0
	while(length(ChosenPeaks) < 4 && cnt < mid/2)
	{
		ChosenPeaks <- minusMid[abs(minusMid) < (10 + cnt)]
		cnt <- cnt + 1
	}
	
	ch <- ChosenPeaks + mid
	ch
}

