# Olivers-R-Scripts
Post-processing and stats for chromatin structures

**Using The Scripts**

At the moment these scripts have not been merged into a formal package.
As I create generalized versions of the scripts I will upload them here.
Please let me know if you find any bugs! (o.j.rimington@gmail.com)

To use them, download the relevant scripts to you your working directory in R.
Then: source("ScriptName.R")

Specific Functions are described below.

___

## Structural Max Alignments

Script: StructureMax.R

R Pre-Reqs: Hmisc, nnet
&nbsp;

&nbsp;

####`StructureMax(MeanSample, PeakL, PeakR, ..)`

Iteratively adjusts the alignment of each chromatin region to obtain the most-structured representation of cumulative nucleosome positions. Each region is first spline-interpolated to a much higher resolution, a range of potential shifts are then tested on that region.
&nbsp;
* MeanSample	-	This is a matrix of levels for feature-oriented regions of chromatin 

* PeakL			-	Integer. This is the user's chosen index for the -1 Peak (Use `FindPeaks` if uncertain)

* PeakR			-	Integer. This is the user's chosen index for the +1 Peak (Use `FindPeaks` if uncertain)

* prealign		-	TRUE/FALSE, Whether or not to pre-align on a single peak before. Only use without if chromatin is already well ordered. **Def=TRUE** 

* binSize		-	Integer, The number of base pairs in the binning method. **Def=10**

* stepSize		- 	Integer. Effects quantity and resolution of the alignment shift attempts, lower=better/more expensive. **Def=5**

* splineSize	-	Integer. Resolution of the spline interpolation. higher=better/more memory. **Def=10000**


Use:

```
myShifts = StuctureMax(meansMatrix, 144, 156, splineSize=5000)
```

Output:

'myShifts' is an object containing the two vectors of alignment shifts, and two matrices of interpolated and shifted regions. One is shifted to left hand maximum, the other shifted to the right hand maximum. This object also contains some contextual information relating to the alignment process which can be used in downstream functions.

The 'myShifts' object is then passed to the following functions:

&nbsp;

&nbsp;


####`PlotStructureMax(MeanSample, Shifts, ..)`

Plots the results of 

* MeanSample	-	The original matrix supplied to `StructureMax`

* Shifts		-	The object returned by the `StructureMax` function

* save			-	TRUE/FALSE, Whether or not to save to file. **Def=FALSE**

* prefix		-	String. Prefix for the saved pdf filename. Format "prefixLRshifts.pdf". **Def=""**


Use:

```
PlotStructureMax(meansMatrix, myShifts, TRUE, "test1_")
```

Output:

A graph of the structural maximum alignments, including the pre-aligned trace. This also shows the relationship between the structural maximum and the central point of the region of interest. Usually this is the 'Transcription Start Site'.

&nbsp;

&nbsp;

####`ApplyStructMax(Shifts, Sample, ..)`

Applies the alignments discovered in `StructureMax` to another sample.

* Shifts		-	The output object from `StructureMax`

* Sample		-	The matrix of feature-oriented chromatin structure regions for a single sample

* Right			-	TRUE/FALSE, flag determines whether the alignement takes place on the left or right hand structural maximum. **Def=FALSE**

Use:

```
Samp1Left = ApplyStructMax(myShifts, Samp1, FALSE)
``` 

Output:

The left or right struct. max aligned version of a sample.

&nbsp;

&nbsp;

####`ApplyDualStructMax(Shifts, Sample)`

Applies both left and right hand structural maximum alignments, and merges the two region-sets at the mid-point. Useful for doing bootstrapping stats tests on the full aligned structures.

* Shifts		-	The output object from `StructureMax`

* Sample		-	The matrix of feature-oriented chromatin structure regions for a single sample

Use:

```
Samp1Both = ApplyDualStructMax(myShifts, Samp1)
```

Output:

A matrix, equal in dimensions to 'Sample', containing the left half of the left-hand structural maximum, and the right half of the right-hand structural maximum

&nbsp;

&nbsp;

####`SaveShifts(Shifts, prefix)`

Utility function to save the alignment shifts object to file, to be re-loaded later. Shifts saved and re-loaded in the way CANNOT be used with `PlotStructureMax`

* Shifts		-	The object output from `StructureMax`

* prefix		-	A prefix to the saved file, can include directory address too. File name format is 'prefix.shifts'

Use:

```
SaveShifts(myShifts, "testSamp")
```

Output:

A file containing alignment shift data that can be loaded later.

&nbsp;

&nbsp;


####`LoadShifts(filename)`

Utility function to re-load the alignment shifts from file. Shifts saved and re-loaded in the way CANNOT be used with `PlotStructureMax`

* filename		-	Path to, and name of, a file saved using `SaveShifts`. This will end in '.shifts'.

Use:

```
myShifts = LoadShifts(myShifts, "Some_Folder/testSamp.shifts")
```

Output:

A Shifts object, that can be used with `ApplyStructMax` and `ApplyDualStructMax`

&nbsp;

&nbsp;


####`FindPeaks(meanSample)`

Utility function to assist in the identification of -1 and +1 peaks for structural maximum alignment. Searches for peaks around the centre of the region.

* meanSample	-	The matrix to be used for `StructureMax`

Use:

```
FindPeaks(myShifts)
```

Output:

A short vector containing the indices of peaks found near the centre of the region. May be helpful in identifying -1 and +1 peak indices for `StructureMax`
&nbsp;

&nbsp;

&nbsp;

___

##Individual Models + Rank Groups and Clustering functions
Script: IndividualModels.R

R Pre-Reqs: plotrix, boot, RColorBrewer, plotrix, fields
&nbsp;

&nbsp;


####`MakeIndividualModels(binSize=10, ...)`

Primary function that Builds a model of variation for each bin from replicates, and processes the statistical significance of the variation found in a single mutant sample. This can take several hours. The efficiency is O(n^2), where n is the number of regions. It is suggested that StructureMax may be used to create optimally structured alignments of regions before running this function.

* binSize		-	The size of the basepair binning system used

* ...			-	2 or more matrices of chromatin levels around a feature of interest. THE MUTANT IS ALWAYS THE LAST IN THE LIST!

Usage:

```
results = MakeIndividualModels(10, Rep1, Rep2, Rep3, Mutant)
```

Output:

A list of output matrices with data associated with the model building process. Within that list:

* results[[1]] - a matrix of distribution means used for tests, equal in dimension to the input matrices
* results[[2]] - like [[1]], except standard deviations
* results[[3]] - a list of matrices representing the -log p-values for each sample tested under the model. In the same order as they were entered into the function.
* results[[4]] - a single value (the bin size)
* results[[5]] - the final matrix of corrected p-values for the mutant sample

This list is then passed to the following functions to visualise and test the model results.

&nbsp;

&nbsp;

####`PlotIndividualModelTest(results, save, prefix)`

This function creates a visualisation of the p-values achieved by all regions when tested under the model. Each graph shows half of the entire vertically sorted set of p-values for each sample. This is useful to look at, it shows which regions had the highest variation in the mutant, and how well the modelling process worked for the given replicates. 

* results		-	The list object produced by `MakeIndividualModels`
* save			-	Whether or not to save the graph to a pdf. **Def=FALSE**
* prefix		-	Filename prefix if save option is used. **Def=""**

Usage:

```
PlotIndividualModelTest(results, TRUE, "filename/indi1")
```

Output:

A graph that shows each of the input samples tested under the model. In the case of the replicates they are being tested under their own model for the purposes of evaluating the suitability of the model, and more importantly, to act as a reference for the level of significance achieved by the mutant sample.

&nbsp;

&nbsp;
####`PlotMultiTestCorrection(results, save, prefix)`

This function shows the level of significance achieved by the entire mutant sample when the maximum significance of the re-tested replicates is used as a multiple testing correction. 

* results		-	The list object produced by `MakeIndividualModels`
* save			-	Whether or not to save the graph to a pdf. **Def=FALSE**
* prefix		-	Filename prefix if save option is used. **Def="mtcorrect"**

Usage:

```
PlotMultiTestCorrection(results, TRUE, "filename/multi1")
```

Output:

A graph showing three distributions of p-values accross the feature region. These are the Replicate Maximum, the Mutant, and the Corrected Mutant. The replicate maximum is essentially the set of levels that are used to correct the mutant significances.

&nbsp;

&nbsp;
####`PlotAllRegionGraphs(prefix, results, ...)`

This function creates graphs for all the regions used in the samples. This is a two-part plot including the levels for each sample, and the sigificance of the mutant deviation. It is highly reccommended to create a new folder for the output. 

* prefix		-	Filename prefix, can include directories
* results		-	The list object produced by `MakeIndividualModels`
* ...			-	2 or more matrices of chromatin levels around a feature of interest. Must be the same as used in `MakeIndividualModels`	

Usage:

```
PlotAllRegionGraphsfunction("Regions/", results, Rep1, Rep2, Rep3, Mutant)
```

Output:

One graph per region. If there are 4000 regions, there will be 4000 graphs - it is reccommended to make a separate directory. This allows closer investigation of the significance levels for any given region. Particularly useful if there are a few regions of interest which bear closer examination. Your input matrix row names will be used for file names, and in the graph titles.

&nbsp;

&nbsp;
####`PlotTopRankMutant(indexes, division, results, ...)`

This function selects a specific range of base-pairs around a feature and ranks regions based on the significance of the changes found in the mutant sample. The top 'rank-group' is then displayed as a graph. The number of rank-groups is determined by 'division'. If it is 5, then the top rank group that is displayed will be the top 20% of the set for the given ranking process. If it is 10, then it will be 10%.
All graphs created by this function are automatically processed via a bootstrap signifiance testing method.

* indexes		-	two integer vector describing the indexes between which to rank the mutant change signifiance
* division		-	single integer, the number by which the set will be divided. Size of top rank group will equal ((Number of regions)/division)
* results		-	The list object produced by `MakeIndividualModels`	
* ...			-	2 or more matrices of chromatin levels around a feature of interest. Must be the same as used in `MakeIndividualModels`	

Usage:

```
PlotTopRankMutant(c(1,120), 10, results, Rep1, Rep2, Rep3, Mutant)
```

Output:

A graph showing the mean levels of most significantly different mutant regions, with the significance of the cumulative differences assessed by a bootstrap test. Vertical bars indicate the area, defined by 'indexes', within which the mutant region p-values were used to rank-order the set. There is no option to automatically save this graph - this can be done using `PlotAllRankGraphs`.

&nbsp;

&nbsp;
####`PlotAllRankGraphs(indexes, division, FilePrefix, results, ...)`

This function is the full version of `PlotTopRankMutant` - the functionality is the same, it's purpose is to generate all rank-group graphs. These are useful to compare the difference between the lowest and highest rank groups.
All graphs created by this function are automatically processed via a bootstrap signifiance testing method.

* indexes		-	two integer vector describing the indexes between which to rank the mutant change signifiance
* division		-	single integer, the number by which the set will be divided. Size of top rank group will equal ((Number of regions)/division)
* results		-	The list object produced by `MakeIndividualModels`	
* FilePrefix	-	The prefix of the generated graph files. Can include directories.
* ...			-	2 or more matrices of chromatin levels around a feature of interest. Must be the same as used in `MakeIndividualModels`

Usage:

```
PlotAllRankGraphs(c(100:170), 10, "RankFile/test_", results, Rep1, Rep2, Rep3, Mutant)
```

Output:

A set of pdfs, equal in number to 'division'. These are in the same format as the graph produced by `PlotTopRankMutant`, except all rank-groups are processed.

&nbsp;

&nbsp;

####`PlotSumSquares(results, MaxNumber)`

This is a utility function to assist in the cluster number determination process for kmeans clustering. The user ought to look for the 'knee' on the graph - this is the number which theoretically represents the optimally structured division of the set in terms of cluster counts. Clusters are created based on regions that share similar patterns of mutant variation.

* results		-	The list object produced by `MakeIndividualModels`
* MaxNumber		-	The highest k value to try (number of clusters) **Def=15**

Usage:

```
PlotSumSquares(results, 10)
```

Output:

A Sum-of-squares graph. Each point is described as x = k (number of clusters), y = Sum of squares within clusters. When k increases, and a large cluster splits evenly into two, (for example), the sum of squares will reduce significantly, creating a 'knee'. If the increase in k only results in another small peripheral cluster, there will be no signifiance drop in the sum of squares.

&nbsp;

&nbsp;
####`ClusterRegions(results, clusters)`

Once the k value has been determined with `PlotSumSquares` this simple wrapper will run the kmeans function and return a 'kmeans' object. 

* results		-	The list object produced by `MakeIndividualModels`
* clusters		-	The number of kmeans clusters to be created.

Usage:

```
kclus = ClusterRegions(results, 6)
```

Output:

A 'kmeans' object. This is the standard output created by R's kmeans function. It contains a cluster index that can be used to split up the sample matrices, and the p-value results matrices.

&nbsp;

&nbsp;
####`PlotOneCluster(kclusters, ClusterNumber, ...)`

This function plots the cumulative levels of chromatin structure from one specific cluster - determined by the `ClusterRegions` function. This should show one particular type of mutant variation - if there is enough structure to the set of mutant variation to justify clustering.
All graphs created by this function are automatically processed via a bootstrap signifiance testing method.

* kclusters		-	The 'kmeans' object returnd by `ClusterRegions`
* ClusterNumer	-	The index of the cluster that you wish to plot.
* ...			-	2 or more matrices of chromatin levels around a feature of interest. Must be the same as used in `MakeIndividualModels`

Usage:

```
PlotOneCluster(kclus, 3, Rep1, Rep2, Rep3, Mutant)
```

Output:

A single graph showing the cumulative variation of the mutant in the selected cluster.

&nbsp;

&nbsp;
####`PlotAllClusters(kclusters, prefix, ...)`

This is the full version of `PlotOneCluster`. It generated comparison graphs for all clusters created in `ClusterRegions`.
All graphs created by this function are automatically processed via a bootstrap signifiance testing method.

* kclusters		-	The 'kmeans' object returnd by `ClusterRegions`
* prefix		-	The prefix of the generated graph files. Can include directories.**Def=""**
* ...			-	2 or more matrices of chromatin levels around a feature of interest. Must be the same as used in `MakeIndividualModels`

Usage:

```
PlotAllClusters(kclus, "ClusGraphs/clus", Rep1, Rep2, Rep3, Mutant)
```

Output:

Graphs showing the cumulative variation of the mutant in each clusters, with the corresponding signifiance test for those variations. All are saved as pdfs.

&nbsp;

&nbsp;
####`SaveIndividualModel(results, prefix)`

Since the time and computing requirements for running `MakeIndividualModels` can be quite high, it is recommended that this function be run to save the results in the case of very large region sets (10,000+).

* results		-	The list object produced by `MakeIndividualModels`
* prefix		-	The prefix for the saved file, can include directories. The file will have the suffix '.indiv'. **Def="model"**

Usage:

```
SaveIndividualModel(results, "Outfile/myModel")
```

Output:

A single text file containing all the modelling results - can be reloaded.

&nbsp;

&nbsp;
####`LoadIndividualModel(fileName)`

This function can load the model output back into R, so long as it was saved with `SaveIndividualModel`.

* fileName		-	The name of the results file to be loaded.

Usage:

```
results = LoadIndividualModel("Outfile/myModel.indiv")
```

Output:

A model output list object, identical to the one produced by `MakeIndividualModels`.

&nbsp;

&nbsp;
