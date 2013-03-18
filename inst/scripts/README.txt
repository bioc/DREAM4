Ten serialized SummarizedExperiment objects are provided with the
DREAM4 package, each of which has 

  * synthetic expression data for a variety of simulated 
    experiments
  * a gold standard matrix

and each of which is assembled from these files

  * dualknockouts.tsv
  * goldStandard.tsv
  * knockdowns.tsv
  * knockouts.tsv
  * multifactorial.tsv
  * timeseries.tsv
  * wildtype.tsv

All ten experiments have seven files, found in the 
extdata/lightlyProcessedDownloadedData.tar.gz.  

These files have indeed been "lightly processed" from the versions
offered by the DREAM project, for which this is the crucial
information:

  http://wiki.c2b2.columbia.edu/dream/data/DREAM4/
  team name: Biocondcutor
  password: 55IQb6h7

 	DREAM4_InSilico_Size10.zip 	37828 	421 	application/zip
	DREAM4_InSilico_Size100.zip 	890813 	352 	application/zip
	DREAM4_InSilico_Size100_Multifactorial.zip 	216787 	317 	application/zip

Gold standard matrices are allegedly available at

   http://wiki.c2b2.columbia.edu/%3C?=$wgScriptName?%3E/data/gold-standards/DREAM4/

but on January 31st 2013, this link is broken.  

The "light processing" creates uniform directories with
identically-named files, for example, for insilico_size10_1/

  dualknockouts.tsv
  goldStandard.tsv
  knockdowns.tsv
  knockouts.tsv
  multifactorial.tsv
  timeseries.tsv
  wildtype.tsv

This light processing is motivated by the simple processing it makes
possible.  File names are made uniform for each directory, and the
separately obtained gold standard matrix files are included with each
directory.

To rebuild the ten serialized SummarizedExperiment objects of this 
package, take these steps:

  * cd <workingDirectory>
  * tar zxf <DREAM4_packagedirectory>/inst/extdata/lightlyProcessedDownloadedData.tar.gz
  * cd lightlyProcessedDownloadedData
  * start R
  * source("<DREAM4_packagedirectory>/inst/scripts/buildRData.R")
  * runTests() 
  * run(dataDir=".", destDir=".")



