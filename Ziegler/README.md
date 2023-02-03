## Prerequisites

The scripts work with following 'expectations':

**basedir** is the basic data / output directory. It should have *RawData*, *FilteredData* and *OutputData* subfolders

**codedir** is the directory containing all code-relevant parts. *git clone* the repository and supply it as the code dir. To be save make sure to end the string with '/'

**Ziegler data** as provided by Ziegler et. al. The data should be place in the RawData of the *basedir*. Following tables are being used: 
    - 20210220_NasalSwab_NormCounts.txt
    - 0210220_NasalSwab_RawCounts.txt
    - 20210701_NasalSwab_MetaData.txt
    - 20210220_NasalSwab_UMAP.txt
**Yoshida data** as provided by on cell atlas. The data should be extracted and place in the RawData of the *basedir*. Only airway data is needed, don't put it into a separate subfolder.  Following tables are being used: 
    - X.mtx
    - obs.csv
    - var.csv



## Ziegler dataset 

### 01_DataPreprocessing.R

**./01_DataPreprocessing.R** *basedir* *codedir*

The output of the script is put into the FilteredData directory.

### 02_ExpressionAnalysis_TcellGeneExpression.R

**./02_ExpressionAnalysis_TcellGeneExpression.R** *basedir* *codedir*

The output of the script is put into OutputData/markerGenes/ directory

### 03_ExpressionAnalysis_IFNPosVSneg.R

**./03_ExpressionAnalysis_IFNPosVSneg.R** *basedir* *codedir*

The output of the script is put into OutputData/markerGenes/ directory

### 04_PathwayAndGOAnalysis_IFNPosVSneg.R

**./04_PathwayAndGOAnalysis_IFNPosVSneg.R** *basedir* *codedir*

The output of the script is put into OutputData/pathwayEnrichment/ directory

### 05_CombineAnalyses.R

This script should be executed after the analysis of Yoshida et. al. dataset was made, as it relies on the output of the analysis

**./05_CombineAnalyses.R** *basedirZiegler* *basedirYoshida* *codedir*

The output of the script is put into OutputData/ directory of *basedirZiegler*

### 06_Figure5.R

**./06_Figure5.R** *basedirZiegler* *basedirYoshida* *codedir*

The full figure is put into *basedirZiegler* OutputData/ directory. Subfigures are in the folders of the respecive dataset.

### 07_Supplementary_Fig6.R

**./07_Supplementary_Fig6.R** *basedirZiegler* *basedirZiegler* *codedir*

The full figure is put into *basedirZiegler* OutputData/ directory. Subfigures are in the folders of the respecive dataset.
