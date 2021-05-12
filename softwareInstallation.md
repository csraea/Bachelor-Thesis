### **Installing Necessary Libraries & Software**

This thesis covers the development of the tool aimed at visualizing different genome properties. The provided guide describes environment adjusting. 
Firstly, the Anaconda distribution is required which can be obtained at https://www.anaconda.com/download. 
The next steps are to install all necessary libraries and packages in the following order.
1) Create a new conda environment called bioinformatics with biopython=1.70:

        conda create -n bioinformatics biopython biopython=1.70
1) Activate the environment:

        source activate bioinformatics
    or

        conda activate bioinformatics
1) Add the _bioconda_ and _conda-forge_ channel to our source list:

        conda config --add channels bioconda
        conda config --add channels conda-forge
1) Install the core packages:

        conda install IPython  scipy matplotlib jupyter-notebook pip pandas
        cython numba scikit-learn seaborn pysam pyvcf simuPOP dendropy rpy2=2.0.x

1) Install _R_ using Anaconda:

        conda install r-essentials r-gridextra r-caret r-lazyeval