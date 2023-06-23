conda create -y -n iflas_toolkit_20210519 python=2.7
conda activate iflas_toolkit_20210519

#conda install -y -c bioconda perl-bioperl perl-app-cpanminus perl-dbi
conda install -y -c bioconda samtools=1.9 hisat2 bedtools=2.29.2 bamtools subread stringtie minimap2 seqkit fastp
#conda install -y -c bioconda subread stringtie minimap2 seqkit fastp
#conda install -y -c bioconda perl-bioperl perl-app-cpanminus perl-dbi subread samtools=1.9 stringtie \
#                  minimap2 seqkit hisat2 bedtools fastp bamtools
conda install -y -c bioconda ucsc-gtftogenepred ucsc-genepredtogtf fmlrc2=0.1.4 nanopolish=0.11.1 regtools=0.5.2
conda install -y -c bioconda isoseq3=3.3 pbccs=4.2 lima pbcoretools rmats=4.0.2 bax2bam pbbam=1.0.6 pbcopper=1.3.0
conda install -y -c conda-forge r-base=3.6.3 rpy2=2.8.6 r-scales
conda install -y -c bioconda bioconductor-deseq2=1.26.0 bioconductor-clusterprofiler=3.14.0 bioconductor-gviz=1.30.0
#conda install -y -c bcbio bx-python
conda install -y r-stringi=1.4.6 r-dplyr=1.0.0 r-tibble=3.0.0 r-gridBase

pip install pandas matplotlib==2.2.3 psutil biopython==1.68 pybedtools PyVCF PyPDF2 bx-python==0.7.3

R -e "install.packages('tidyr', repos='https://mirrors.tuna.tsinghua.edu.cn/CRAN/')"
R -e "install.packages('valr', repos='https://mirrors.tuna.tsinghua.edu.cn/CRAN/')"
#R -e "install.packages(c('gapminder','ggplot2','scales','reshape', 'gridBase','gridExtra','grid','rpart','ROCR','caret','optparse','lattice','foreach','e1071','randomForest','partykit','ipred','rpart.plot','doMC','nnet','ROSE','pROC','MLmetrics'), repos = 'https://mirrors.tuna.tsinghua.edu.cn/CRAN/')"


git clone https://hub.fastgit.org/CrazyHsu/cDNA_Cupcake.git
cd cDNA_Cupcake
git checkout Py2_v8.7.x
python setup.py build && python setup.py install
cd ../

git clone https://hub.fastgit.org/CrazyHsu/SpliceGrapher_packages.git
cd SpliceGrapher_packages
tar -xf PyML-0.7.14.tar.gz
cd PyML-0.7.14
python setup.py build && python setup.py install
cd ../
tar -xf SpliceGrapher-0.2.7.tgz
cd SpliceGrapher-0.2.7
python setup.py build && python setup.py install
cd ../
cd ../