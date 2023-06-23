#FROM miseyu/docker-ubuntu16-python3.6:latest
FROM ubuntu:16.04
MAINTAINER xufeng <crazyhsu9527@gmail.com>

# install basic dependencies
#RUN sed -i s@/archive.ubuntu.com/@/mirrors.tuna.tsinghua.edu.cn/@g /etc/apt/sources.list
RUN apt-get clean && apt-get update && apt-get install -y --no-install-recommends \
	build-essential \
	cmake \
	git \
	wget \
	apt-transport-https \
	libz-dev \
	gcc \
	ca-certificates \
	vim \
	less \
	sudo \
	apt-utils

#RUN python3.6 -m pip install pip --upgrade && \
#    python3.6 -m pip install wheel

#pip3 install pip==19.3.1 numpy==1.16 pkgconfig

ARG UID=1001
ARG GID=1001
RUN groupadd --gid $GID iflas && useradd --uid $UID --gid iflas --shell /bin/bash --create-home iflas
RUN echo "root:123456" | chpasswd
USER iflas

ENV IFLAS_HOME /home/iflas
WORKDIR $IFLAS_HOME

RUN mkdir $IFLAS_HOME/data && mkdir $IFLAS_HOME/software

# install miniconda for python2.7
RUN wget --no-check-certificate https://mirrors.tuna.tsinghua.edu.cn/anaconda/miniconda/Miniconda2-latest-Linux-x86_64.sh -O ~/miniconda2.sh
#RUN wget --no-check-certificate https://repo.anaconda.com/miniconda/Miniconda2-latest-Linux-x86_64.sh -O ~/miniconda2.sh
RUN /bin/bash ~/miniconda2.sh -b -p $IFLAS_HOME/software/conda

ENV PATH=$IFLAS_HOME/software/conda/bin:$IFLAS_HOME/software:$PATH
RUN conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free/ && \
	conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main/ && \
	conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/ && \
	conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda/ && \
	conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/r/
#
# RUN conda update conda -y
#
# # create isoseq3 environment with python2.7
RUN conda create -n isoseq3_pipeline python=2.7

# set some environment variables
ENV CONDA_DEFAULT_ENV=isoseq3_pipeline
ENV CONDA_PREFIX=$IFLAS_HOME/software/conda/envs/$CONDA_DEFAULT_ENV
ENV CONDA_AUTO_UPDATE_CONDA=false

RUN echo ". $IFLAS_HOME/software/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
	echo "conda activate isoseq3_pipeline" >> ~/.bashrc

RUN conda install -y -c conda-forge r-base=3.6.3
RUN conda install -y -c conda-forge numpy=1.16.5 pandas=0.24.2 matplotlib=2.2.3 scipy=1.2.1 scikit-image=0.14.3 psutil pulp networkx=2.2
RUN conda install -y -c bioconda biopython pysam blast subread samtools=1.9 pybedtools stringtie minimap2 seqkit hisat2 mash pyvcf
RUN conda install -y -c bioconda ucsc-bedtogenepred ucsc-gtftogenepred ucsc-genepredtogtf
RUN conda install -y -c bioconda isoseq3=3.3 pbccs=4.2 rmats=4.0.2 bioconductor-deseq2 bioconductor-clusterprofiler=3.14.0
RUN conda install -y r-ggplot2 r-stringi=1.4.6
RUN R -e "install.packages('valr', repos='https://mirrors.tuna.tsinghua.edu.cn/CRAN/')" && R -e "install.packages('tidyr', repos='https://mirrors.tuna.tsinghua.edu.cn/CRAN/')"


conda install -y -c bioconda samtools=1.9 hisat2 bedtools=2.29.2 subread stringtie minimap2 seqkit fastp
conda install -y -c bioconda ucsc-gtftogenepred ucsc-genepredtogtf fmlrc2=0.1.4 nanopolish=0.11.1 regtools=0.5.2
conda install -y -c bioconda isoseq3=3.3 pbccs=4.2 lima pbcoretools rmats=4.0.2 bax2bam pbbam=1.0.6 pbcopper=1.3.0
conda install -y -c conda-forge r-base=3.6.3 rpy2=2.8.6 r-scales
conda install -y -c bioconda bioconductor-deseq2=1.26.0 bioconductor-clusterprofiler=3.14.0 bioconductor-gviz=1.30.0
conda install -y bx-python==0.7.3
conda install -y r-stringi=1.4.6 r-dplyr=1.0.0 r-tibble=3.0.0 r-gridBase

pip install pandas matplotlib==2.2.3 psutil biopython==1.68 pybedtools==0.8.2 PyVCF PyPDF2 bx-python==0.7.3

#pysam==0.15.3 numpy==1.16.5
pip install pandas==0.24.2 matplotlib==2.2.3 psutil biopython==1.68 pybedtools==0.8.2 PyVCF==0.6.8 PyPDF2 bx-python==0.7.3 pulp networkx==2.2 scipy==1.2.1 scikit-image==0.14.3 scikit-learn==0.20.4 yattag==1.14.0 PyMuPDF==1.18.0

# generic installation
RUN git clone https://github.com/Magdoll/cDNA_Cupcake.git $IFLAS_HOME/software/cDNA_Cupcake && \
    cd $IFLAS_HOME/software/cDNA_Cupcake && git checkout Py2_v8.7.x && \
    $CONDA_PREFIX/bin/python setup.py build && $CONDA_PREFIX/bin/python setup.py install

RUN git clone https://hub.fastgit.org/CrazyHsu/SpliceGrapher_packages.git $IFLAS_HOME/software/SpliceGrapher_packages
RUN cd $IFLAS_HOME/software/SpliceGrapher_packages && tar -xf PyML-0.7.14.tar.gz && rm -rf PyML-0.7.14.tar.gz
WORKDIR $IFLAS_HOME/software/SpliceGrapher_packages/PyML-0.7.14
RUN $CONDA_PREFIX/bin/python setup.py build && $CONDA_PREFIX/bin/python setup.py install

RUN cd $IFLAS_HOME/software/SpliceGrapher_packages && tar -xf SpliceGrapher-0.2.7.tgz && rm -rf SpliceGrapher-0.2.7.tgz
WORKDIR $IFLAS_HOME/software/SpliceGrapher_packages/SpliceGrapher-0.2.7
RUN $CONDA_PREFIX/bin/python setup.py build && $CONDA_PREFIX/bin/python setup.py install

WORKDIR $IFLAS_HOME