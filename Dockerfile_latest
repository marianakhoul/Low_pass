FROM r-base:latest

RUN apt-get update
RUN apt-get -y install curl libcurl4-gnutls-dev apt-transport-https \
  ca-certificates gnupg lsb-release libxml2-dev

RUN R --vanilla -e 'install.packages("BiocManager", repos = "http://cran.us.r-project.org")'

# needs curl
RUN R --vanilla -e 'BiocManager::install("QDNAseq")'
RUN R --vanilla -e 'BiocManager::install("QDNAseq.hg19")'
RUN R --vanilla -e 'BiocManager::install("BSgenome")'
RUN R --vanilla -e 'BiocManager::install("Biobase")'
RUN R --vanilla -e 'BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")'
RUN R --vanilla -e 'install.packages("optparse", repos = "http://cran.us.r-project.org")'
RUN R --vanilla -e 'install.packages("ggplot2", repos = "http://cran.us.r-project.org")'
RUN R --vanilla -e 'install.packages("dplyr", repos = "http://cran.us.r-project.org")'
RUN R --vanilla -e 'install.packages("future", repos = "http://cran.us.r-project.org")'

COPY low_pass_cna_tool.R /low_pass_cna_tool.R

RUN echo "alias ulp-tool='Rscript --vanilla /low_pass_cna_tool.R'" >> /root/.bashrc

