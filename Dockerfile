FROM ubuntu
ENV RENV_VERSION 0.9.2-16
ARG DATA="./data/"
RUN apt-get update && \
	apt-get --assume-yes install apt-transport-https software-properties-common && \
	apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 && \
	add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/' && \
	apt-get update && \ 
	apt-get --assume-yes install wget git && \
	apt-get --assume-yes install r-base && \
	apt-get install -y wget git tabix vcftools && \
	apt-get install python3-pip && \
	pip3 install snakemake
RUN mkdir /home/analysis/
WORKDIR /home/analysis
COPY renv.lock renv.lock
COPY R/ R/
COPY ${DATA} data/
COPY snakefiles/ snakefiles/
RUN R -e "setwd("home/analysis/")
RUN R -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org'))"
RUN R -e "remotes::install_github('rstudio/renv@${RENV_VERSION}')"
RUN R -e 'renv::restore()'