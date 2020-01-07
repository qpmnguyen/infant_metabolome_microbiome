FROM rocker/r-ver:3.6.0
ENV RENV_VERSION 0.9.2-16
RUN apt-get update && \
	apt-get install tree && \
	apt-get install wget && \
	wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
	bash Miniconda3-latest-Linux-x86_64.sh && \
	conda install -c bioconda snakemake 
RUN mkdir /home/analysis/
COPY renv.lock /home/analysis/renv.lock
COPY R/ /home/analysis/R/
COPY snakefiles/ /home/analysis/snakefiles 
COPY data/ /home/analysis/data
RUN R -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org'))"
RUN R -e "remotes::install_github('rstudio/renv@${RENV_VERSION}')"
RUN R -e 'renv::restore()'