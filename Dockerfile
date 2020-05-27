FROM rocker/rstudio:3.6.0
ENV RENV_VERSION 0.10.0-26
RUN mkdir home/analysis
WORKDIR /home/analysis

COPY renv.lock /home/analysis

RUN apt-get -y update
RUN apt-get -y upgrade
RUN apt -y install libz-dev
RUN apt -y install libpng-dev
RUN apt -y install libjpeg-dev
RUN apt -y install libmagick++-dev

RUN R -e "install.packages('remotes',repos = c(CRAN = 'https://cloud.r-project.org'))"
RUN R -e "remotes::install_github('rstudio/renv@${RENV_VERSION}')"
RUN R -e "renv::restore()"

COPY ./R/ /home/analysis/R

