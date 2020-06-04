FROM rocker/rstudio:3.6.0
ENV RENV_VERSION 0.10.0-26
ENV RENV_PATHS_CACHE /main/renv/cache
#RUN mkdir home/analysis
WORKDIR /main/

COPY renv.lock /main/

#RUN apt-get -y update 
#RUN apt -y install zlib1g-dev
#RUN apt -y install libpng-dev
#RUN apt -y install libjpeg-dev
#RUN apt -y install libmagick++-dev

RUN echo "options(renv.consent = TRUE)" >> .Rprofile
RUN ls /main/
RUN ls /main/renv/cache/
RUN ls /main/
RUN R -e "install.packages('remotes',repos = c(CRAN = 'https://cloud.r-project.org'))"
RUN R -e "remotes::install_github('rstudio/renv@${RENV_VERSION}')"
RUN R -e "Sys.setenv(RENV_PATHS_CACHE = '/main/renv/cache')"
RUN R -e "print(renv:::renv_paths_cache())"
#RUN R -e "renv::restore(confirm=FALSE)"
#RUN R -e "setwd('/home/analysis')"

COPY ./R/ /home/analysis/R

