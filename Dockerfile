FROM rocker/rstudio:3.6.0
ENV RENV_VERSION 0.10.0-26
RUN mkdir home/analysis
WORKDIR /home/analysis

COPY . home/analysis

RUN R -e "install.packages('remotes',repos = c(CRAN = 'https://cloud.r-project.org'))"
RUN R -e "remotes::install_github('rstudio/renv@${RENV_VERSION}')"
RUN R -e "renv::restore()"


