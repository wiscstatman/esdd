#!/bin/bash

# untar your R installation
tar -xzf R.tar.gz

# make sure the script will use your R installation
export PATH=$(pwd)/R/bin:$PATH
export RHOME=$(pwd)/R

# run R, with the name of your  R script
Rscript condor.R $1 $2
