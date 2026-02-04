# spec_analysis
Spectral analysis of galaxy clusters

## Abstract

This is the automation of analyzing multiple spectra of galaxy clusters from X-ray observations made by the Chandra X-Ray Observatory. 

To best utilize the software needed to analyze x-ray spectra, CIAO from the Chandra X-Ray Observatory, and HEASoft from NASA have been combined into a Docker container. This removes the difficulty of installing and configuring the software, facilitates repeatability, stability of results, allows the full analysis pipline to be agnostic of the host machine, and allows it to be run on a compute cluster if needed.

## Quickstart

When you clone this repo there are 3 main sections inside spec_analysis alongside some scripts/files that Docker needs:

-Dockerfile
-entrypoint.sh
-Makefile
-requirements.txt
-README.md

-/data
-/docs
-/tools

To run the entire process, assuming you have Docker working on your machine: In the same directory as the Docker file, run this command in a terminal, "make run". This will build the docker image, run container from the image, and finally, run every command in the entrypoint.sh file.

## Dockerfile
This is the Dockerfile for the containter, this is what is used to build the image that you will use to create/run a container. This file dictates the behavior of the container, to make changes to the container you change the Dockerfile to build a new image with those changes, then create/run a container from that point. You do not save containers, they are created/run in the state you intend them to be, and the code you run can be repeated reliably. For more info look at the docker links in the Notes section.

## entrypoint.sh
This script is used to contain all the scripts/commands to be run in succession in the container. When the container runs it will begin by running all the code in this file. These are all the steps to be automated in one place for this project. To make changes to what scripts run make changes to this file.

## Makefile
This is a Makefile for the convience of creating the docker images, creating the container/running them, automating and parameterizing the docker commands and process for ease of use. When wanting to create an image/container, make changes to docker commands, or anything with docker, it is all managed in this file, and you can just run the make commands to run docker easily.

## requirements.txt

This file is where all the python packages to be installed in the container are managed. To make changes to what python packages you want the container to have, edit this file.

## README.md

The file you are reading right now, where general documentation is managed. To make changes to the documentaion edit this file.

## Data

This is the directory where you can put as many observation files as you want to analyse. This is where all the input data will be referenced for the use of all the scripts. This directory is mounted on locally on your machine for the container to use. So no data is ever kept in the Docker contianer, just pulled from your machine on an identically named directory. 
When finished, any output data will automatically be placed in an "out" subdirectory on your machine(do not delete this directory). This way you can run the container as many times as you want and can keep the input and output data in one place seperate from the process.

## Docs

Relevant papers, instructions on the Data Prep and Spectral Extraction steps, and some reference manuals for the analysis tools used.

## Tools

This is where all the scripts used to complete the Data Prep and Spectral Extraction steps are kept.

## Notes

https://docs.docker.com/reference/dockerfile/

https://docs.python-guide.org/writing/structure/

https://cda.harvard.edu/chaser/mainEntry.do

https://heasarc.gsfc.nasa.gov/lheasoft/ubuntu.html

https://cxc.cfa.harvard.edu/ciao/download/custom.html
# open_clusters
