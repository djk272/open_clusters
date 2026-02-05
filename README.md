# open_clusters
Proper-motion membership for open clusters


## Quickstart

When you clone this repo there are 3 main sections inside open_clusters alongside some scripts/files that Docker needs:

-Dockerfile
-Makefile
-requirements.txt
-pm_membership.py
-README.md

-/data
-/data/out


To run the entire process, assuming you have Docker working on your machine: In the same directory as the Docker file, run this command in a terminal, "make run". This will build the docker image, run container from the image, and run the python script in it's entirety.

## Dockerfile
This is the Dockerfile for the containter, this is what is used to build the image that you will use to create/run a container. This file dictates the behavior of the container, to make changes to the container you change the Dockerfile to build a new image with those changes, then create/run a container from that point. You do not save containers, they are created/run in the state you intend them to be, and the code you run can be repeated reliably. For more info look at the docker links in the Notes section.


## Makefile
This is a Makefile for the convience of creating the docker images, creating the container/running them, automating and parameterizing the docker commands and process for ease of use. When wanting to create an image/container, make changes to docker commands, or anything with docker, it is all managed in this file, and you can just run the make commands to run docker easily.

## requirements.txt

This file is where all the python packages to be installed in the container are managed. To make changes to what python packages you want the container to have, edit this file.

## README.md

The file you are reading right now, where general documentation is managed. To make changes to the documentaion edit this file.

## Data

This is the directory where you can put as many observation files as you want to analyse. This is where all the input data will be referenced for the use of all the scripts. This directory is mounted on locally on your machine for the container to use. So no data is ever kept in the Docker contianer, just pulled from your machine on an identically named directory. 
When finished, any output data will automatically be placed in an "out" subdirectory on your machine(do not delete this directory). This way you can run the container as many times as you want and can keep the input and output data in one place seperate from the process.

## Notes

https://docs.docker.com/reference/dockerfile/

https://docs.python-guide.org/writing/structure/
