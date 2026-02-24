# open_clusters
Automate the the identification of the blue straggler population, subsubgiant population, variables, etc in open clusters. We identify variability and the spectral energy distributions for sources that are members of these clusters.

Search a catalog of old open clusters (> 1 Gyr) that fall within the catalog of the Zwicky transient factory (ZTF). Because ZTF has multiple filters and epochs of photometry for the stars associated with these clusters. 

Then retrieve all the objects that have GAIA information within the area of the cluster.
Then define members to the clusters based on the GAIA proper motions, and define the membership probability for the stars associated with the clusters. 

Then match those objects back with the photometry from ZTF, and other surveys that could add information.



# Quickstart

When you clone this repo there are 3 main sections inside open_clusters alongside some scripts/files that Docker needs:

Dockerfile,
Makefile,
requirements.txt,
pm_membership.py,
README.md,
/data,
/data/out


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

## pm_membership.py

1) Load and isolate motion data:
We read our catalog and extract only:
pmrapmdecIf either is missing, we drop the star because we cannot judge its kinematic allegiance without both velocity components

2) Assume two overlapping populations:
We assume every star belongs to one of two kinematic species:
A tight cluster distributionA broad Galactic field distributionWe model both as 2D Gaussians in proper-motion space

3) Define and freeze the field:
We estimate the global mean and covariance of the full sample, inflate the dispersion, and freeze this as our field model. We lock it to prevent the field from reshaping itself around the cluster during fitting.

4) Seed the cluster center:
We find the global median motion, select the densest ~15% of stars near it, and use their median as our initial cluster center. This gives us a robust starting guess for the kinematic peak.

5) Estimate initial cluster dispersion:
From stars nearest that center, we estimate a characteristic velocity spread. This defines our first cluster Gaussian

6) Compute membership probabilities:
For each star, we compute:
Cluster likelihood and Field likelihood
We combine them using Bayes’ theorem with a fixed prior (π = 0.15) to obtain: P(cluster∣motion)
Each star now has a soft membership weight between 0 and 1

7) Refit the cluster (weighted update):
We use the probabilities as weights to recompute:
The cluster centerThe covariance matrixHigh-probability stars pull harder. Low-probability stars barely tug

8) Iterate to convergence:
We repeat probability calculation and refitting until the cluster center stops shifting beyond a tolerance. At that point, the model and its members are self-consistent

9) Measure Mahalanobis distance:
Using the final covariance, we compute each star’s elliptical distance from the cluster center. This measures separation in units of the cluster’s natural velocity shape

10) Assign membership tiers:
We discretize the continuous probabilities:
≥ 0.9 → Core0.5–0.9 → Probable0.2–0.5 → Candidate< 0.2 → FieldThese are confidence labels layered on top of the Bayesian posterior

11) Validate with CMD:
We plot a color–magnitude diagram using our tiers
If our kinematic model is correct, core members trace a clean main sequence while field stars scatter chaotically

# Notes

https://docs.docker.com/reference/dockerfile/

https://docs.python-guide.org/writing/structure/
