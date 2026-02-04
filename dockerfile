#Use base Python image

FROM python:latest

#Copy python code into container(code into working dir and data mounted to volume on PC)

RUN mkdir /root/open_clusters/
WORKDIR /root/open_clusters/
COPY requirements.txt .
COPY pm_membership.py .

#Install any requirements

RUN python -m pip install --upgrade pip
RUN python -m pip install -r requirements.txt

#Run python program

CMD [ "python", "./pm_membership.py" ]


##For every change to this file you must 
#docker rm --force <ContainerName>
#Make (rebuilds image to reflect changes)

#to run lastest docker img with mounted volume: docker run -v ~/open_clusters/data:/root/open_clusters/data -it $(docker build -q .)

