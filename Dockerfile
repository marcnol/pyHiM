# set base image (host OS)
FROM python:3.7

# updates libraries so that opencv2 does not complain
RUN apt-get update 
RUN apt-get install ffmpeg libsm6 libxext6  -y

# set the working directory in the container
WORKDIR /code

# copy the dependencies file to the working directory
COPY requirements.txt .

# install dependencies

RUN pip install -r requirements.txt
RUN pip install --upgrade scikit-image
RUN pip install --upgrade tensorflow
RUN pip install dask distributed

# copy the content of the local src directory to the working directory
COPY src/ .
#RUN cp ./*py .
#RUN cp -rf matrixOperations imageProcessing fileProcessing .

# command to run on container start
ENTRYPOINT [ "python", "./pyHiM.py"]

