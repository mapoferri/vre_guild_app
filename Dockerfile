FROM python:3.7
#Define dependencies and Bash shell
ENV DEBIAN_FRONTEND noninteractive

WORKDIR  /home

RUN apt-get update -qq ; apt-get upgrade ; \
    apt-get install -y gridengine-client git ; \
    apt-get install -y build-essential tar curl ; \
    apt-get autoremove -y && apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

#RUN chsh -s /bin/bash
#SHELL ["/bin/bash", "-c"]
#RUN source venv/bin/activate
RUN pwd

#Install python and other dependencies for the Tool

RUN apt-get update ; apt-get upgrade ; \
    apt-get install -y python3-pip python-dev build-essential; \
    apt-get autoremove -y && apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

#RUN update-alternatives --set python /usr/bin/python3.7
#RUN update-alternatives --install /usr/bin/python python /usr/bin/python3.7 1

RUN python3 --version
RUN pip3 install --upgrade pip

RUN pip3 install biobb-analysis==3.7.0 biobb-common==3.7.0 biobb-io==3.7.0 biobb-model==3.7.0 biopython==1.79 pandas opencv-python networkx biobb-guild==1.2
    

#Define VRE tool download 
WORKDIR /home
WORKDIR /home
RUN git clone --branch main https://github.com/mapoferri/guild_app.git
RUN pwd
RUN cd /home/guild_app
WORKDIR /home/guild_app
RUN ls


COPY disc4all-data/ /home/disc4all-data/
COPY datasets/ /home/datasets/
#RUN tar -xzvf datasets.tar.gz
RUN ls

#RUN conda env create -f environment.yml 
#RUN conda activate biobb_guild
#COPY disc4all-data/ /home/vre_app_3dshaper/disc4all-data/
#COPY datasets/ /home/vre_app_3dshaper/datasets/


RUN pip3 install --upgrade wheel
RUN pip3 install -r requirements.txt
RUN python /home/guild_app/setup.py build
RUN python /home/guild_app/setup.py install


#RUN python /home/vre_app_guild/biobb_guild build
RUN pip3 install --upgrade wheel
RUN pip3 install -r requirements.txt
#RUN python /home/vre_app_guild/biobb_guild install

RUN ls

CMD ./VRE_RUNNER --config tests/basic/config.json --in_metadata tests/basic/in_metadata.json --out_metadata out_metadata.json --log_file VRE_RUNNER.log
