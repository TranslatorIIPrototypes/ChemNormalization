# ChemNormalization
Service that produces Translator compliant nodes for a .

## Installation
##### Clone the repo
    
    > cd <to directory for the cloned repo>
    > git clone https://github.com/TranslatorIIPrototypes/ChemNormalization.git
     
##### Install and configure Anaconda (if is is not already installed)

    > cd <to a directory for the install>
    > wget http://repo.continuum.io/archive/Anaconda3-4.0.0-Linux-x86_64.sh
    > bash Anaconda3-4.0.0-Linux-x86_64.sh
    > conda update

##### Add the conda-forge package install channel

    > conda config --add channels conda-forge

##### Create a conda virtual environment

    > conda create -n chemNodeVenv python=3.7

##### Activate the virtual environment (linux)

    > source activate chemNodeVenv

##### Install package requirements 

    > conda install --yes --file requirements.txt
         
##### Informational conda commands
    To remove a conda virtual environment
    
        > conda env remove -n chemNodeVenv
     
     To list all conda virtual environments
     
        > conda env list

## Loading Redis

#### Starting the redis server 
We recommend using 
[R3 (Redis-REST with referencing)](https://github.com/TranslatorIIPrototypes/r3) to create, install and start a Docker container serving a Redis instance. 

##### Configure the loader runtime settings

    Configurable run settings in ./config.json
    
    {
        "redis_port": "redis-port>,
        "redis_host": "<redis-host>",
        "redis_password": "<redis-password>",
    
        "neo4j_uri": "<neo4j-host>,
        "neo4j_user": "<neo4j-user>",
        "neo4j_password": "<neo4j-password",
    
        "debug_record_limit": "<limit=#>",
        "debug_messages": <0 or 1>
    }   

##### Initiate the load

Once the configuration is completed, run the load file. 
 
    > python src/load-redis.py
    

