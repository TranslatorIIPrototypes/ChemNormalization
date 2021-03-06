# ChemNormalization

[![Build Status](https://travis-ci.com/TranslatorIIPrototypes/ChemNormalization.svg?branch=master)](https://travis-ci.com/TranslatorIIPrototypes/ChemNormalization)

A ROBOKOP Translator compliant lookup service that provides a list of similar chemical substances given a standardized chemical identifier.

This project uses RDKit (https://rdkit.org/) to calculate a simplified SMILES value for each chemical substance in a ROBOKOP database. The simplified SMILES are then used as a mechanism to group similar chemical substances whose relationships are loaded into a Redis database. Access to this lookup information is provided by a Swagger web UI/service which accesses the Redis database directly. 

This project relies on a redis/swagger project called r3. To learn more on this product please refer to its' Github repo located at: https://github.com/TranslatorIIPrototypes/r3.git

Note: This project uses an Anaconda virtual environment because RDKit (used to calculate simplified SMILES) which is currently unavailable in pip package repositories.
## Installation

Note: This environment expects Python version 3.7 and, due to RDKit, uses conda package repositories. 

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

##### Change your working directory

    > cd <ChemNormalization cloned repo location>

##### Activate the virtual environment (linux)

    > source activate chemNodeVenv
    
##### Configure the loader runtime settings

    Configurable run settings in ./config.json
    
    {
        "redis_port": "redis-port>",
        "redis_host": "<redis-host>",
        "redis_password": "<redis-password>",
    
        "neo4j_uri": "<neo4j-host>",
        "neo4j_user": "<neo4j-user>",
        "neo4j_password": "<neo4j-password",
    
        "debug_record_limit": "<leave empty or limit #>",
        "debug_messages": <0 (no msgs shown) or 1>
                
        "do_KGX": <0 (do not create KGX files) or 1>,
        "do_Redis": <0 (do not load Redis) or 1>
            
        "output_node_file": "<KGX node file path>.csv",
        "output_edge_file": "<KGX edge file path>.csv",
    
        "node_norm_endpoint": "https://nodenormalization-sri.renci.org/get_normalized_nodes",
        "node_norm_chunk_size": <record batch size for node normalizaton>,
    
        "do_curie_update": <0 (no curie update processing) or 1>,
        "do_node_normalization": <0 (no node normalization processing) or 1>,
    }   

##### Once the configuration is completed, initiate the load. 
 
    > python src/load.py

## Creating a Neo4j database for PLATER (loaded with KGX files)

    This project has docker configuration files to create and load KGX data into a Neo4j database. 
    
    The Neo4j database can then be consumed by PLATER.
    
### Docker setup and commands
    Change your working directory to "ChemNormalization/graph" and edit the .env file. 
    
    Configure with appropriate settings for the local implementation of the Neo4j container.

        NEO4J_HOST=0.0.0.0
        NEO4J_HTTP_PORT=7474
        NEO4J_BOLT_PORT=7687
        NEO4J_HTTPS_PORT=7473
        NEO4J_PASSWORD=<password>

    The docker configuration creates 3 directories on the host for data, authentication and ssl certificates. 
    
    Please insure these directories have been created prior to running the container. these paths are relative to your working directory.
    
        ../../../neo4j_data:/data
        ../../../neo4j_logs:/logs
        ../../../neo4j_ssl:/ssl
    
    Once configured, launch the following docker commands:
    
        > docker-compose build
        > docker-compose up -d 
