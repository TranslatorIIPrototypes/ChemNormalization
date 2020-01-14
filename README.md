# NodeNormalization
Service that produces Translator compliant nodes given a curie.

## Installation

Create a virtual environment

    > python -m venv chemNormalization-env

Activate the virtual environment

    # on Linux
    > source chemNodemaization-env/bin/activate
    # on Windows
    > source chemNormalization-env/Scripts/activate 

Install requirements 

    > pip install -r requirements.txt

## Loading Redis


### Starting redis server 
The Load script can be used to put data to a running Redis instance. Inline with this we recommend using 
[R3 (Redis-REST with referencing)](https://github.com/TranslatorIIPrototypes/r3). 
### Config
Once we have a running
redis-server we can modify our config file located at `./config.json` as the following.

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
### Loading

Once the configuration is completed, move to the "src" directory and python run the load file. 
 
    > cd src
    > python load-redis.py
    

