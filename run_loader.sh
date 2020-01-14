#!/bin/sh
read -p "Have you entered your configuration parameters into config.json [y/n]?" answer

if [ $answer = 'n' ]
then
  echo "Please input your settings into config.json and try again. Exiting..."
  exit 1
fi

# Activate the virtual environment (linux)
bash source activate chemNodeVenv

# Once the configuration is completed, initiate the load.
#python src/load-redis.py

