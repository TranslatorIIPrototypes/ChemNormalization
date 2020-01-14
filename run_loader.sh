#!/bin/sh
read -p "Have you entered your configuration parameters into config.json [y/n]? " answer

if [ $answer = 'n' ]
then
  echo "Please input your settings into config.json and try again. Exiting..."
  exit 1
fi

echo "Activating virtual python environment..."

# Activate the virtual environment (linux)
source activate chemNodeVenv

# Once the configuration is completed, initiate the load.

read -p "Are you ready to start the loading [y/n]? " answer

if [ $answer = 'n' ]
then
  echo "Loading aborted. Exiting..."
  exit 1
fi

python src/load-redis.py

