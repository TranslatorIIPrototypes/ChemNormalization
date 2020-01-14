#!/bin/sh
read -p "Have you entered your configuration parameters into config.json [y/n]? " answer

if [ $answer = 'n' ]
then
  echo "Please input your settings into config.json and try again. Exiting..."
else
  echo "Activating virtual python environment..."

  # Activate the virtual environment (linux)
  source activate chemNodeVenv

  read -p "Are you ready to start the loading [y/n]? " answer

  if [ $answer = 'n' ]
  then
    echo "Loading aborted. Exiting..."
  else
    echo "Loading..."
    python src/load-redis.py
  fi
fi


