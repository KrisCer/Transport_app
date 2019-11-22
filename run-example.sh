#!/bin/bash

echo "Installing tools-barebone requirements.."
pip3 install -r requirements.txt

echo "Installing user requirements.."
pip3 install -r user_requirements.txt

echo "Creating logs/requests.log file which if missing doesnt let the app start"
mkdir webservice/logs
touch requests.log

echo "Running v0.1-transport-app application.."
python3 webservice/run_app.py
