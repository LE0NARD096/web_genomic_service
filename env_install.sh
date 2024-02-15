#!/bin/bash
python3 -m venv .venv
source .venv/bin/activate
pip install django
pip install biopython
pip install loader
pip install plotly
pip install django-phonenumber-field
pip install django-phonenumbers