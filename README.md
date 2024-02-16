# ProjectWeb

**This project was developed as part of the Web Programming course of the M2 AMI2B of Paris-Saclay University.**

The app is called **BactaHub**, it was created to facilitate the exchange of information in the biological domain between three different figures who approach the world of gene annotations: users, annotators and validators. 

The application allows people with biological backgrounds to annotate bacterial genomes and specific genes or proteins.
This web-application offers the possibility to make queries on already annotated genomes and genes, to download the results in a fasta format and to upload new sequences, already noted or needing to be cataloged, in the BactAHub database. 
To increase and make more transparent communication between users, a forum is made available to the BactAHub community.

## Project guidelines

To access our project, it has to be cloned on your computer with the following command :

`git clone <url>` 
"url" has to be replaced by the url found on the first page of this git repository

To developp this project, we used a **python3 virtual environment** by executing the followwing command:

`python3 -m venv .venv`

`source .venv/bin/activate`

Then, to create the correct work environment for this project, the needed **installations** are in the file env_install.sh :  
`source env_install.sh`

## Django command

Populate the database with fasta files and users :

`python3 manage.py import-my-data file/path/`

Options :

`--annotated`

Takes into account the fully annotated fasta files

`--populatewithusers`

Includes in the database some test users for the site

Once your database is filled, you can run the following command to **launch the site** :

`python3 manage.py runserver`

It will give you an **http address** to paste in your browser to access the site

If you want to use the site has an **admin**, you have to execute the following command:

`python3 manage.py createsuperuser`



