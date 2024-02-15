# ProjectWeb

## Project guidelines

To create the correct work environment for this project, the needed installations are in the file env_install.sh :  
`source env_install.sh`

## Django command

Populate the database with fasta files and users

**python3 manage.py import-my-data file/path/**

Options

**--annotated**

Takes into account the fully annotated fasta files

**--populatewithusers**

Includes in the database some test users for the site



## Git command 

Cloner un repository git en local : **git clone url**

Bring modifications on local from the git repo:
**git pull**

Bring modifications from your computer to the git repo : 

**git add .** to push everything 
**git add nameFile** to push specific file(s)

**git commit -m « … »**

**git push**

Vizualise on what branch we are currently on :
**git status**

Forcer à pull même quand erreur :
**git reset --hard origin/main**

Change branch : 
**git checkout NomBranche**

Update a branch with main : 
**git merge main**

