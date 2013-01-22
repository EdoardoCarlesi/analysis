#!/bin/bash

# Commit changes locally:
git add *.file
git commit 

git clone git@comodo.ft.uam.es:/home/git/PROJECTS/gadget.git DIR

git remote add REMOTENAME url

# Remotename = origin 
# url = git@comodo.uam.es:/home/edoardo/git/PROJECTS/project.git
 
#make the files editable
cd url 
chmod -R 777 *

# Branch = master 

# commit changes
git push REMOTE BRANCH

# get changes on another machine
git pull REMOTE BRANCH
