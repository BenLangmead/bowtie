#!/bin/sh

if [[ -z $SFUSER ]]
then
    read -p "Username: " SFUSER
fi

rsync -avz --exclude='.*' --exclude='*.sh' "$PWD/" "$SFUSER@web.sourceforge.net:/home/project-web/bowtie-bio/htdocs"
