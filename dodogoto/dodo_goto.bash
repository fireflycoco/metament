#!/bin/bash -f

dest=`dodo_goto_parser ${*}`
if [ "${dest}" != "." ] 
then
	cd "$dest"
fi
