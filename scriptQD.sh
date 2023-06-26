#!/bin/bash

PYTHON=/usr/bin/python
PYTHON=python

#run executable on, in order, hourly, 5min, 1min
cd /n/space_data/MagModelInputs/QinDenton
year=$(date -d "$(date +%Y-%m-15) -2 month" +"%Y") #run for whatever year it was two months ago...
echo $year
DIRECTORY=/n/space_data/MagModelInputs/QinDenton/$year
if [ ! -d "$DIRECTORY" ]; then
  mkdir $DIRECTORY
fi
echo "60 $year 1" | /n/space_data/MagModelInputs/QinDenton/runQD
echo "5 $year 1"  | /n/space_data/MagModelInputs/QinDenton/runQD
echo "1 $year 1"  | /n/space_data/MagModelInputs/QinDenton/runQD

#now add headers and clean dates that are known to be fill

# prepend header to all the source files 
for f in */QinDenton*hour; do cat QDheader.txt $f > $f.txt; rm $f; done
for f in */QinDenton*5min; do cat QDheader.txt $f > $f.txt; rm $f; done
for f in */QinDenton*1min; do cat QDheader.txt $f > $f.txt; rm $f; done

${PYTHON} /n/space_data/MagModelInputs/QinDenton/removeQDAfterDate.py -a

#copy files to folder needed for MagEphem

# DIRECTORY_ME=/n/space_data/QinDenton/$year
# if [ ! -d "$DIRECTORY_ME" ]; then
#   mkdir $DIRECTORY_ME
# fi
# cp $DIRECTORY/*.txt $DIRECTORY_ME
