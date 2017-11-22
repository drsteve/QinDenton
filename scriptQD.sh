#!/bin/bash

export PYTHONPATH=/packages/lib/python-new/site-python:/packages/lib/site-python
PYTHON=/usr/bin/python2.6

#first make sure the KpDst.lst file doesn't have its weird auto-header
${PYTHON} /n/space_data/OMNI/cleanKpDstFile.py

#then update the Kp from the World Data Centre in Potsdam
${PYTHON} /n/space_data/OMNI/updateKpDst.py

#run executable on, in order, hourly, 5min, 1min
cd /n/space_data/MagModelInputs/QinDenton
year=$(date -d "$(date +%Y-%m-15) -2 month" +"%Y") #run for whatever year it was two months ago...

DIRECTORY=/n/space_data/MagModelInputs/QinDenton/$year
if [ ! -d "$DIRECTORY" ]; then
  mkdir $DIRECTORY
fi
echo "60 $year 1" | /n/space_data/MagModelInputs/QinDenton/runQD
echo "5 $year 1"  | /n/space_data/MagModelInputs/QinDenton/runQD
echo "1 $year 1"  | /n/space_data/MagModelInputs/QinDenton/runQD

#now add headers and clean dates that are known to be fill

source add_header
${PYTHON} /n/space_data/MagModelInputs/QinDenton/removeQDAfterDate.py -a
