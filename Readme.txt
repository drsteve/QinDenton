
 To download the OMNI data
 ------------------------
 Two options:
 1. The preferred update method is to use the update script on the ISR Scheme
 (/n/space_data/OMNI/updateOMNI.sh). This script is run automatically once per
 week.
 2. If you are running off of the scheme, use the makefile provided and run
 "make omniupdate"
 
  The 1963-todate hourly data are downloaded from 
  spdf.gsfc.nasa.gov/pub/data/omni/low_res_omni/omni2_all_years.dat  
   with output in file hour/omni2_hour.dat
   then the 1995-todate KpDst data with 1h cadence are extracted from that download
  The 1981-todate 1 min and 5 min data are downloaded from
  spdf.gsfc.nasa.gov/pub/data/omni/high_res_omni/
   with output in 1min/omni_1min.asc and 5min/omni_5min.asc
 
 To compile the code
 -------------------
 The makefile contains a "build" command.
 "make build" will compile the code only (no data download, no running code)
 Alternatively, type "g77 -o runQD Code.f"

 To calculate the input parameters for magnetic field models
 -----------------------------------------------------------
 The runQD executable will calculate/make the Qin-Denton data files, but it
 requires three inputs from the standard input device. To run this on the
 command line without interactive input,
 type "echo [cadence] [starting year] [no of years] | ./runQD"
  setting the OMNI data cadence to be used, 
  a starting year (not before 1995), 
  and a number of years (last year should not be after 2015).
 Output files are in yearly directories: year/QinDenton_[year][mo(nth)][da(y)]_[cade(nce)]

 The ViRBO code calculates averages of Kp, interpolates hourly Kp/Dst data, 
  calculates the 1 day average of all parameters, interpolates across data gaps, 
  and outputs the W, Bz, and G parameters for Tsyganenko B-field models in daily files,
  labeled by date and input cadence.
 The W parameters, calculated recursively, are initialized at the values preceding the 
 starting epoch. 

 To add header to mag field pars files 
 ---------------------------------
 type "bash add_header"
 Output files are in yearly directories as above and with suffix .txt
 
 To remove daily files for which there are no OMNI data
 -----------------------------------------------------
 ./removeQDAfterDate.py -a

 To calculate mag field parameters AND add header
 ------------------------------------------------
 1. The hard way:
  edit in file Makefile the cadence, starting year (not before 1995), and number of years (not after 2015)
  on line "echo [cadence] [starting year] [no of years]" | ./runQD 
  then type "make ViRBO"
  Output files are in yearly directories: year/QinDenton_[year][mo(nth)][da(y)]_[cade(nce)].txt
 2. The easy way (on the ISR Scheme):
  Run the script "scriptQD.sh", if necessary, modifying the start year, etc. It will run the runQD executable for all cadences of input data and process the output.

 To complain about the above
 ---------------------------
 please send polite email to alin@lanl.gov and explain your problem
 Note: Alin hisself don't speaks English that goodly, so why don't you send that email to noreply@nowehere.never ?
