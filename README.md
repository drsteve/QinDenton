# Qin-Denton Data Files for Magnetic Field Modeling
This code is modified from Richard Denton's original code
(described in [Qin et al., 2007](https://doi.org/10.1029/2006SW000296)),
as the original required all input data
to be read at-once and processed from the start date to
present in a single block.

This modified version:
- Processes the data in yearly blocks
  - This will try to read the initial coefficients from the previously processed year. If those do not exist, default parameters will be used and it may take a few weeks to obtain exact results, while the time-weighted parameters evolve to the correct solution.
- Writes the output to daily files with a JSON header, as used by [LANLGeoMag](https://github.com/drsteve/LANLGeoMag) and [SpacePy](https://github.com/spacepy/spacepy)

Using the code has several steps as described below. The
brief description is that the Qin-Denton files perform
post-processing on the NASA OMNI data set. The script `scriptQD.sh`
is intended to run the Qin-Denton code at hourly, 5-minute,
and 1-minute resolution, assuming an updated OMNI database.

## Download the OMNI data
Two options:
1. The method employed at LANL is to use the update script `updateOMNI.sh` on the `space_data` volume. This script is run automatically once per week.
2. Use the makefile provided (with any appropriate path updates) and run `make omniupdate`

The 1963 - present hourly data are downloaded from 
`spdf.gsfc.nasa.gov/pub/data/omni/low_res_omni/omni2_all_years.dat`
 with output in file hour/omni2_hour.dat
 then the 1995-todate KpDst data with 1h cadence are extracted from that download

The 1981-todate 1 min and 5 min data are downloaded from
`spdf.gsfc.nasa.gov/pub/data/omni/high_res_omni/`
 with output in `1min/omni_1min.asc and 5min/omni_5min.asc`
 
## Compilation
The code requires a Fortran compiler, and the code itself is F77.
The makefile contains a `build` command.
`make build` will compile the code only (no data download, no running code)
Alternatively, use either `g77 -o runQD Code.f` or `gfortran -std=legacy -o runQD Code.f`

## Calculate the "Qin-Denton" parameters
The runQD executable will calculate/make the Qin-Denton data files, but it
requires three inputs from the standard input device. To run this on the
command line without interactive input,

```echo [cadence] [starting year] [no of years] | ./runQD```

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

### To add header to mag field pars files 
type `bash add_header`
Output files are in yearly directories as above and with suffix .txt
 
### To remove daily files for which there are no OMNI data
`./removeQDAfterDate.py -a`

### To calculate mag field parameters AND add header
1. The hard way:
 edit in file Makefile the cadence, starting year (not before 1995), and number of years (not after 2015)
 on line `echo [cadence] [starting year] [no of years] | ./runQD`
 then type `make ViRBO`
 Output files are in yearly directories: year/QinDenton_[year][mo(nth)][da(y)]_[cade(nce)].txt
2. The easy way:
 Run the script `scriptQD.sh`, if necessary modifying the output pathing, start year, etc.
 It will run the runQD executable for all cadences of input data and process the output.
