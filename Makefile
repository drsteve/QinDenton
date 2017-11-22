
# Location to place downloaded OMNI data from NSSDC ftp site
DIRTMP=./tmp

all:
	make omniupdate
	make ViRBO
	make header

ViRBO:
	make build
	echo "60 2015 1" | ./runQD
	bash add_header

build:
	g77 -o runQD Code.f

omniupdate:
	rm -f hour/omni2_hour.dat
	rm -f 5min/omni_5min.asc
	rm -f 1min/omni_min.asc
	rm -f 5min/kpdst.lst
	make hour/omni2_hour.dat
	make 5min/omni_5min.asc
	make 1min/omni_min.asc
	make 5min/kpdst.lst

hour/omni2_hour.dat:
	mkdir -p $(DIRTMP)
	cd $(DIRTMP); wget -m -l1 ftp://spdf.gsfc.nasa.gov/pub/data/omni/low_res_omni/omni2*
	cd $(DIRTMP)/spdf.gsfc.nasa.gov/pub/data/omni/low_res_omni/; mv omni2* /n/space_data/OMNI/hour/
	cd hour; ln -s /n/space_data/OMNI/hour/omni2_all_years.dat omni2_hour.dat 

5min/kpdst.lst: hour/omni2_hour.dat
	cut -c 1-11,219-221,226-231 hour/omni2_hour.dat > tmp.dat
	more +280513 tmp.dat > 5min/kpdst.lst
	rm -f tmp.dat
	cd 1min ; ln -s -f ../5min/kpdst.lst 

5min/omni_5min.asc:
	mkdir -p $(DIRTMP)
	cd $(DIRTMP); wget -m -l1 ftp://spdf.gsfc.nasa.gov/pub/data/omni/high_res_omni/omni_5min*
	find $(DIRTMP)/spdf.gsfc.nasa.gov/pub/data/omni/high_res_omni/ -name "omni_5min*.asc" | sort | xargs -i cat {} > /n/space_data/OMNI/5min/omni_5min.asc	
	cd $(DIRTMP)/spdf.gsfc.nasa.gov/pub/data/omni/high_res_omni/; mv omni_5min* /n/space_data/OMNI/5min/
	cd 5min; ln -s /n/space_data/OMNI/5min/omni_5min.asc

1min/omni_min.asc:
	mkdir -p $(DIRTMP)
	cd $(DIRTMP); wget -m -l1 ftp://spdf.gsfc.nasa.gov/pub/data/omni/high_res_omni/omni_min*
	find $(DIRTMP)/spdf.gsfc.nasa.gov/pub/data/omni/high_res_omni/ -name "omni_min*.asc" | sort | xargs -i cat {} > /n/space_data/OMNI/1min/omni_min.asc	
	cd $(DIRTMP)/spdf.gsfc.nasa.gov/pub/data/omni/high_res_omni/; mv omni_min* /n/space_data/OMNI/1min/
	cd 1min; ln -s /n/space_data/OMNI/1min/omni_min.asc

