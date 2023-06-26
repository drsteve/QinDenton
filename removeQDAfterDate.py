#!/usr/bin/python2.6

import os, glob, re
import urllib.request
import datetime as dt
from optparse import OptionParser
import dateutil.parser as dup
import spacepy.time as spt


def parserSetup():
    # Define a command-line option parser and add the options we need
    parser = OptionParser(  usage="%prog [options]",\
                            version="%prog Version 1.0 (August 6, 2015)"  )

    parser.add_option("-d", "--Date",      dest="Date",
                            help="Last valid date in ISO format",
                            metavar="YYYY-MM-DD")

    parser.add_option("-a", "--Auto", action="store_true", dest="Auto",
                            help="Automatically determine last valid date")

    return parser


def removeFilesAfterDate(indate):
    #find directories of input year to end of current year
    start_year = indate.year
    this_year = dt.datetime.now().year
    if this_year<start_year: raise Exception
    check_years = [start_year+n for n in range(1+this_year-start_year)]
    for year in check_years:
        filelist = sorted(glob.glob('{0}/QinDenton_*.txt'.format(year)))
        for fn in filelist:
            fdate = dup.parse(re.search('\d{8}', fn).group())
            if fdate > indate:
                os.unlink(fn)


def findLastValidDate():
    '''Parse the OMNI webpage to find the last valid date for IMF and plasma data'''
    proxies = {'http': 'http://proxyout.lanl.gov:8080/'}
    #opener = urllib.FancyURLopener(proxies)
    opener = urllib.request.FancyURLopener()
    siteurl = 'http://omniweb.gsfc.nasa.gov/html/omni_min_data.html'
    out = opener.open(siteurl)
    dat = out.readlines()
    candidates = [line.decode('latin1') for line in dat if 'IMF and Plasma' in line.decode('latin1')]
    if len(candidates)>=1:
        try:
            last_valid = re.findall('\d{4}-\d{2}-\d{2}', candidates[0])
            last_valid[1]
        except:
            raise Exception("Can't uniquely determine a valid end date.  (try)")            
    else:
        raise Exception("Can't uniquely determine a valid end date.")
    #now parse last valid (yyyy DOY)
    #yy, doy = last_valid.split()
    #mm, dd = spt.doy2date(int(yy), int(doy))
    vdate = dt.datetime.strptime(last_valid[1],'%Y-%m-%d')
    #print vdate
    return vdate


if __name__=='__main__':
    parser = parserSetup()
    # Parse the args that were (potentially) given to us on the command line
    (options, in_args) = parser.parse_args()
    
    #roll defaults into options dict and do any checking
    valid_opt = 0
    for key in options.__dict__.keys():
        if options.__dict__[key] != None:
            valid_opt += 1

    #check for any input arguments, if not, print help
    if not valid_opt:
        parser.print_help()
        exit()

    if not options.Auto:
        indate = dup.parse(options.Date)
    else:
        indate = findLastValidDate()
    removeFilesAfterDate(indate)
