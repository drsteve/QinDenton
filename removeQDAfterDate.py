#!/usr/bin/python2.6

import os, glob, urllib.request, urllib.parse, urllib.error, urllib.request, urllib.error, urllib.parse, re
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

    parser.add_option("-t", "--Test", action="store_true", dest="Test",
                            help="Don't delete any files. Dry run only.")

    return parser


def removeFilesAfterDate(indate, test):
    #find directories of input year to end of current year
    fnames = []
    start_year = indate.year
    this_year = dt.datetime.now().year
    if this_year<start_year: raise Exception
    check_years = [start_year+n for n in range(1+this_year-start_year)]
    for year in check_years:
        filelist = sorted(glob.glob('{0}/QinDenton_*.txt'.format(year)))
        for fn in filelist:
            fdate = dup.parse(re.search('\d{8}', fn).group())
            if fdate > indate:
                if not test:
                    os.unlink(fn)
                else:
                    fnames.append(fn)
    return fnames


def findLastValidDate():
    '''Parse the OMNI webpage to find the last valid date for IMF and plasma data'''
    siteurl = 'https://omniweb.gsfc.nasa.gov/html/omni_min_data.html'
    proxies = urllib.request.ProxyHandler({'https': 'http://proxyout.lanl.gov:8080/'})
    opener = urllib.request.build_opener(proxies)
    urllib.request.install_opener(opener)
    response = urllib.request.urlopen(siteurl)
    dat = response.readlines()
    candidates = [line for line in dat if 'IMF and Plasma' in line]
    if len(candidates)>=1:
        try:
            last_valid = re.search('(\d{4}-\d{2}-\d{2} \(\d{3}\)) - (\d{4}-\d{2}-\d{2} \(\d{3}\))', candidates[0]).group(2)
            last_valid = re.search('(\d{4}-\d{2}-\d{2})', last_valid).group(1)
        except:
            raise Exception("Can't uniquely determine a valid end date.")            
    else:
        raise Exception("Can't uniquely determine a valid end date.")
    #now parse last valid (yyyy DOY)
    vdate = dt.datetime.strptime(last_valid, '%Y-%m-%d')
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

    flist = removeFilesAfterDate(indate, options.Test)
    if options.Test:
        print('Last valid data determined to be {0}'.format(indate))
        print('Files to be removed:')
        for ff in flist:
            print(ff)
