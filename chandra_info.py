#! /usr/bin/env python

#Author: Kari Frank
#Date: June 11, 2013
#Purpose: Print out basic information about a Chandra
#          observation and event file
#
#Usage: python chandra_info.py [EVT_IN] [--region REGION]

#Input:
#
# evt_in: provide the evt2 file from which to read 
#          information. if not provided will choose 
#          the first file (in pwd) it finds with 
#          'evt' in the file name.
#
# region:  optionally provide a region file to return only
#          counts and chips in that region
#
#Output:
# - prints to information to screen
#
#Usage Notes:
# - CIAO must be initialized prior to running this script (heainit, setciao)
# - This script uses the CIAO installation of python
# - input evt2 file must be in the working directory unless full path is 
#   given in evt_in

#----Import Modules----
import argparse,os, sys
import ciao_contrib.runtool as crt #all functions should be prefixed with crt
import home_grown as hg #import my own custom functions
import chandra #import my own chandra related functions
from math import *
from numpy import *

#----Parse arguments and set defaults---

parser = argparse.ArgumentParser(description='Print basic information about Chandra evt2 file.',epilog='NOTE: Must initialize CIAO before running this script ("heainit" then "setciao").')

pwd = os.getcwd()+'/'

parser.add_argument('evt_in',help='Input event file',default='default')

parser.add_argument('--region',help='Optional region file for determining counts and chip coverage.',default='none')

args = parser.parse_args()

#--rename args--

if args.evt_in == 'default':
  evt_infile = hg.fetch_file(pwd,pat='evt',prompt='yes')
else:
  evt_infile = args.evt_in
print "\nInput File: "+evt_infile#+'\n'

#----Observation Information----
#obsid, target name, pointing, observation mode, read mode,
# data mode, detector names

#--get obsid, object name, dates, and modes--
obsid = chandra.get_obsid(args.evt_in)
objectname = chandra.get_objectname(args.evt_in)

start_date = crt.dmkeypar(args.evt_in,"DATE-OBS","echo+")

obs_pointing = chandra.get_pointing(args.evt_in)
obs_mode = crt.dmkeypar(args.evt_in,"OBS_MODE","echo+")

data_mode = crt.dmkeypar(args.evt_in,"DATAMODE","echo+")

det_name = crt.dmkeypar(args.evt_in,"DETNAM","echo+")
common_name = ''
if det_name == 'ACIS-0123':
  common_name = '(ACIS-I)'
if det_name == 'ACIS-456789':
  common_name = '(ACIS-S)'

if det_name[:4]=='ACIS': 
  read_mode = crt.dmkeypar(args.evt_in,"READMODE","echo+")
#focal plane temperature
  fp_temp_name = crt.dmkeypar(args.evt_in,"FP_TEMP","echo+")
  fp_temp = str(float(fp_temp_name)-272.15) #convert Kelvin to Celsius
  frametime = crt.dmkeypar(args.evt_in,"EXPTIME","echo+")
  
grating_name = crt.dmkeypar(args.evt_in,"GRATING","echo+")


#--if sn1987a, calculate day number from the supernova--

age='0'
if  objectname == 'SNR1987A' or objectname == 'SN1987A' or objectname == 'SN 1987A' :

  shortdate = start_date.split('T')
#  sys.path.insert(0,'/astro/research/kaf33/Dropbox/Research/Python_Programs/sn1987a')
#  from sn1987a_time import convert_time
#  age = str(convert_time(shortdate[0]))

####old way####
#elapsed_obs_date = 0.0
  
  #get julian date of observation
#  julian_obs_date = float(crt.dmkeypar(args.evt_in,"MJD_OBS","echo+"))

  #get julian date of explosion (day 0)
  #Note this date is only accurate to within ~24 hours!
  #calculated assuming start of obsid 12539 observation was day 8796
  # (Helder2012)
  #from obsid12539 evt file header MJD_OBS value
#  julian_ref_date = float(5.5645617042927)*10.0**4.0
#  elapsed_ref_date = 8796
#    sn_time = Time('1987-02-23','iso',scale='utc')

  #calculate day since day 0
#  elapsed_obs_date = int(julian_obs_date - julian_ref_date + elapsed_ref_date)

print "-----------------------------------------------"
print "OBSERVATION INFO"
print "Object Name   : "+objectname
print "Obs ID        : "+obsid
print "Start Date    : "+start_date
#if elapsed_obs_date != 0.0: print "SN1987A Day   : "+str(elapsed_obs_date) 
if age != '0': print "SN1987A Age   : "+age 
print "Obs RA,DEC    : "+obs_pointing[0]+', '+obs_pointing[1] 
print "Obs Roll      : "+obs_pointing[2]
print "Obs Mode      : "+obs_mode
if det_name[:4]=='ACIS': 
  print "Read Mode     : "+read_mode
print "Data Mode     : "+data_mode
print "Detector      : "+det_name+' '+common_name
print "Grating       : "+grating_name
if det_name[:4]=='ACIS': 
  print "FP Temperature: "+fp_temp
if det_name[:4]=='ACIS':
  print "Frametime     : "+frametime
print "-----------------------------------------------"
#print ""


#----Event File Information----
#number of counts, chips, exposure time

#--Get number of counts--
if args.region!='none':
  evt_reg_str = args.evt_in+'[sky=region('+args.region+')]'
else:
  evt_reg_str = args.evt_in

num_counts = crt.dmlist(evt_reg_str,"COUNTS")
if det_name[:4]=='ACIS': chip_ids = crt.dmkeypar(evt_reg_str,"CCD_ID","echo+")
#crt.dmstat.punlearn()
#chip_ids = crt.dmstat(evt_reg_str+'[cols ccd_id]')

#--Get exposure time--
exp = crt.dmkeypar(args.evt_in,"EXPOSURE","echo+")

count_rate = str(float(num_counts)/float(exp))


print "-----------------------------------------------"
print "EVENT FILE INFO"
print "Exp. Time  : "+exp+" s"
print "Counts     : "+num_counts
print "Count rate : "+count_rate+" counts/s"
if det_name[:4]=='ACIS': 
  print "CCD_ID     : "+chip_ids
print "-----------------------------------------------"
print ""
