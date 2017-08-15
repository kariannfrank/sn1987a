#----------------------------------------------------------
# Misc. Functions for handling dates and times related to sn1987a
#
#
#
#
# Contains the following functions:
#
# get_obs_time
# convert_time
# convert_units
#
#----------------------------------------------------------


#----------------------------------------------------------
#----------------------------------------------------------
# get_obs_time()
#
# Author: Kari A. Frank
# Date: June 16, 2014
# Purpose: Look up sn1987a age for the given observation.
# Usage: get_obs_time(obsid)
#
# Input:
# obsid -- string or integer of observation id number, may include leading zeroes.
#
# get   -- string specifying the format of the returned value.
#          allowed values: 'age' (default), 'year','month','day','date'
#          all but 'date' will return an integer (date returns as a string)
#
# Output:
# - returns int(day#),int(year),int(month),int(day),or str(date)
#
# Usage Notes:
# - only works for observations which the age has already been calculate
#   for.  to add more observations, add them to the list below.

def get_obs_time(obsid, get='age'):

    #---------------------------------------
    #           Import Modules
    #---------------------------------------

    #---------------------------------------
    #   Define Look-up Lists
    #---------------------------------------

    obsid = int(obsid)

#    observations = [1387,122,1967,1044,2831,2832,3829,3830,4614,4615,5579,5580,6668,6669,7636,7637,9142,9144,10855,10926,11090,13131,12539,12540,13735,14697,14698,15809,14417,10222]

#    days = [4608,4711,5036,5175,5406,5559,5789,5978,6157,6358,6529,6713,6913,7094,7270,7445,7524,7799,8000,8232,8433,8617,8796,8975,9165,9523,9523,9713,9885,9168,8169]

#    info = [(1387,4608),(122,4711),(1967,5036),(1044,5175),(2831,5406),(2832,5559),(3829,5789),(3830,5978),(4614,6157),(4615,6358),(5579,6529),(5580,6713),(6668,6913),(6669,7094),(7636,7270),(7637,7445),(9142,7624),(9144,7799),(10855,8000),(10926,8232),(11090,8433),(13131,8617),(12539,8796),(12540,8975),(13735,9165),(14697,9523),(14698,9713),(15809,9885),(14417,9168),(10222,8169),(10130,7986),(9143,7801)]

    info = [
        (124,
         4608,
         '1999-10-60'),
        (1387,
         4608,
         '1999-10-06'),
        (122,
         4711,
         '2000-01-17'),
        (1967,
         5036,
         '2000-12-07'),
        (1044,
         5175,
         '2001-04-25'),
        (2831,
         5406,
         '2001-12-12'),
        (2832,
         5559,
         '2002-05-15'),
        (3829,
         5789,
         '2002-12-31'),
        (3830,
         5978,
         '2003-07-08'),
        (4614,
         6157,
         '2004-01-02'),
        (4615,
         6358,
         '2004-07-22'),
        (5579,
         6529,
         '2005-01-09'),
        (5580,
         6713,
         '2005-07-11'),
        (6668,
         6913,
         '2006-01-28'),
        (6669,
         7094,
         '2006-01-28'),
        (7636,
         7270,
         '2007-01-19'),
        (7637,
         7445,
         '2007-07-13'),
        (9142,
         7624,
         '2008-01-09'),
        (9144,
         7799,
         '2008-07-01'),
        (10855,
         8000,
         '2009-01-18'),
        (10926,
         8232,
         '2009-09-08'),
        (11090,
         8433,
         '2010-03-28'),
        (13131,
         8617,
         '2010-09-28'),
        (12539,
         8796,
         '2011-03-25'),
        (12540,
         8975,
         '2011-09-21'),
        (13735,
         9165,
         '2012-03-28'),
        (14697,
         9523,
         '2013-03-21'),
        (14698,
         9713,
         '2013-09-28'),
        (15809,
         9885,
         '2014-03-19'),
        (14417,
         9168,
         '2012-04-01'),
        (10222,
         8169,
         '2009-07-06'),
        (10130,
         7986,
         '2009-01-05'),
        (9143,
         7801,
         '2008-07-04'),
        (6178,
         6534,
         '2005-01-13'),
        (6345,
         6718,
         '2005-07-16'),
        (9806,
         7627,
         '2008-01-11'),
        (10852,
         7994,
         '2009-01-12'),
        (10221,
         7995,
         '2009-01-13'),
        (10853,
         7997,
         '2009-01-17'),
        (10854,
         7999,
         '2009-01-17'),
        (12125,
         8423,
         '2010-03-17'),
        (12126,
         8423,
         '2010-03-17'),
        (11091,
         8619,
         '2010-09-29'),
        (14344,
         8977,
         '2011-09-22'),
        (15810,
         10071,
         '2014-09-20'),
        (17415,
         10068,
         '2014-09-17')]

    observations, days, dates = zip(*info)

#---------------------------------------
#  Find Matching Obsid and Return Value
#---------------------------------------

    if get == 'age':
        age = days[observations.index(obsid)]
        return age

    if get != 'age':

        date = dates[observations.index(obsid)]
        parts = date.split('-')
        year = parts[0]
        month = parts[1]
        day = parts[2]

        if get == 'date':
            return date
        if get == 'year':
            return int(year)
        if get == 'month':
            return int(month)
        if get == 'day':
            return int(day)

#----------------------------------------------------------


#----------------------------------------------------------
#----------------------------------------------------------
# convert_units()
# Author: Kari A. Frank
# Date: July 24, 2014
# Purpose: Convert between sn1987a related units (arcsec, km, velocities)
# Usage: convert_sn1987a_time(val,fromunit,tounit)
#
# Input:
#
# val  -- value to be converted
#
# fromunit -- string specifying units of the input val, 'r' (units output from
#         image fitting),'pc' (in plane of sky),'arcsec','km','cm','r/yr' (velocity),
#         'arcsec/yr','km/yr','km/s','arcsec/day'
#
# tounit   -- string specifying the units to convert val into. takes
#         same values as from.
#
# Output:
# - returns the converted val
#
# Usage Notes:
# -
#

def convert_units(val, fromunit, tounit):

    #---------------------------------------
    #           Import Modules
    #---------------------------------------

    import numpy as np

#---------------------------------------
#           Set Constants
#---------------------------------------

    #--constant to convert radius measurements into arcseconds
    r_to_arcsec = 0.0147928866

    #--constants to convert velocity into km/s--
    kpc_to_km = 1000.0 * 3.0857 * 10.0 ** 13.0
    pc_to_km = kpc_to_km / 1000.0
    distance_kpc = 51.4  # Panagia2003
    distance_km = distance_kpc * kpc_to_km
    arcsec_to_radian = np.pi / (3600.0 * 180.0)
    arcsec_to_km = arcsec_to_radian * distance_km
    day_to_s = 24.0 * 3600.0
    yr_to_s = 365.0 * day_to_s
    km_to_cm = 100000

    arcsecdays_to_kms = arcsec_to_km / day_to_s


#---------------------------------------
#        Determine Input Type
#---------------------------------------


#---------------------------------------
#          Convert Lengths
#---------------------------------------

    if fromunit == 'r':
        out = val * r_to_arcsec
        if tounit == 'arcsec':
            return out
        if tounit == 'km':
            return out * arcsec_to_km
        if tounit == 'pc':
            return out * 1000.0 * arcsec_to_km / kpc_to_km
        if tounit == 'cm':
            return out * arcsec_to_km * km_to_cm

    if fromunit == 'pc':
        if tounit == 'r':
            return val * pc_to_km / arcsec_to_km / r_to_arcsec
        if tounit == 'arcsec':
            return val * pc_to_km / arcsec_to_km
        if tounit == 'km':
            return val * pc_to_km / arcsec_to_km
        if tounit == 'cm':
            return val * pc_to_km / arcsec_to_km * km_to_cm

    if fromunit == 'arcsec':
        if tounit == 'r':
            return val / r_to_arcsec
        if tounit == 'pc':
            return val * arcsec_to_km / pc_to_km
        if tounit == 'km':
            return val * arcsec_to_km
        if tounit == 'cm':
            return val * arcsec_to_km * km_to_cm

    if fromunit == 'km':
        if tounit == 'r':
            return val / arcsec_to_km / r_to_arcsec
        if tounit == 'pc':
            return val / pc_to_km
        if tounit == 'arcsec':
            return val / arcsec_to_km
        if tounit == 'cm':
            return val * km_to_cm

    if fromunit == 'cm':
        out = val / km_to_cm
        if tounit == 'r':
            return out / arcsec_to_km / r_to_arcsec
        if tounit == 'pc':
            return out / pc_to_km
        if tounit == 'arcsec':
            return out / arcsec_to_km
        if tounit == 'cm':
            return out * km_to_cm
        if tounit == 'km':
            return out

#---------------------------------------
#          Convert Velocities
#---------------------------------------

    if fromunit == 'r/yr':
        # convert length units
        if tounit == 'km/yr' or tounit == 'km/s':
            out = convert_units(val, 'r', 'km')
        if tounit == 'arcsec/yr' or tounit == 'arcsec/day':
            out = convert_units(val, 'r', 'arcsec')
        # convert time units
        if tounit == 'km/s':
            return out / yr_to_s
        if tounit == 'arcsec/yr' or tounit == 'km/yr':
            return out
        if tounit == 'arcsec/day':
            return (out / yr_to_s) * day_to_s

    if fromunit == 'arcsec/yr':
        # convert length units
        if tounit == 'r':
            out = val / r_to_arcsec
        if tounit == 'km/yr' or tounit == 'km/s':
            out = convert_units(val, 'arcsec', 'km')
        if touniti == 'arcsec/day':
            out = val
        # convert time units
        if tounit == 'r/yr' or tounit == 'km/yr':
            return out
        if tounit == 'km/s':
            return out / yr_to_s
        if tounit == 'arcsec/day':
            return (out / yr_to_s) * day_to_s

    if fromunit == 'km/yr':
        # convert length units
        if tounit == 'r/yr':
            out = convert_units(val, 'km', 'r')
        if tounit == 'arcsec/yr' or 'arcsec/day':
            out = convert_units(val, 'km', 'arcsec')
        if tounit == 'km/s':
            out = val
        # convert time units
        if tounit == 'r/yr' or tounit == 'arcsec/yr':
            return out
        if tounit == 'km/s':
            return out * yr_to_s
        if tounit == 'arcsec/day':
            return (out / yr_to_s) * day_to_s

    if fromunit == 'km/s':
        # convert length units
        if tounit == 'r/yr':
            out = convert_units(val, 'km', 'r')
        if tounit == 'arcsec/yr' or tounit == 'arcsec/day':
            out = convert_units(val, 'km', 'arcsec')
        if tounit == 'km/yr':
            out = val
        # convert time units
        if tounit == 'arcsec/day':
            return out * day_to_s
        else:
            return out / yr_to_s

#----------------------------------------------------------


#----------------------------------------------------------
# convert_time()
# Author: Kari A. Frank
# Date: July 9, 2014
# Purpose: Convert between sn1987a age (days since explosion) and date (year-month-day)
# Usage: convert_sn1987a_time(age,get='date')
#
# Input:
# intime -- sn1987a age (day#) or year-month-day. age may be either
#           int or string (it will be converted to string).
#           year-month-day must be a string of format '1001-12-26'
#
# get   -- string specifying the format of the returned value, if intime is age.
#          allowed values: 'year','month','day','date' (default='date')
#          all but 'date' will return an integer (date returns as a string)
#
# Output:
# - returns int(day#),int(year),int(month),int(day),or str(date)
#
# Usage Notes:
# -
#

def convert_time(intime, get='date'):

    #---------------------------------------
    #           Import Modules
    #---------------------------------------

    import datetime as dt

#---------------------------------------
#        Determine Input Type
#---------------------------------------

    intime = str(intime)
    times = intime.split('-')

    inage = 0
    indate = 0
    if len(times) == 1:
        inage = 1
    if len(times) == 3:
        indate = 1
        get = 'age'
    #-if intime is wrong format-
    if (inage == 0) and (indate == 0):
        print "ERROR: input time must be either age in days or 'year-month-day'. Quitting."
        return 0

    #---Set zero age = Feb. 23, 1987---
    zero_date = dt.datetime(1987, 0o2, 23, 0, 0, 0)


#---------------------------------------
#      Convert from Date to Age
#---------------------------------------

    if indate == 1:

        #--convert to datetime objects--

        #-provided date-
        thedate = dt.datetime.strptime(intime, '%Y-%m-%d')

        #--calculate days in between--
        age = (thedate - zero_date).days

        #--Return age--
        # returned as int
        return age

#---------------------------------------
#      Convert from Age to Date
#---------------------------------------

    if inage == 1:

        #--add age to zero_date--

        newdate = zero_date + dt.timedelta(days=int(intime))

        if get == 'date':
            return newdate.strftime('%Y-%m-%d')
        if get == 'year':
            return newdate.year
        if get == 'month':
            return newdate.month
        if get == 'day':
            return newdate.day

#----------------------------------------------------------
