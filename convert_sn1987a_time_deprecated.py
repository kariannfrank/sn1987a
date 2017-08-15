#Author: Kari A. Frank
#Date: July 9, 2014
#Purpose: Convert between sn1987a age (days since explosion) and date (year-month-day)
#Usage: convert_sn1987a_time(age,get='date')
#
#Input:
# intime -- sn1987a age (day#) or year-month-day. age may be either 
#           int or string (it will be converted to string). 
#           year-month-day must be a string of format '1001-12-26'
# 
# get   -- string specifying the format of the returned value, if intime is age.
#          allowed values: 'year','month','day','date' (default='date')
#          all but 'date' will return an integer (date returns as a string)
#
#Output:
# - returns int(day#),int(year),int(month),int(day),or str(date)
#
#Usage Notes:
# - 
#   

def convert_sn1987a_time(intime,get='date'):

#---------------------------------------
#           Import Modules
#---------------------------------------

#---------------------------------------
#        Determine Input Type   
#---------------------------------------

    intime = str(intime)
    times = intime.split('-')

    #-intime = age
    if len(times) == 1:
        inage = 1
    if len(times) == 3:
        indate = 1
    #-if intime is wrong format-
    if (inage == 0) and (indate == 0): 
        print "ERROR: input time must be either age in days or 'year-month-day'. Quitting."
        return 0
        

#---------------------------------------
# 
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

        if get == 'date': return date
        if get == 'year': return int(year)
        if get == 'month': return int(month)
        if get == 'day': return int(day)
