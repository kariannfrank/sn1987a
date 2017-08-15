#! /usr/bin/env python

###! /usr/astro/ciao-4.6/bin/python

#Author: Kari A. Frank
#Date: June 19, 2016
#Purpose: Use MCMC to fit a broken linear model to the 
#         SN1987A expansion curve.
#Usage: expansion_mcmc.py [--niterations niterations] [--burnin burnin] [--outroot outroot] [--clobber CLOBBER]
#
#Input:
#
# niterations -- specify number of MCMC iterations. default = 10000
#
# burnin      -- specify the number of iterations to discard before
#                calculating the reported final values and plotting 
#                posteriors. default = 1000
#
# outroot     -- string containing optional prefix to output file names,
#                optionally including path (same as in ciao tools).
#                default = pwd+'/expansion_mcmc_burnin-niterations_'
#
# indir       -- directory containing the radius vs age file. default=pwd
#
# save        -- switch to save mcmc results from all iterations
#                to text file. default='no'
#
# read        -- switch to read previous MCMC results from a file instead
#                of running a new MCMC. default='no'.  must also set file.
#
# mcmcfile    -- mcmc results file from previous run (saved with save=yes)
#
# clobber:    specifies whether files should be overwritten
#             if they already exist (same as in ciao tools,
#             default = 'no')
#
#Output:
# - plots of the parameters vs iteration and plots of the parameter posteriors
#
#Usage Notes:
# - assumes that the errors are gaussian and symmetric (takes the larger
#   of the upper and lower errors on each point).  switching to asymmetric
#   errors, e.g. by using piecewise truncated gaussians, would be done by
#   changing the likelihood function, but has not been implemented.

#---------------------------------------
#           Import Modules
#---------------------------------------
import argparse,os,sys,string
import random as rd
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import fitting as fit
sys.path.append('/astro/research/kaf33/Dropbox/Research/Python_Programs/xmc_analysis')
from xmcplot import plot_param_scatters

#---------------------------------------
#   Parse arguments and set defaults
#---------------------------------------

parser = argparse.ArgumentParser(description='Fits a broken linear model to the SN1987A expansion curve.')#,epilog='NOTE: Extra note.')

pwd = os.getcwd()

#parser.add_argument('required_arg',help='Description and usage for required_arg.',default='default for required_arg')

parser.add_argument('--niterations',help='Number of MCMC iterations.',default='10000')

parser.add_argument('--burnin',help='Number of MCMC iterations to discard when calculating final statistics.',default='1000')

parser.add_argument('--outroot',help='Prefix for output plot files.',default='default')

parser.add_argument('--indir',help='Directory containing the radius vs age file.',default=pwd)

parser.add_argument('--save',help='Save parameters from all iteration to file.',default='no')

parser.add_argument('--readonly',help='Read parameters from previous MCMC file.  Must also supply mcmcfile.',default='no')
parser.add_argument('--mcmcfile',help='Read parameters from previous MCMC saved to file.  Ignored if read=no.')


#this argument does not do anything at the moment!
parser.add_argument('--clobber',help='Overwrite existing files.',default='no')

args = parser.parse_args()

args.niterations = int(args.niterations)
args.burnin = int(args.burnin)

#----Check that niterations > burnin----
if args.niterations <= args.burnin:
    print 'Warning: niterations < burnin.  Setting burnin = 0.10*niteraions.'
    args.burnin = np.ceil(0.10*args.niterations)

#----Set file paths----

#sn1987a_dir = '/export/bulk/rice3/kaf33/Chandra_Observations/SN1987A/comparison_dir/imaging/'
#sn1987a_dir = '/Users/kafrank/Dropbox/Research/SN1987A/'
#sn1987a_dir = '/astro/research/kaf33/Dropbox/Research/SN1987A/'
#outdir = sn1987a_dir
#indir = sn1987a_dir
outdir = pwd+'/'
#indir = pwd+'/'

itstr = str(args.burnin)+'-'+str(args.niterations)

if args.outroot == 'default':
    args.outroot = outdir+'expansion_mcmc_'
else:
    args.outroot = args.outroot+'_'

result_file = args.outroot+itstr+'_results.txt'
model_plot_file = args.outroot+itstr+'_fit.pdf'
iteration_plot_file = args.outroot+itstr+'_iterations.pdf'
posterior_plot_file = args.outroot+itstr+'_posteriors.pdf'
scatter_plot_file = args.outroot+itstr+'_scatter.pdf'
postfile = args.outroot+itstr+'_chains.txt'

infile = args.indir+'/sn1987a_radii_fits.txt'

#----Set Useful Constants----

#--convert fitted radii and errors to arcsec--
r_to_arcsec = 0.0147928866

#--unit conversions--
kpc_to_km = 1000.0*3.0857*10.0**13.0
distance_kpc = 51.4 #Panagia2003
distance_km = distance_kpc*kpc_to_km
arcsec_to_radian = np.pi/(3600.0*180.0)
arcsec_to_km = arcsec_to_radian*distance_km
days_to_s = 24.0*3600.0
arcsecdays_to_kms = arcsec_to_km/days_to_s


#---------------------------------------
#           Read in Data
#---------------------------------------

obsids,ages,r0,r0err_low_raw,r0err_upp_raw,counts,deconcounts = np.loadtxt(infile,unpack=True,comments='#',skiprows=1)

#----Calculate Errors----
r0err_low = np.sqrt( r0err_low_raw**2.0 + r0**2.0/counts  )
r0err_upp = np.sqrt( r0err_upp_raw**2.0 + r0**2.0/counts  )

r0err_max = [max(r0err_low[j],r0err_upp[j]) for j in range(len(r0))]

#----Convert Radii to arcsec----
r0 = r0*r_to_arcsec
r0err_low = r0err_low*r_to_arcsec
r0err_upp = r0err_upp*r_to_arcsec
r0err_max = np.array(r0err_max)*r_to_arcsec

#---------------------------------------
# If Read Only, Read from Previous MCMC
#---------------------------------------
if args.readonly == 'yes':
    model_params = np.loadtxt(args.mcmcfile,comments='#')

else: #if readonly=no

#---------------------------------------
#             MCMC Setup
#---------------------------------------

#----Set Parameter Limits----

# to fix a parameter at its initial value, set its 
#  upper limit = lower limit; mcmc() knows
#  what to do with it.

#--slopes in km/s--
    a1_limits = [5000.0,15000.0]
    a2_limits = [1000.0,3000.0]

#--fixed at Judith's values (within errors)--
    #a1_limits = [7539.0-1524.0,7539.0+2139.0]
    #a2_limits = [1592.0-566.0,1592.0+589.0]

#--convert to arcsecs/day--
    a1_limits = [a/arcsecdays_to_kms for a in a1_limits]
    a2_limits = [a/arcsecdays_to_kms for a in a2_limits]

#--intercepts--
#not sure about these
    #b1_limits = [-1.0,1.0]
    b1_limits = [-0.5,0.5]
    #b1_limits = [0.0,0.0]

#--change/break point--
    tb_limits = [5000,7000] 
    #tb_limits = [4500,10000]

#note that the change point tb = (b2-b1)/(a1-a2)

#--combine limits into a single array for passing to mcmc--
    limits = np.array([[0,0],a1_limits,a2_limits,b1_limits,tb_limits])
#extra [0,0] is standin for iteration column, to make all arrays match
# such that a given parameter always has the same column number

#--calculate number of free parameters--
    m = 0
    for p in [1,2,3,4]:
        if limits[p,0] != limits[p,1]: 
            m = m + 1

#----Initialize Empty Arrays----
#use pandas here?

#----Choose Initial Values----
# !! Make sure initial values are inside the limits!!

    a1_zero = 7539.0/arcsecdays_to_kms #values from Judith's paper
    a2_zero = 1592.0/arcsecdays_to_kms 
    b1_zero = np.average(b1_limits)
    tb_zero = 6060

#--combine into array for passing to mcmc--
    zero_params = np.array([0,a1_zero,a2_zero,b1_zero,tb_zero])

#---------------------------------------
#              Run MCMC
#---------------------------------------

    model_params = fit.mcmc(ages,r0,r0err_max,limits,zero_params,args.niterations)

#-------end if readonly=no loop

#---Write All Data to File----

parheads = ['iteration','earlyslope["/day]','lateslope["/d]','intercept["]','changepoint[day]']

if args.save == 'yes':
    np.savetxt(postfile,model_params,header=string.join(parheads))

#---------------------------------------
#          Get Final Values
#---------------------------------------

#----Cut Out Iterations Before Burn-In----
print model_params.shape
final_params = model_params[args.burnin:,:]


#----Calculate Statistics----

#-array to hold statistics (rows = [median,average,stdev])
stats = np.zeros((3,4))

for par in [1,2,3,4]:
    stats[0,par-1] = np.median(final_params[:,par])
    stats[1,par-1] = np.mean(final_params[:,par])
    stats[2,par-1] = np.std(final_params[:,par])

#-convert slopes to velocities-
stats_vel = np.copy(stats)
stats_vel[:,0] = stats_vel[:,0]*arcsecdays_to_kms
stats_vel[:,1] = stats_vel[:,1]*arcsecdays_to_kms

print 'a1,a2,b1,tb'
print stats_vel

if args.readonly == 'no':
    dof = r0.shape[0]-m #number of data points - number of free parameters
    print 'Number free parameters = ',m
    print 'Number data points = ',r0.shape[0]
    print 'Degrees of Freedom = ',dof

#----Write to File----

rfile = open(result_file,'w')
rfile.write('iterations\t{ilow:d}\t{ihigh:d}\n'.format(ilow=args.burnin,ihigh=args.niterations))
rfile.write('stat\tv_early\tv_late\tintercept\tchangepoint\n')
rfile.write('median\t{one:}\t{two}\t{three}\t{four}\n'.format(one=stats_vel[0,0],two=stats_vel[0,1],three=stats_vel[0,2],four=stats_vel[0,3]))
rfile.write('mean\t{one}\t{two}\t{three}\t{four}\n'.format(one=stats_vel[1,0],two=stats_vel[1,1],three=stats_vel[1,2],four=stats_vel[1,3]))
rfile.write('stdev\t{one}\t{two}\t{three}\t{four}\n'.format(one=stats_vel[2,0],two=stats_vel[2,1],three=stats_vel[2,2],four=stats_vel[2,3]))            
rfile.close()


#---------------------------------------
#            Validation
#---------------------------------------


#---------------------------------------
#           Plot Results
#---------------------------------------

#----Plot Expansion Curve with Best Fit Model----
pdffile = PdfPages(model_plot_file)
fig,ax=plt.subplots()

ax.set_xticks(ages,minor=False)
ax.xaxis.grid(True,which='major')
for tick in ax.xaxis.get_major_ticks():
    tick.label.set_fontsize(10)
    tick.label.set_rotation('vertical')
plt.subplots_adjust(bottom=0.13)

plt.xlim(4000,11000)
plt.ylim(0.5,0.9)

#--plot data--
plt.errorbar(ages,r0,yerr=[r0err_low,r0err_upp],fmt='.')
plt.scatter(ages,r0,marker='o')

#--add labels--
plt.xlabel('SN1987A Age [days]')
plt.ylabel('Radius [arcsec]')

#--plot model with median parameters--
medpars = np.array([0,stats[0,0],stats[0,1],stats[0,2],stats[0,3]])
modely = fit.arr_brokenlinear(range(10000),medpars)
plt.plot(range(10000),modely)

#--plot best from mpfit--
mpfile = '/export/bulk/rice3/kaf33/Chandra_Observations/SN1987A/comparison_dir/imaging/expansion_fits/mpfit_results.txt'
mp_results = np.loadtxt(mpfile,comments='#',skiprows=1,usecols=(1,2,3,4))
mppars = np.array([0,mp_results[0,0]/arcsecdays_to_kms,mp_results[0,1]/arcsecdays_to_kms,mp_results[0,2],mp_results[0,3]])
modely = fit.arr_brokenlinear(range(10000),mppars)
plt.plot(range(10000),modely,linestyle='--',color='green')

pdffile.savefig()
plt.close()
pdffile.close()

#----Plot Parameters vs. Iteration (Trace Plots)----
# ONLY USES ITERATIONS AFTER BURNIN

par_names = ['Iteration','Early Slope (km/s)','Late Slope (km/s)','Intercept (arcsec)','Change Point (days)']

pdffile = PdfPages(iteration_plot_file)

#-only plot every markit-th iteration--
ilog = np.log10(args.niterations)
markit= 10.0**(ilog-1)

#--one plot per page--
x = model_params[:,0]

for p in [1,2,3,4]:
    
    if (p == 1) or (p == 2): #convert slopes to velocities
        y = model_params[:,p]*arcsecdays_to_kms
    else:
        y = model_params[:,p]

    fig,ax=plt.subplots()        
    plt.plot(x,y,'.',rasterized=True)
    plt.xlabel(par_names[0])
    plt.ylabel(par_names[p])
    ymin,ymax = ax.get_ylim()
    xmin,xmax = ax.get_xlim()
    #print ymin,ymax,xmin,xmax
    # plot vertical bar indicating burnin
    plt.plot([args.burnin,args.burnin],[ymin,ymax],linestyle='-',color='green')

    # plot horizontal bar indicating the median value
    med = stats_vel[0,p-1]
    plt.plot([xmin,xmax],[med,med],linestyle='-',color='red')

    # fill range corresponding to 1sigma
    sig = stats_vel[2,p-1]
    plt.fill_between([xmin,xmax],[med-sig,med-sig],[med+sig,med+sig],color='red',alpha=0.1)
        
    pdffile.savefig()
    plt.close()

pdffile.close()

#----Plot Posteriors----
pdffile = PdfPages(posterior_plot_file)
nbins = 30

for p in [1,2,3,4]:

    if (p == 1) or (p == 2): #convert slopes to velocities
        y = final_params[:,p]*arcsecdays_to_kms
    else:
        y = final_params[:,p]
        
    y_hist, x_bin = np.histogram(y,bins=nbins)
    fig,ax=plt.subplots()
    plt.bar(x_bin[:-1],y_hist,width=x_bin[1]-x_bin[0])
    plt.xlabel(par_names[p])
    ymin,ymax = ax.get_ylim()

    #-plot median and 1sigma range-
    med = stats_vel[0,p-1]
    sig = stats_vel[2,p-1]
    plt.plot([med,med],[ymin,ymax],color='red',linestyle='-')
    plt.fill_betweenx([0.0,ymax],[med-sig,med-sig],[med+sig,med+sig],color='red',alpha=0.2)

    pdffile.savefig()
    plt.close()

pdffile.close()

#----Pair Plots (Parameter Scatter Plots)----

scatter_names = par_names[1:4]
scatter_cols = [1,2,3,4]
scatter_tup = zip(scatter_cols,scatter_names)
pdffile = PdfPages(scatter_plot_file)
#plot_param_scatters(scatter_tup,final_params,pdffile)

for p in scatter_cols: #one plot per page

    if min(final_params[:,p]) != max(final_params[:,p]): #skip if fixed

        if (p == 1) or (p == 2): #convert slopes to velocities
            x = final_params[:,p]*arcsecdays_to_kms
        else:
            x = final_params[:,p]

        fig,ax=plt.subplots()
        for py in scatter_cols:
            if (py != p) and (min(final_params[:,py]) != max(final_params[:,py])):
                if (py == 1) or (py == 2): #convert slopes to velocities
                    y = final_params[:,py]*arcsecdays_to_kms
                else:
                    y = final_params[:,py]
                plt.plot(x,y,'.',rasterized=True)#,markevery=markit)
                plt.xlabel(par_names[p])
                plt.ylabel(par_names[py])
                pdffile.savefig()
                plt.close()

pdffile.close()

#---------------------------------------
#       Print out final status
#---------------------------------------



