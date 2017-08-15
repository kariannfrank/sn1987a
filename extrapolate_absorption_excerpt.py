##this section cut from spectra_results.py (comes just before the plots section)

#---------------------------------------
#---------------------------------------
# Extrapolate Contamination Absorption
#---------------------------------------
#---------------------------------------

#--get indices of first and last observations to include--
start_i = list(b_age).index(8000)
stop_i =  list(b_age).index(9165) + 1

if args.dofit == 'yes':

#---------------------------------------
#      Set Up X,Y
#---------------------------------------

#b_abs_soft = b_cflux_soft/b_flux_soft
#b_abs_soft12 = b_cflux_soft12/b_flux_soft12
#b_abs_soft23 = b_cflux_soft23/b_flux_soft23
#b_abs_soft13 = b_cflux_soft13/b_flux_soft13
#b_abs_hard = b_cflux_hard/b_flux_hard

    print 'first to last ages in extrapolation: ',b_age[start_i:stop_i]

#----Fit 0.5-2.0 keV----

#--set initial values--
    soft_initial = [0,-8.6*10.0**(-5.0),1.1,0,0]

#--set up limits--
    soft_lims = [[0,0],[11*soft_initial[1],0.0],[0.7,5.0],[0,0],[0,0]]

#--run mcmc--
    soft_params = fitting.mcmc(b_age[start_i:stop_i],b_abs_soft[start_i:stop_i],b_abs_soft_low[start_i:stop_i],soft_lims,soft_initial,niterations=10000000,func='linear')

#--cut burnin--
    burnin = 1000

#--calculate statistics--
    soft_stats = np.zeros((3,2)) #rows=median,average,stdev,columns=params

    for par in [1,2]:
        soft_stats[0,par-1] = np.median(soft_params[burnin:,par])
        soft_stats[1,par-1] = np.mean(soft_params[burnin:,par])
        soft_stats[2,par-1] = np.std(soft_params[burnin:,par])

        print 'soft_stats = ',soft_stats

#--plot posteriors--
        pnames = ['iterations','slope','intercept']
        pdffile2 = PdfPages(posteriorfile)
        fitting.plot_posteriors(soft_params[burnin:,:],pdffile2,names=pnames)
        pdffile2.close()


#---fit with polyfit---

if dofit == 'yes':
    poly_x = range(int(b_age[start_i]),10001)
    soft_poly_a,soft_poly_b = np.polyfit(b_age[start_i:stop_i],b_abs_soft[start_i:stop_i],1)
    soft_poly_y = [soft_poly_a*x+soft_poly_b for x in poly_x]

    soft12_poly_a,soft12_poly_b = np.polyfit(b_age[start_i:stop_i],b_abs_soft12[start_i:stop_i],1)
    soft12_poly_y = [soft12_poly_a*x+soft12_poly_b for x in poly_x]

    soft13_poly_a,soft13_poly_b = np.polyfit(b_age[start_i:stop_i],b_abs_soft13[start_i:stop_i],1)
    soft13_poly_y = [soft13_poly_a*x+soft13_poly_b for x in poly_x]

#---print extrapolated ratios----
    i15809 = list(b_obsid).index(15809)
    i14697 = list(b_obsid).index(14697)
    ext_ages = b_age[i14697:i15809+1]

    ext_polyratio_soft = fitting.arr_linear(b_age[i14697:i15809+1],[0,soft_poly_a,soft_poly_b])
    ext_polyratio_soft12 = fitting.arr_linear(b_age[i14697:i15809+1],[0,soft12_poly_a,soft12_poly_b])
    ext_polyratio_soft13 = fitting.arr_linear(b_age[i14697:i15809+1],[0,soft13_poly_a,soft13_poly_b])

    ext_flux_soft = b_cflux_soft[i14697:i15809+1]/ext_polyratio_soft
    ext_flux_soft12 = b_cflux_soft12[i14697:i15809+1]/ext_polyratio_soft12
    ext_flux_soft13 = b_cflux_soft13[i14697:i15809+1]/ext_polyratio_soft13

    print 'new soft fluxes = ',ext_flux_soft
    print 'new soft12 fluxes = ',ext_flux_soft12
    print 'new soft13 fluxes = ',ext_flux_soft13

#obs15809_soft_polyratio = fitting.linear(b_age[i15809],[0,soft_poly_a,soft_poly_b])
#obs15809_soft12_polyratio = fitting.linear(b_age[i15809],[0,soft12_poly_a,soft12_poly_b])
#obs15809_soft13_polyratio = fitting.linear(b_age[i15809],[0,soft13_poly_a,soft13_poly_b])
#print 'new 15809 ratios = ',obs15809_soft_polyratio,obs15809_soft12_polyratio,obs15809_soft13_polyratio
#new_15809_soft = b_cflux_soft[i15809]/obs15809_soft_polyratio
#new_15809_soft12 = b_cflux_soft12[i15809]/obs15809_soft12_polyratio
#new_15809_soft13 = b_cflux_soft13[i15809]/obs15809_soft13_polyratio
#print 'new 15809 fluxes = ',new_15809_soft,new_15809_soft12,new_15809_soft13
