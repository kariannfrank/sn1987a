def flux_function(tR,eta=7.3,F0=0.2,D=0.018):
    """
    Return flux (as function time according to the function in Park2005

    Parameters
    
    tR : tuple of numerical arrays or
      (t,R_t), where t and R_t are as defined below. these are the 
      independent variables

    t : numerical array or list
      days since t0, beginning of interaction with inner HII region
    
    t0 : float
      day (since SN) of beginning of interaction with inner HII region

    R_ring : float
      inner radius of the equatorial ring (where the filling factor
      of the clumps goes to unity)

    eta : float
      ratio nr/n0, where nr is density of the ring and n0 is density
      of the HII region.  should be greater than 1.  eta = 100 is good 
      initial guess.

    D : float
      characteristic scale height of the exponential filling factor of
      the clumps. has same units as R_ring and R_t

    R_t : np.array or numerical list
      array of radius values for each time t. must be same length as t
      and have the same order (e.g. R[5] must correspond to t[5])

    F0 : float
      flux at time t0


    Notes
    
    - Fluxes are typically in units of 10e-13 ergs/s/cm^2
    - Best fit values from Park2005 are a good place to start for initial 
      guesses: D = 0.018", F0 = 0.2e-13 ergs/s/cm^2, eta = 7.3
    - For now, all parameters except D, F0, and eta are fixed, as done in 
      Park2005.  However, I use a more up-to-date value for R_ring.  This
      requires R_t and t to be in units of arcsec and days since explosion,
      respectively.
    - Default values are from the best fit in Park2005

    """

    #--import modules--
    import numpy as np

    #--set fixed parameters--
    t0 = 1200.0
    R_ring = 0.74

    #--calculate fluxes--
    t = tR[0].astype(float)
    R_t = tR[1]
    tau = t/t0
    f = F0 * (R_t**3.0) * (tau-1.0)**0.4 * ( 1 + (eta**2.6 - 1) 
                                             * np.exp((R_t-R_ring)/D) )

    return f


