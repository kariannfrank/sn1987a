def power_law(x,pars)
    """
    Return Y(x), where y = a + b*x^beta

    Parameters
    
    X : numerical array or list
      array of X values

    pars : numerical list
      list of parameters for the power law function:
      pars = [a,b,beta]

    """

    #--import modules--
    import numpy as np

    slope2=slope1
    slope3=slope1

    #--warn if changepoint1 > changepoint2--
    if changepoint1 > changepoint2:
        print 'Warning: changepoint1 > changepoint2. Switching their order.'
        c1temp = changepoint1
        changepoint1 = changepoint2
        changepoint2 = c1temp

    pars = [intercept,slope1,changepoint1,slope2,changepoint2,slope3]
    print 'pars = ',pars

    #--calculate the second and third segment's intercepts--
    b2 = pars[2]*(pars[1]-pars[3]+pars[0])
    b3 = pars[4]*(pars[3]-pars[5]+b2)

    #--calculate y--
    
    y = np.empty_like(x)
    y[ x <= pars[2] ] = pars[1]*x[ x<= pars[2] ]+pars[0]
    #y[ np.where((x > pars[2]) & (x <= pars[4])) ] = pars[3]*x[ np.where((x > pars[2]) & (x <= pars[4])) ]+b2
    y[ (x > pars[2]) & (x <= pars[4]) ] = pars[3]*x[ (x > pars[2]) & (x <= pars[4]) ]+b2
    y[ x > pars[4] ] = (pars[5]*x[ x > pars[4] ]+b3)

    return y


