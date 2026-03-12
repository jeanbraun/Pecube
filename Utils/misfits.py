"""
    This file gather all functions and classes related to misfit calculation.

    @author: Maxime Bernard

"""

import math

def compute_misfit(flag,obs,pred,error,errorpred):
    """
    This function computes the misfit value for observed vs predicted data.
    obs and pred are 1D array containing single thermochronometer data.
    flag: flag for misfit function
    
    """
    # Check if number of observation matches number of predictions
    LL = 0
    sigma = 0
    N = 0
    # if errorpred > 0:
    #     error = error + pred*(errorpred/100)
    print("Misfit function = ", flag)
    if len(obs) == len(pred):
        if flag == 1 or flag == 2:
            # Compute with chi-square
            for i in range(len(pred)):
                if obs[i] > 0 and pred[i]>0 and error[i] > 0:
                    N+=1
                    error[i] += 1e-10
                    # Sum differences between predicted and observed ages
                    # and divide by Observed error
                    try:
                        sigma += (obs[i]-pred[i])**2/(error[i]**2+errorpred[i]**2)
                    except ZeroDivisionError:
                        print("Division by zero in misfit calculation.")
                        LL = N = sigma = 0
                        return
                    # val = math.log(2*math.pi)/2 + math.log(error[i]) + 0.5 * ((obs[i]-pred[i])/error[i])**2 # Braun et al. (2012)
                    val = 0.5 * ((obs[i]-pred[i])**2/ (error[i]**2 + errorpred[i]**2)) # Braun et al. (2012)
                    LL += val
                    # print(obs[i],pred[i],error[i],val,abs(obs[i]-pred[i]))
            LL = -LL
        elif flag == 3:
            # Compute with L1-norm
            for i in range(len(pred)):
                if obs[i] > 0 and pred[i]>0 and error[i] > 0:
                    N+=1
                    error[i] += 1e-10
                    # Sum differences between predicted and observed ages
                    # and divide by Observed error
                    try:
                        sigma += (abs(obs[i]-pred[i])/ math.sqrt(error[i]**2 + errorpred[i]**2))
                        LL +=  -math.log(2*error[i]) - (abs(obs[i]-pred[i])/(math.sqrt(error[i]**2 + errorpred[i]**2)))
                    except ZeroDivisionError:
                        print("Division by zero in misfit calculation.")
                        LL = N = sigma = 0
                        return
    else:
        print('Number of observations does not match number of prediction. Failed to compute misfit.')
    
    return LL, N, sigma