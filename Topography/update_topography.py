# -*- coding: utf-8 -*-
"""
Created on Fri Nov 24 10:09:29 2023

This module gather all functions/classes related to the topography evolution for
Pecube.

@author: maxime - maxime.bernard@uni-potsdam.de
"""
from math import *
import numpy as np




#------------------------------------------------------------
def get_initial_topo(topoRef, Zpresent, AmpArray, OffsetArray, TimeArray, Start_incision):
    """ Get the elevations at the timing of incision start."""
    
    # When does incision start
    index = [i  if TimeArray[i] < Start_incision else None for i in range(len(TimeArray))]
    index = np.where(TimeArray < Start_incision)
    index = index[0][0] # first index < time start

    # Get amplification and offset values
    if index == 0:
        amplification = AmpArray[0]
        Offset = OffsetArray[0]
    else:
        rate = (AmpArray[index] - AmpArray[index-1]) / (TimeArray[index-1] - TimeArray [index])
        amplification = AmpArray[index-1] + rate * (TimeArray[index-1] - Start_incision)
        # Offset
        rate = (OffsetArray[index-1] - OffsetArray[index]) / (TimeArray[index-1] - TimeArray [index])
        offset = OffsetArray[index-1] - rate * (TimeArray[index-1] - Start_incision)
    
    Initial_topo = topoRef - (amplification * (topoRef - Zpresent)) + offset
    
    Zmin = np.min(Initial_topo)
    Zmax = np.max(Initial_topo)

    return Initial_topo, Zmin, Zmax, amplification, offset

#------------------------------------------------------------
def Headward_propagation(topoRef, Zpresent, Zinit, amplification, offset, istep, timeArray,
                Start_incision,Stop_incision,Depth_incision,Tauh, Zmax, Zmin):
    """ Do headward propagation of erosion.
        amplification and offset are the array provides by the user.

    """ 
    # Parameters
    AmpArray = Zpresent * 0
    OffsetArray = Zpresent * 0
    Di = Depth_incision / 100
    time_diff = Start_incision - Stop_incision
    
    # Propagation rate
    W = (Zmax - Zmin) / time_diff
    # Loop time until istep
    for t in range(istep+1):
        time = timeArray[t]
        if t == 0:
            AmpArray = Zpresent * 0 + amplification[t]

        else:
            time_prev = timeArray[t-1]
            if time_prev == Start_incision:
                time_prev = time_prev + 1e-8
        
            # Handle non-linear change of wave propagation
            ftime = 1 - (time - Stop_incision) / time_diff
            if time <= Stop_incision:
                ftime = 1.0
            elif time >= Start_incision:
                ftime = 0.0
            # previous time step
            ftime_prev = 1 - (time_prev - Stop_incision) / time_diff
            if time_prev <= Stop_incision:
                ftime_prev = 1.0
            elif time_prev >= Start_incision:
                ftime_prev = 0.0
                
            # non-linear change
            if Tauh != 0.0:
                ftime = (1 - exp(-ftime * Tauh / time_diff)) / (1 - exp(-Tauh / time_diff))
                ftime_prev = (1 - exp(-ftime_prev * Tauh / time_diff)) / (1 - exp(-Tauh / time_diff))
            
            # Wave front elevation
            wave_front = Zmin + W * time_diff * ftime
            wave_front_prev = Zmin + W * time_diff * ftime_prev
            
            # Do propagation
            if time == 0.0:
                AmpArray = Zpresent * 0 + 1.0
                
            elif time < Start_incision and time >= Stop_incision:
                # find elevations in the wave front area
                indexes =  np.where((Zinit >= wave_front_prev) & (Zinit < wave_front))
                AmpArray[indexes] = AmpArray[indexes] + (Di * (1 - AmpArray[indexes]))
                # find elevation above wave front
                indexes =  np.where(Zinit >= wave_front)
                AmpArray[indexes] = amplification[t]
                # find elevations below wave front prev
                indexes = np.where(Zinit < wave_front_prev)
                rate = (amplification[t] - amplification[t-1]) / (time_prev - time)
                AmpArray[indexes] = AmpArray[indexes] +  rate * (time_prev - time)
            elif time < Stop_incision:
                rate = (amplification[t] - amplification[t-1]) / (time_prev - time)
                AmpArray = AmpArray +  rate * (time_prev - time)
                indexes = np.where(AmpArray > 1)
                AmpArray[indexes] = 1
            
            else:
                AmpArray = Zpresent * 0 + amplification[t]
                
            
    # Get new topography
    z = topoRef - (AmpArray * (topoRef - Zpresent)) + offset[istep]

    return z