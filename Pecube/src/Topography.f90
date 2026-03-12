subroutine Synthetic_topography(topomag,topoRefflag, z, topooffset, istep, nstep, nsurf,&
                    X,lambda,Final_offset,phase_array,Final_amp,Final_phase,topo_ref_custom,nx,ny)
                    
    ! Generate a synthetic topography (sinusoidal) for each time step.
    ! topography is controlled by:
    !  - topomag: amplification of relief at time(istep)
    !  - topoRefflag: flag for reference elevation from which to apply amplification
    !  - topooffset: vertical offset on z (m) at time(istep)
    !  - X: distance along longitude (m)
    !  - lambda: topography wavelength (m)
    
    ! output is z: the topography at time(istep)
     implicit none
     
     integer nsurf, istep, nstep, ij, topoRefflag,nx,ny
     double precision topoRef, pi, lambda
     double precision topooffset(nstep+1), topomag(nstep+1), phase_array(nstep+1)
     double precision Zpresetn(nsurf), z(nsurf), X(nsurf),Final_offset,Final_amp,Final_phase
     double precision FinalTopo(nsurf), topo_ref_custom
     
     pi = 3.1415927
     ! Final topography
     do ij=1,nsurf
        FinalTopo(ij) = (Final_amp/2*sin(2*pi* (X(ij)+ Final_phase) / lambda) + Final_offset)
     enddo  
     ! Get reference topography from final topography
     if (topoRefflag.eq.0) then !from sea level
      topoRef = 0
    elseif (topoRefflag.eq.1) then !from min elevation
      topoRef = minval(FinalTopo(:))
    elseif (topoRefflag.eq.2) then !from max elevation
      topoRef = maxval(FinalTopo(:))
    elseif (topoRefflag.eq.3) then !from mean elevation
      topoRef = sum(FinalTopo(:))/(nx*ny)
    elseif (topoRefflag.eq.4) then !from custom elevation
      topoRef = topo_ref_custom
    elseif (topoRefflag.eq.5) then !from mean elevation but varying amplitude
      topoRef = sum(FinalTopo(:))/(nx*ny)
    endif !end Maxime
            
            
     ! topography at time(istep)
     do ij=1,nsurf
        z(ij) = topoRef - (topomag(istep) * (topoRef - (Final_amp/2*sin(2*pi* (X(ij)+ Final_phase + phase_array(istep)) / lambda) +&
        Final_offset))) + topooffset(istep)
     enddo  
    
    
                
end subroutine Synthetic_topography


subroutine Headward_propagation(topomag,topoRef, z, topooffset, istep, nstep, amplification, nsurf, Di_percent,&
                    time_incision_start,time_incision_stop, Zinit, Zpresent, time, time_prev, Zmin, Zmax,tauH)

    ! Routine to get a wave propagating upvalley

    ! W : propagation rate [m/Myr]
    ! Di : relative depth of incision (%)
    ! time_incision_start : Onset timing of wave propagation [Ma]
    ! time_incision_stop : ending time of wave propagation [Ma]
    ! Zinit: initial elevation before incision 
    ! Zpresent: present topography (elevation)

    ! Zmin: min elevation at time of wave propagation start
    ! Zmax: max elevation at time of wave propagation start
    ! amplification = array of size nx*ny where to store spatial variation of topo amplitude
    ! z: new topography
    ! time_prev > time

    ! Authors: Maxime Bernard

    implicit none

    integer nsurf, istep, nstep, ij
    double precision Di, W, time_incision_start,time_incision_stop, time_diff, tauH
    double precision topoRef, time, time_prev, Di_percent, wave_front, wave_front_prev
    double precision minR, maxR, Zmin, Zmax,  check, Amp, dt, ftime, ftime_prev
    double precision topooffset(nstep+1), topomag(nstep+1), tol
    double precision amplification(nsurf), rate(nsurf), Zpresent(nsurf), z(nsurf), Zinit(nsurf)
    double precision time_after_wave, time_after_wave_array(nsurf)

    !!!!! Parameters !!!!
    Di = Di_percent/100
    tol = 1e-8
    time_diff = time_incision_start - time_incision_stop
    if (time_prev.eq.time_incision_start) time_prev = time_incision_start+1e-8
    ! Get wave propagation rate (km/Myr)
    W = (Zmax - Zmin) / (time_incision_start - time_incision_stop) ! for constant wave propagation rate (i.e., linear increase of wave front)
    
    !!!!! Handle non-linear change of wave propagation rate !!!!!!
    ftime = 0.0
    if (time.ge.time_incision_start) then ! if the wave propagation did not start yet the wave front stay at minimum elevation
        ftime = 0.0
    elseif (time.le.time_incision_stop) then ! if the wave propagation finished  the wave front is at maximum elevation
        ftime = 1.0
    else
        ftime = 1 - (time - time_incision_stop) / ( time_incision_start - time_incision_stop)
    endif
    ! for previous time step
    ftime_prev = 0.0
    if (time_prev.ge.time_incision_start) then ! if the wave propagation did not start yet the wave front stay at minimum elevation
        ftime_prev = 0.0
    elseif (time_prev.le.time_incision_stop) then ! if the wave propagation finished  the wave front is at maximum elevation
        ftime_prev = 1.0
    else
        ftime_prev = 1 - (time_prev - time_incision_stop) / ( time_incision_start - time_incision_stop)
    endif
    ! Apply non-constant propagation rate
    if (tauH.ne.0.0) then
        ftime = ( 1 - exp(-ftime * tauH / ( time_incision_start - time_incision_stop))) / &
                (1 - exp(-tauH / ( time_incision_start - time_incision_stop)))
        
        ftime_prev = ( 1 - exp(-ftime_prev * tauH / ( time_incision_start - time_incision_stop))) / &
                (1 - exp(-tauH / ( time_incision_start - time_incision_stop)))
    endif
    
    !!!!! Get the elevation of wave front !!!!
    wave_front = Zmin + W*(time_incision_start - time_incision_stop) * ftime
    wave_front_prev = Zmin + W*(time_incision_start - time_incision_stop) * ftime_prev
    
    ! print *,'time = ', time
    ! print *,'time_prev = ', time_prev
    ! print *,'time_incision_start = ', time_incision_start
    ! print *,'time_incision_stop = ', time_incision_stop
    ! print *,'W = ', W
    ! print *,'Di = ', Di
    ! print *,'zmin = ', Zmin
    ! print *,'tauH = ', tauH
    ! print *,'wave_front = ', wave_front
    ! print *,'wave_front_prev = ', wave_front_prev
    ! print *,'ftime = ', ftime
    ! print *,'ftime_prev = ', ftime_prev
    ! print *,' '
    
    ! Do propagation (get new amplitude)
    if (time.eq.0.0) then
        amplification(:) = 1.0 ! this should be equal to amplification at the ending time of wave propagation
    elseif (time.lt.time_incision_start.and.time.ge.time_incision_stop) then
        minR = wave_front_prev
        maxR = wave_front
        do ij=1,nsurf
            ! Topo in wave mask
            if (Zinit(ij).ge.minR.and.Zinit(ij).lt.maxR) then
                ! amplification(ij) = amplification(ij) + (Di*(1-amplification(ij)))
                ! We assume the incision during this time step is only due to the wave
                ! applied to the amp at time = time_prev
                ! Amplification now
                amplification(ij) = amplification(ij) + (Di*(1-amplification(ij)))
            ! topo above wave
            elseif (Zinit(ij).ge.maxR) then
                amplification(ij) = topomag(istep+1)
            ! topo below wave
            elseif (Zinit(ij).lt.minR) then !Below wave, rate of incision = regional rate of incision
                !time_after_wave = (Zinit(ij) - Zmin) / W + time_incision_start
                rate(ij) = (topomag(istep+1) - topomag(istep)) / (time_prev - time)
                !dt = time_after_wave - time
                amplification(ij) = amplification(ij) + rate(ij)* (time_prev - time)
            endif
        enddo
    elseif (time.lt.time_incision_stop) then
        ! All elevation are below the wave front and incised at regional rate
        !time_after_wave_array = (Zinit - Zmin) / W + time_incision_start
        rate = (topomag(istep+1) - topomag(istep)) / (time_prev - time)
        amplification = amplification + rate * (time_prev - time)
        ! Ensure Amplification <= 1
        do ij=1,nsurf
            if (amplification(ij).gt.1) amplification(ij) = 1
        enddo
    else ! if time >= time_incision start
        amplification(:) = topomag(istep+1)
    endif

    ! Get new elevation
    z(:) = topoRef - (amplification * (topoRef - Zpresent)) + topooffset(istep+1)

end subroutine Headward_propagation


subroutine get_initial_topo(topomag, topooffset, timek, nstep, topoRef, Zpresent, nsurf,&
                            time_incision_start,time_incision_stop, NewTopo)

    ! Routine to get initial topography before wave propagation

    ! topoRef: Reference topography
    ! Zpresent: Present-day topography
    ! time_incision : Onset timing of wave propagation [Ma]
    ! timek: model time array (time increases)
    ! maxTime: Starting time of modelling (Ma)

    implicit none

    integer nsurf,nstep,istep,ind
    double precision topoRef, topomag(nstep), topooffset(nstep), timek(nstep+1), Zpresent(nsurf)
    double precision topomag2(nstep+1), topooffset2(nstep+1)
    double precision time, rate, NewTopo(nsurf), amplification, Toffset
    double precision time_incision_stop, time_incision_start, maxTime, tol
    

    tol = 1e-8
    maxTime = maxval(timek(:))
    ! Consider last time step for amplification and offset
    topomag2 = 0.0
    topomag2(1:nstep) = topomag
    topomag2(nstep+1) = 1.0
    topooffset2 = 0.0
    topooffset2(1:nstep) = topooffset
    topooffset2(nstep+1) = 0.0

    ! First find index time incision start
    ind = nstep +1 !Default is last time step (i.e., t = 0)
    do istep=1,nstep
        time = timek(istep)
        ! print *, 'timek = ', timek(istep)
        ! print *, 'maxTime = ', maxTime
        if (time.lt.time_incision_start) then ! We take the first index that is lower than tstart
            ind = istep
            goto 999
        endif
    enddo

    999 continue

    ! Amplification
    if (ind.eq.1) then 
        ! If we start wave propagation at the start of the modelling 
        ! then amplification equal the user input amplification for the initial topography
        rate = 0.0
        amplification = topomag2(1)
        ! Offset
        Toffset = topooffset2(1) 
    else
        ! We search for the amplification value at time of wave propagation start
        rate = (topomag2(ind)+tol - topomag2(ind-1)+tol) / (timek(ind-1) - timek(ind))
        amplification = topomag2(ind-1) + rate * (timek(ind-1) - time_incision_start)
        ! Offset
        rate = (topooffset2(ind-1)+tol - topooffset2(ind)+tol) / (timek(ind-1) - timek(ind))
        Toffset = topooffset2(ind-1) - rate * (timek(ind-1) - time_incision_start)
    endif
    
    ! Get topo at the onset of incision
    NewTopo(:) = topoRef - (amplification *(topoRef-Zpresent)) + Toffset
    
end subroutine get_initial_topo