subroutine calculate_ages (nrec,nsurf,nz,istep,age1,age2,age3,age4,age5,age6,age7,age8, &
                ftdist,fmean,dfmean,tprevious,ftlflag,ageflag,iproc,nd,saveTt,grainsize,&
                                 ASize,Uppm,Thppm,kinFTL_AFT,kinFTL_AHe,Init_FTL_Model,Init_FTL_value,kinFTLID,rhoST,&
                                 D0_spec,Ea_spec,RDmodel,Alpha_Ejec_Flag,data_ADAM_array,D0z,Eaz,RDmodelz,&
                                 ZUppm,ZThppm,ZSize,sampleID,EdgeAge,duration,temp_heat,He43_flag,D0k,Eak,&
                                 D0b,Eab,D0m,Eam,D0h,Eah)
 !__________________________________________________________________________
 !This file has been changed to use ketcham program. 
 !original version is found in calculate_ages_old.f90 in this same folder
 !by VKP
 ! Calculate thermochronological ages. Called by Pecube.                                
 !__________________________________________________________________________

  implicit none

  integer nsurf,nrec,nz,istep,i,ij,irec,k,ftlflag,ageflag(9)
  real*4 age1(nsurf),age2(nsurf),age4(nsurf),age5(nsurf),age3(nsurf) 
  real*4 age6(nsurf),age7(nsurf),age8(nsurf)
  real*4 EdgeAge(nsurf)
  real*4,dimension(:),allocatable::depth,ztime,ztemp,ztime_3,ztemp_3    
  real*4 zjunk,timenow,ftld(17),ftldmean,ftldsd,grainsize(nsurf)
  real*4 Uppm(nsurf),Thppm(nsurf),ASize(nsurf)
  real*4 ZUppm(nsurf),ZThppm(nsurf), ZSize(nsurf)
  character *50 sampleID(nsurf)
  real*8 D0_spec,Ea_spec !Maxime
  real*8 D0z,Eaz
  double precision Init_FTL_value, kinFTL_AFT(nsurf),kinFTL_AHe(nsurf)
  real*8 D0k,Eak,D0b,Eab,D0m,Eam,D0h,Eah
  integer RDmodel,Alpha_Ejec_Flag,RDmodelz !Maxime
  double precision fdist(200),oldest_age,final_age,fmean(nsurf),counter,dfmean(nsurf)
  double precision ftdist(20,nsurf),tprevious, rhoST, fmeanp,dfmeanp
  integer aftmodel, Init_FTL_Model, kinFTLID
  integer iproc,nd
  double precision,dimension(:,:),allocatable::mem
  logical saveTt
  character*4 c(9999)
  character*10000 line    
  double precision,dimension(:),allocatable::ztime_double,ztemp_double,dummy 
  real ageTime_a(nrec), ageTime_z(nrec)   
  ! Maxime - 4He/3He thermochronometer
  double precision,dimension(:),allocatable:: He_ratio,TotHe3_released,He_ratioEjec
  double precision,dimension(:),allocatable:: StepAge,He43obs,SF3Heobs
  integer ndur, ns43,He43_flag
  double precision duration(20),temp_heat(20)
  double precision,dimension(:),allocatable:: released, agereleased
  real*4 FT,Fedge,S, RsRbAlpha
  ! kinetic model parameters
  double precision,dimension(3,1000000)::data_ADAM_array
  real,dimension(:),allocatable::rho_red
      
      
  do i = 1, 9999
    write (c(i),'(i4)') i
    if (i.lt.10) c(i)(1:3)='000'
    if (i.lt.100) c(i)(1:2)='00'
    if (i.lt.1000) c(i)(1:1)='0'
  enddo

  aftmodel=ftlflag
      
  allocate (ztime(nrec),ztemp(nrec),depth(nrec),ztime_3(nrec),ztemp_3(nrec))

  if (saveTt) allocate (mem(nrec,0:nsurf)) 
  age1=0.
  age2=0.
  age3=0.
  age4=0.
  age5=0.
  age6=0.
  age7=0.
  age8=0.
  ndur = 20

  ftdist = 0. ! bug fix by Yuntao May 4, 2020                                        

  !-------------------- Get Thermal histories ---------------------------
  do i=1,nsurf

    if (iproc.eq.0.and.nd.eq.0) call screen_counter (i,nsurf,0)

    ij=(i-1)*nz+1
    ! Read thermal history
    do irec=1,nrec
      read (100+istep,rec=irec) ztime(irec),(zjunk,k=1,i-1),ztemp(irec), &
                                    (zjunk,k=i+1,nsurf), &
                                    (zjunk,k=1,i-1),depth(irec)
    enddo

    timenow=ztime(nrec) ! timenow is minimum time (or time of output)
    ztime=ztime-timenow

    if (saveTt) then
      mem(:,i)=ztemp
      mem(:,0)=ztime
    endif

    ! ---------------- Apatite Fission Track ages -------------------------
    if (aftmodel.ge.1) then
      ! reverse time and temperature array (from 0 to x Ma)
      do irec=1,nrec
        ztime_3(irec)=ztime(nrec-irec+1) !MOD VKP
        if (ztemp(nrec-irec+1).gt.500) ztemp(nrec-irec+1) = 500
        ztemp_3(irec)=ztemp(nrec-irec+1) !MOD VKP
      enddo
          
      final_age = 0.0
      fmeanp = 0.0
      dfmeanp = 0.0
      oldest_age  = 0.0

      if (ageflag(3).eq.1.and.aftmodel.ge.1.and.kinFTL_AFT(i).ne.-9999) then
        call ketch_main(nrec,ztime_3,ztemp_3,aftmodel,rhoST,Init_FTL_Model,kinFTL_AFT(i),&
            kinFTLID,Init_FTL_value,final_age,oldest_age,fmeanp,dfmeanp,fdist)

        if (final_age.le.1e-4) then
          final_age = 1e-4
        endif

        age3(i)=real(final_age,4)
        fmean(i) = fmeanp
        dfmean(i) = dfmeanp
      endif

      if (ageflag(9).eq.2) then
        do k = 1, 20
          !ftdist(k,i) = sum(fdist(1+(199*(k-1))/20:1+(199*k)/20))/(200/20)
          ftdist(k,i) = sum(fdist(1+10*(k-1):10*k))/(200/20)
        enddo
      else
        ftdist(:,i)=0.
      endif

    else ! Deprecited..
      ! Maxime
      if (ageflag(3).eq.1.and.aftmodel.gt.2.and.kinFTL_AFT(i).ne.-9999) then
        allocate (rho_red(nrec))
        allocate (ztime_double(nrec),ztemp_double(nrec))
        ztime_double = ztime
        ztemp_double = ztemp
              
        call FissionTrackAge(ztime_double,ztemp_double, nrec,rho_red,age3(i),1,rhoST,kinFTL_AFT(i),fdist)
        deallocate (rho_red)
        deallocate (ztime_double, ztemp_double)
      endif             

      if (ageflag(3).eq.1.and.aftmodel.eq.0) call Mad_Trax (ztime,ztemp,nrec,1,2, &
                                            age3(i),ftld,ftldmean,ftldsd)
      if (ageflag(9).eq.1) then
        ftdist(:,i)=0.d0
        ftdist(1:17,i)=ftld(1:17)/100.d0
      else
        ftdist(:,i)=0.
      endif
    endif
    ! ------------------- End Apatite Fission track ages ------------------
        
    ! ------------------ Zircon Fission track ages ------------------------
    if (ageflag(4).eq.1) call Mad_Zirc (ztime,ztemp,nrec,0,2, &
                                            age4(i),ftld,ftldmean,ftldsd)
    ! ------------------------ End zircon FT ages -------------------------
        
    ! ------------------------Apatite helium ages -------------------------
          
    allocate (He_ratio(ndur),He_ratioEjec(ndur),TotHe3_released(ndur),StepAge(ndur))
    allocate (He43obs(ndur),SF3Heobs(ndur))

    if (ageflag(1).eq.1.and.ASize(i).gt.0.and.RDmodel.gt.1) then ! Maxime
      
      allocate (ztime_double(nrec),ztemp_double(nrec),dummy(nrec))
              
      ztime_double = ztime
      ztemp_double = ztemp
      do k=1,ndur
       SF3Heobs(k) = dble(k)/dble(ndur)  ! Initiate cumulative fraction of He3 array
      enddo
      
      if (kinFTL_AHe(i).ne.-9999) then
        ! Compute rmr0
        call Compute_rmr0(kinFTL_AHe(i), kinFTLID)
        ! print *, 'kinFTL = ', kinFTL(i)
        ! Helium diffusion
        
        call Hediff(ztime_double, ztemp_double, nrec, duration, temp_heat, ndur,&
                age1(i), ASize(i), Alpha_Ejec_Flag, RDmodel,Uppm(i),&
                Thppm(i),D0_spec,Ea_spec,kinFTL_AHe(i),ageTime_a,He43_flag,He_ratio,&
                TotHe3_released,He_ratioEjec,He43obs,SF3Heobs,0,data_ADAM_array,1,0)
      endif
      deallocate (ztime_double, ztemp_double,dummy)

      if (age1(i).le.1e-4) then
        age1(i) = 1e-4
      endif
      ! Edge Age
      S = 20 ! Average stopping distance (um)
      FT = 1 - (3*S/(4*ASize(i))) + S**3/(16*S**3) ! alpha-ejection correction factor
      Fedge = - 0.25*((S-2*ASize(i))/ASize(i)) ! Fraction of alphas retained at the edge
      RsRbAlpha =  Fedge/FT * 1 ! we assume no diffusivity fractionation between 4He and 3He
      ! EdgeAge(i) = He_ratio(2)/RsRbAlpha * age1(i)
      EdgeAge(i) = He_ratio(2)/He_ratioEjec(2) * age1(i)
      !##### Calculate Step Age profile ####
      StepAge = He_ratio(:)/He_ratioEjec(:) * age1(i)
              
    elseif (ageflag(1).eq.1.and.RDmodel.lt.2) then
      allocate (dummy(nsurf))
      allocate (released(ndur),agereleased(ndur))
      allocate (ztime_double(nrec),ztemp_double(nrec))
      ztime_double = ztime
      ztemp_double = ztemp
      
      call Mad_He (ztime,ztemp,nrec,age1(i),1,ASize(i),D0_spec,Ea_spec)
  
      if (He43_flag.eq.1) then
          call Diffusion1D (ztime_double,ztemp_double,nrec,temp_heat,duration,ndur, &
                            dummy(i),released,agereleased,ASize(i))
      endif

      if (age1(i).le.1e-4) then
                age1(i) = 1e-4
      endif
      !###### Calculate average 4He/3He #####
      EdgeAge(i) = 0.
      ns43 = 0
      ! Edge Age
      S = 20 ! Average stopping distance (um)
      FT = 1 - (3*S/(4*ASize(i)))  + S**3/(16*S**3)! alpha-ejection correction factor
      Fedge = - 0.25*((S-2*ASize(i))/ASize(i)) ! Fraction of alphas retained at the edge
      RsRbAlpha =  Fedge/FT * 1 ! we suppose no diffusivity fractionation between 4He and 3He
      EdgeAge(i) = agereleased(2)/RsRbAlpha * age1(i)
              
      deallocate (dummy,released,agereleased,ztime_double,ztemp_double)
    endif   
          
    deallocate (He_ratio,He_ratioEjec,TotHe3_released,StepAge)
    deallocate (He43obs,SF3Heobs)
          
    ! --------------------End Apatite helium ages -------------------------
        
    ! ----------------------- Zircon helium ages --------------------------
    if (ageflag(2).eq.1.and.ZSize(i).gt.0) then ! Maxime
      allocate (ztime_double(nrec),ztemp_double(nrec))
            allocate (He43obs(ndur),SF3Heobs(ndur))
      ztime_double = ztime
      ztemp_double = ztemp
      if (RDmodelz.lt.6) then
        call Mad_He (ztime,ztemp,nrec,age2(i),2,ZSize(i),D0z,Eaz)
      else 
                do k=1,ndur
                    SF3Heobs(k) = k/ndur ! Initiate cumulative fraction of He3 array
                enddo
        call Hediff(ztime_double, ztemp_double, nrec, duration, temp_heat, ndur,&
              age2(i), ZSize(i), Alpha_Ejec_Flag, RDmodelz,ZUppm(i),&
              ZThppm(i),D0z,Eaz,0.d0,ageTime_z,0,He_ratio,&
              TotHe3_released,He_ratioEjec,He43obs,SF3Heobs,0,data_ADAM_array,1,1)
      endif

            if (age2(i).le.1e-4) then
                age2(i) = 1e-4
      endif


      deallocate (ztime_double, ztemp_double)
      deallocate (He43obs,SF3Heobs)
    endif   
    ! --------------------End Zircon helium ages -------------------------- 

    if (ageflag(5).eq.1) call Mad_He (ztime,ztemp,nrec,age5(i),3,grainsize(i),D0k,Eak)
        if (age5(i).le.1e-4) age5(i) = 1e-4

    if (ageflag(6).eq.1) call Mad_He (ztime,ztemp,nrec,age6(i),4,grainsize(i),D0b,Eab)
        if (age6(i).le.1e-4) age6(i) = 1e-4

    if (ageflag(7).eq.1) call Mad_He (ztime,ztemp,nrec,age7(i),5,grainsize(i),D0m,Eam)
        if (age7(i).le.1e-4) age7(i) = 1e-4

    if (ageflag(8).eq.1) call Mad_He (ztime,ztemp,nrec,age8(i),6,grainsize(i),D0h,Eah)
        if (age8(i).le.1e-4) age8(i) = 1e-4

        if (timenow.lt.0) timenow = 0.0
    age1(i)=age1(i)+timenow
    age2(i)=age2(i)+timenow
    age3(i)=age3(i)+timenow
    age4(i)=age4(i)+timenow
    age5(i)=age5(i)+timenow
    age6(i)=age6(i)+timenow
    age7(i)=age7(i)+timenow
    age8(i)=age8(i)+timenow
    
  enddo
  ! print *, 'Timenow: ', timenow
  if (iproc.eq.0.and.nd.eq.0) print*,''

  deallocate (ztime,ztime_3,ztemp,ztemp_3,depth)

  if (saveTt.and.nd.eq.0) then
      write (line,'(a,1023(",",a))') "Time",(trim(sampleID(k)),k=1,nsurf)
    write (82,'(a)') line(:len(trim(line))-1)           
    do ij=1,nrec
            write (line,'(g12.6,1023(",",g12.6))') (mem(ij,i),i=0,nsurf)
        write (82,'(a)') line(:len(trim(line))-1)                                    
    enddo
  endif

  if (saveTt) deallocate (mem)

  return
end
