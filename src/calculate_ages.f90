      subroutine calculate_ages (nrec,nsurf,nz,istep,age1,age2,age3,age4,age5,age6,age7,age8, &
                                 ftdist,tprevious,ftlflag,ageflag,iproc,nd,saveTt,grainsize)
!__________________________________________________________________________
!This file has been changed to use ketcham program. 
!original version is found in calculate_ages_old.f90 in this same folder
!by VKP
!__________________________________________________________________________

      implicit none

      integer nsurf,nrec,nz,istep,i,ij,irec,k,ftlflag,ageflag(9)
      real*4 age1(nsurf),age2(nsurf),age4(nsurf),age5(nsurf),age3(nsurf) 
      real*4 age6(nsurf),age7(nsurf),age8(nsurf)
      real*4,dimension(:),allocatable::depth,ztime,ztemp,ztime_3,ztemp_3    
      real*4 zjunk,timenow,ftld(17),ftldmean,ftldsd,grainsize(nsurf)
      double precision fdist(200),alo,oldest_age,final_age,fmean
      double precision ftdist(20,nsurf),tprevious
      integer aftmodel
      integer iproc,nd
      double precision,dimension(:,:),allocatable::mem
      logical saveTt
      character*4 c(9999)
      character*10000 line

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
      ftdist = 0. ! bug fix by Yuntao May 4, 2020

        do i=1,nsurf
        if (iproc.eq.0.and.nd.eq.0) call screen_counter (i,nsurf,0)
        ij=(i-1)*nz+1
            do irec=1,nrec
            read (100+istep,rec=irec) ztime(irec),(zjunk,k=1,i-1),ztemp(irec), &
                                    (zjunk,k=i+1,nsurf), &
                                    (zjunk,k=1,i-1),depth(irec)
            enddo
        timenow=ztime(nrec)
        ztime=ztime-timenow
          if (saveTt) then
          mem(:,i)=ztemp
          mem(:,0)=ztime
          endif

        if (aftmodel.eq.1) then
        do irec=1,nrec
            ztime_3(irec)=ztime(nrec-irec+1) !MOD VKP
            if (ztemp(nrec-irec+1).gt.500) ztemp(nrec-irec+1)=500
            ztemp_3(irec)=ztemp(nrec-irec+1) !MOD VKP
        enddo
        alo = 16.3
        final_age = 0
        fmean = 0
        oldest_age  = 0
          if (ageflag(3).eq.1) then
          call ketch_main(nrec,ztime_3,ztemp_3,alo,final_age,oldest_age,fmean,fdist)
          age3(i)=real(final_age,4)
          endif
          if (ageflag(9).eq.1) then
            do k = 1, 20
            !ftdist(k,i) = sum(fdist(1+(199*(k-1))/20:1+(199*k)/20))/(200/20)
            ftdist(k,i) = sum(fdist(1+10*(k-1):10*k))/(200/20)
            enddo
          else
          ftdist(:,i)=0.
          endif

        else
        if (ageflag(3).eq.1) call Mad_Trax (ztime,ztemp,nrec,1,2, &
                                            age3(i),ftld,ftldmean,ftldsd)
          if (ageflag(9).eq.1) then
          ftdist(:,i)=0.d0
          ftdist(1:17,i)=ftld(1:17)/100.d0
          else
          ftdist(:,i)=0.
          endif
        endif

        if (ageflag(4).eq.1) call Mad_Zirc (ztime,ztemp,nrec,0,2, &
                                            age4(i),ftld,ftldmean,ftldsd)

        if (ageflag(1).eq.1) call Mad_He (ztime,ztemp,nrec,age1(i),1,grainsize(i))
        if (ageflag(2).eq.1) call Mad_He (ztime,ztemp,nrec,age2(i),2,grainsize(i))
        if (ageflag(5).eq.1) call Mad_He (ztime,ztemp,nrec,age5(i),3,grainsize(i))
        if (ageflag(6).eq.1) call Mad_He (ztime,ztemp,nrec,age6(i),4,grainsize(i))
        if (ageflag(7).eq.1) call Mad_He (ztime,ztemp,nrec,age7(i),5,grainsize(i))
        if (ageflag(8).eq.1) call Mad_He (ztime,ztemp,nrec,age8(i),6,grainsize(i))
        age1(i)=age1(i)+timenow
        age2(i)=age2(i)+timenow
        age3(i)=age3(i)+timenow
        age4(i)=age4(i)+timenow
        age5(i)=age5(i)+timenow
        age6(i)=age6(i)+timenow
        age7(i)=age7(i)+timenow
        age8(i)=age8(i)+timenow
        enddo

      if (iproc.eq.0.and.nd.eq.0) print*,''

      deallocate (ztime,ztime_3,ztemp,ztemp_3,depth)

      if (saveTt.and.nd.eq.0) then
      write (line,'(a,1023(",",a))') "Time",("point"//c(k),k=1,nsurf)
      write (82,'(a)') line(:len(trim(line))-1)
        do ij=1,nrec
        write (line,'(g12.6,1023(",",g12.6))') (mem(ij,i),i=0,nsurf)
        write (82,'(a)') line(:len(trim(line))-1)
        enddo
      endif

      return
      end
