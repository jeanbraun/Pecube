      subroutine Mad_He (time,temperature,ntime,Apatite_age,iflag,grainsize)

! iflag=1 He in Apatite
! iflag=2 He in Zircon
! iflag=3 Ar in K-feldspar
! iflag=4 Ar in Biotite
! iflag=5 Ar in Muscovite
! iflag=6 Ar in Hornblende

      real time(ntime),temperature(ntime),grainsize
      real,dimension(:),allocatable::age,diag,sup,inf,f

!      D0a2=10.**7.7*3600.*24.*365.25e6
!      Ea=36.2e3*4.184

      if (iflag.eq.1) then
! He in Ap
      D0a2=5.e5*3600.*24.*365.25e6*(100./grainsize)**2
      Ea=138.e3
      elseif (iflag.eq.2) then
! He in Zi
      D0a2=4.6e3*3600.*24.*365.25e6
      Ea=168.e3
      elseif (iflag.eq.3) then
! Ar in Ks
      D0a2=5.6*3600.*24.*365.25e6
      Ea=120.e3
      elseif (iflag.eq.4) then
! Ar in Bi 500mic
      D0a2=160.*3600.*24.*365.25e6
      Ea=211.e3
      elseif (iflag.eq.5) then
! Ar in Mu 500mic
      D0a2=13.2*3600.*24.*365.25e6
      Ea=183.e3
      elseif (iflag.eq.6) then
! Ar in Ho 500mic
      D0a2=24.*3600.*24.*365.25e6
      Ea=276.e3
      else
      stop 'option not supported in Mad_He'
      endif

      R=8.314
      dt0=0.1

      n=100

      allocate (age(n),diag(n),sup(n),inf(n),f(n))

      age=0.

      do itime=1,ntime-1

      dt=dt0
1     nstep=max(1,int((time(itime)-time(itime+1)+tiny(dt))/dt))
      dt=(time(itime)-time(itime+1))/nstep
      alpha=0.5
      dr=1./(n-1)

! beta determines the geometry: 2 = spherical
!                               1 = cylindrical
      beta=2.

      temps=temperature(itime)
      tempf=temperature(itime+1)

      Da2now=D0a2*exp(-Ea/R/(temps+273.))

        do istep=1,nstep

        fstep=float(istep)/(nstep)
        temp=temps+(tempf-temps)*fstep
        Da2then=Da2now
        Da2now=D0a2*exp(-Ea/R/(temp+273.))
        f1=alpha*dt*Da2now/dr**2
        f2=(1.-alpha)*dt*Da2then/dr**2

          do i=2,n-1
          diag(i)=1.+2.*f1
          sup(i)=-f1*(1.+beta/(i-1)/2.)
          inf(i)=-f1*(1.-beta/(i-1)/2.)
          f(i)=age(i)+f2* &
               ((age(i+1)-2.*age(i)+age(i-1)) + &
                beta*(age(i+1)-age(i-1))/(i-1)/2.)+dt
          enddo

        diag(1)=1.
        sup(1)=-1.
        f(1)=0.
        diag(n)=1.
        inf(n)=0.
        f(n)=0.

        call tridag (inf,diag,sup,f,age,n)

        agei=0.
          do i=1,n
          fact=1.
          if (i.eq.1.or.i.eq.n) fact=0.5
          agei=agei+age(i)*fact*dr**3*(i-1)*(i-1)
          enddo
        agei=3.*agei

        enddo

      enddo

      Apatite_age=agei

      return
      end
