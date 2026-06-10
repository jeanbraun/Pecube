cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                           c
c  MMM        MMM                TTTTTTTTTT                 c 
c  MMMM      MMMM                TTTTTTTTTT                 c 
c  MMMMM    MMMMM                    TT                     c
c  MM  MM  MM  MM                    TT                     c
c  MM   MMMM   MM             d      TT                     c
c  MM    MM    MM             d      TT                     c
c  MM          MM  aaa a  ddd d      TT r rrr  aaa a x   x  c
c  MM          MM a   aa d   dd      TT rrr   a   aa  x x   c
c  MM          MM a   aa d   dd      TT r     a   aa   x    c
c  MM          MM a   aa d   dd      TT r     a   aa  x x   c
c  MM          MM  aaa a  ddd d      TT r      aaa a x   x  c
c                                                           c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c

      subroutine Mad_Trax (time_i,temp_i,n,out_flag,param_flag,
     &                     fta,ftld,ftldmean,ftldsd)

c subroutine Mad_Trax to calculate fission track age, track length
c distribution and statistics from a given thermal history
c
c in input:
c real*4  time_i(n)  :   the time values (in Myr) in descending order
c                        at which the thermal history is given
c                        (ex: 100,50,20,10,0); the last value
c                        should always be 0; the first value
c                        should be smaller than 1000.
c real*4  temp_i(n)  :   the thermal history in degree Celsius
c integer n          :   the number of time-temperature pairs used
c                        to describe the temperature history
c integer out_flag   :   =0 only calculate fission track age
c                        =1 also calculate track length distribution
c                           and statistics
c integer param_flag :   =1 uses Laslett et al, 1987 parameters
c                        =2 uses Crowley et al., 1991 Durango parameters
c                        =3 uses Crowley et al., 1991 F-apatite parameters
c
c in output:
c real*4  fta        :   fission track age in Myr
c real*4  ftld(17)   :   normalised track length distribution where
c                        ftld(k) is the percentage of track with
c                        length between k-0.5 and k+0.5 microns
c real*4  ftldmean   :   mean track length in microns
c real*4  ftldsd     :   track length standard deviation in microns
c
c This subroutine is based on the subroutine "ftmod.pas" provided by
c Peter vanderBeek in December 1995. The algorithm is explained in
c Peter's PhD thesis and is based on the work by Lutz and Omar (1991)
c
c References:
c
c VanderBeek, P., 1995.Tectonic evolution of continental rifts, PhD Thesis,
c       Faculty of Earth Sicences, Free University, Amsterdam.
c
c Lutz, T.M. and Omar, G.I., 1991. An inverse method of modeling thermal
c       histories from apatite fission-track data. EPSL, 104, 181-195.
c
c Laslett, G.M., Green, P.F., Duddy, I.R. and Gleadow, A.J.W., 1987. Thermal
c       annealing of fission tracks in apatite 2. A quantitative analysis.
c       Chem. Geol. (Isot. Geosci. Sect.) 65, 1-13.
c
c Crowley, K.D., Cameron, M. and Schaefer, R.L., 1991. Experimental studies
c       of annealing of etched fission tracks in fluorapatite. Geochim.
c       Cosmochim. Acta, 55, 1449-1465.
c


      real*4   temp_i(n),time_i(n)
c MOD Jean 31/10/08
c r(100) modified into r(1000)
      real*4   r(1000),prob(101),ftld(*)
c end MOD
      integer  out_flag,param_flag

      if (param_flag.eq.1) then

c from Laslett et al., 1987

      a=0.35
      b=2.7
      c0=-4.87
      c1=0.000168
      c2=0.00472416
      c3=0.

      elseif (param_flag.eq.2) then

c from Crowley et al., 1991 (Durango)

      a=0.49
      b=3.0
      c0=-3.202
      c1=0.0000937
      c2=0.001839
      c3=0.00042

      else

c from Crowley et al., 1991 (F-apatite)

      a=0.76
      b=4.16
      c0=-1.508
      c1=0.00002076
      c2=0.0002143
      c3=0.0009967

      endif

c unannealed fission track length

      xind=16.

c mean length of spontaneous tracks in standards

      xfct=14.5

c calculate the number of time steps assuming 1My time step length
c if run time > 100 My, else take 100 time steps

      nstep=int(time_i(1))
c     if (nstep.gt.8000) stop 'Your final time is greater than '//
c    &                        'the age of the Universe...'
c     if (nstep.gt.4500) stop 'Your final time is greater than '//
c    &                        'the age of the Earth...'
c     if (nstep.gt.1000) stop 'Fission track does not work very well '//
c    &                        'for time spans greater than 1Byr...'
      if (nstep.gt.1000) stop 
      time_interval=1.
      if (nstep.lt.100) then
        nstep=100
        time_interval=time_i(1)/100.
      endif
      deltat=time_interval*1.e6*365.24*24.*3600.

c calculate final temperature

      tempp=temperature(temp_i,time_i,0.,n)+273.
      rp=.5

c begining of time stepping

        do i=1,nstep

        time=float(i)*time_interval

c calculate temperature by linear interpolation

        temp=temperature(temp_i,time_i,time,n)+273.

c calculate mean temperature over the time step

        tempm=(temp+tempp)/2.

c calculate the "equivalent time", teq

        teq=exp((-c2/c1)+((g(rp,a,b)-c0)/c1)*(1./tempm-c3))
        if (i.eq.1) teq=0.

c check if we are not getting too close to r=0
c in which case r remains 0 for all following time steps

          if (alog(teq+deltat).ge.
     &     (expos(1./b,a)-a*c0-1.)/a/c1*(1./tempm-c3)-c2/c1) then

            do j=i,nstep
            r(j)=0.
            enddo
          nstep=i    
          goto 90

          endif

c otherwise calculte reduction in length, r, over the time step, dt

        dt=teq+deltat
        gr=c0+((c1*alog(dt)+c2)/((1./tempm)-c3))
        r(i)=xinv(gr,a,b)

c update variables for next time step

        tempp=temp
        rp=r(i)

c       print*,i,time,temp,r(i)

        enddo

90    continue

c all reduction factors for all time steps have been calculated
c now estimate the fission track age by simple summation
c (here it helps to use 1Myr time steps)

      sumdj=0.

        do i=1,nstep

          if (r(i).le.0.35) then
          dj=0.
          elseif (r(i).le.0.66) then
          dj=2.15*r(i)-0.76
          else
          dj=r(i)
          endif

        sumdj=sumdj+dj

        enddo

      fta=(xind/xfct)*sumdj*time_interval

c now (if out_flag.ne.0) let's do some statistics

      if (out_flag.eq.0) return

c first, calculate probability density function using Luts and Omar (1991)
c method and assuming a Gaussian distribution

      sumprob=0.

        do j=1,101

        rr=float(j-1)/100.

          if (rr.le..43) then
          h=2.53
          elseif (rr.le..67) then
          h=5.08-5.93*rr
          else
          h=1.39-.61*rr
          endif

        fr=0.

          do i=1,nstep

          x=(rr-r(i))*xind/h
          dfr=xk(x)/h
          fr=fr+dfr

          enddo

        prob(j)=fr/nstep
        sumprob=sumprob+prob(j)

        enddo

c now let's rescale the track length distribution, its mean and standard
c deviation

      ftld(17)=100.
      imin=1

        do l=1,16

        imax=int(l*100./xind)
        ftld(l)=0.

          do i=imin,imax
          ftld(l)=ftld(l)+prob(i)
          enddo

        ftld(l)=(ftld(l)*100./sumprob)
        ftld(17)=ftld(17)-ftld(l)

        imin=imax+1

        enddo

      sumftld=0.

        do l=1,17
        sumftld=sumftld+ftld(l)*(float(l)-0.5)
        enddo

      ftldmean=sumftld/100.

      devftld=0.

        do l=1,17
        devftld=devftld+ftld(l)*(float(l)-0.5-ftldmean)**2
        enddo

      ftldsd=sqrt(devftld/100.)

      return
      end

c---
      function g(r,a,b)

      g=(expos((1.-expos(r,b))/b,a)-1.)/a

      return
      end

c---
      function xinv (gr,a,b)

      xinv=expos(1.-b*expos(a*gr+1.,1./a),1./b)

      return
      end

c---
      function temperature (temp,time,t,n)

      real*4   temp(n),time(n)

      temperature=temp(1)

        do i=1,n-1
          if ((t-time(i))*(t-time(i+1)).le.0.) then
          rat=(t-time(i))/(time(i+1)-time(i))
          temperature=temp(i)+rat*(temp(i+1)-temp(i))
          return
          endif
        enddo

      temperature=temp(n)

      return
      end

c---
      function expos (x,a)

      expos=exp(a*alog(x))

      return
      end

c---
      function xk (x)

      xk=0.39894228*exp(-(x**2)/2.)

      return
      end
