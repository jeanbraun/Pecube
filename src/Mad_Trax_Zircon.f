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

      subroutine Mad_Zirc (time_i,temp_i,n,out_flag,param_flag,
     &                     fta,ftld,ftldmean,ftldsd)

c subroutine Mad_Trax to calculate zircon fission track age from a given thermal history
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
c NB integer out_flag=1 not implemented yet for Zircon !

c integer param_flag :   =1 uses parameters for alpha-damaged zircon
c                        =2 uses parameters for zero-damage zircon
c                        (parameters from Rahn et al. 2004)
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
c Peter van der Beek in December 1995. The algorithm is explained in
c Peter's PhD thesis and is based on the work by Lutz and Omar (1991)
c This adaptation of the program for zircon fission-track annealing was
c written by Peter in August/September 2006 and is based on algorithms
c given by Galbraith & Laslett (1997), Tagami et al. (1998) and
c Rahn et al. (2004)
c
c References:
c
c van der Beek, P., 1995.Tectonic evolution of continental rifts, PhD Thesis,
c       Faculty of Earth Sicences, Free University, Amsterdam.
c
c Lutz, T.M. and Omar, G.I., 1991. An inverse method of modeling thermal
c       histories from apatite fission-track data. EPSL, 104, 181-195.
c
c Galbraith, R. F., and G. M. Laslett (1997), Statistical modelling of thermal
c       annealing of fission tracks in zircon, Chemical Geology, 140, 123-135.
c
c Tagami, T., et al. (1998), Revised annealing kinetics of fission tracks in 
c        zircon and geological implications, in Advances in Fission-Track Geochronology, 
c        edited by P. Van den haute and F. De Corte, pp. 99-112, Kluwer Academic 
c        Publishers, Dordrecht, Netherlands.
c
c Rahn, M. K., et al. (2004), A zero-damage model for fission track annealing in zircon,
c        American Mineralogist, 89, 473-484.


      real*4   temp_i(n),time_i(n),a,b,c
      real*4   r(1000),prob(101),ftld(*)
      integer  out_flag,param_flag

      if (param_flag.eq.1) then

c alpha-damaged zircon

      a=-10.77
      b=2.599E-4
      c=1.026E-2

      else

c zero-damage zircon

      a=-11.57
      b=2.755E-4
      c=1.075E-2

      endif

c unannealed fission track length

      xind=10.8

c mean length of spontaneous tracks in standards

      xfct=10.8

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

      tempp=temperaturezr(temp_i,time_i,0.,n)+273.
      rp=.5

c beginning of time stepping

        do i=1,nstep

        time=float(i)*time_interval

c calculate temperature by linear interpolation

        temp=temperaturezr(temp_i,time_i,time,n)+273.

c calculate mean temperature over the time step

        tempm=(temp+tempp)/2.

c calculate the "equivalent time", teq

        if (i.eq.1)then
          teq=0.
        else
          teq=exp((gzr(rp)-a-(c*tempm))/(b*tempm))
        endif
        

c check if we are not getting too close to r=0
c in which case r remains 0 for all following time steps

c        if (param_flag.lt.4) then
c           rcrit=(expos(1./b,a)-a*c0-1.)/a/c1*(1./tempm-c3)-c2/c1
c        else
c           rcrit=(expos(1./b,a)-a*c0-1.)/a/c1*(alog(1./tempm)-c3)-c2/c1
c        endif

c          if (alog(teq+deltat).ge.rcrit) then
c            do j=i,nstep
c            r(j)=0.
c            enddo
c          nstep=i    
c          goto 90

c          endif

c otherwise calculate reduction in length, r, over the time step, dt

        dt=teq+deltat
        
        gr=a+((b*tempm)*alog(dt))+(c*tempm)
        
        r(i)=xinvzr(gr)
c stop calculation if r<0.4
        if (r(i).le.0.4) then
          nstep=i
          goto 90
        endif

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

          if (r(i).le.0.4) then
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
          dfr=xkzr(x)/h
          fr=fr+dfr

          enddo

        prob(j)=fr/nstep
        sumprob=sumprob+prob(j)

        enddo

c now let's rescale the track length distribution, its mean and standard
c deviation

      ftld(11)=100.
      imin=1

        do l=1,10

        imax=int(l*100./xind)
        ftld(l)=0.

          do i=imin,imax
          ftld(l)=ftld(l)+prob(i)
          enddo

        ftld(l)=(ftld(l)*100./sumprob)
        ftld(11)=ftld(11)-ftld(l)

        imin=imax+1

        enddo

      sumftld=0.

        do l=1,11
        sumftld=sumftld+ftld(l)*(float(l)-0.5)
        enddo

      ftldmean=sumftld/100.

      devftld=0.

        do l=1,11
        devftld=devftld+ftld(l)*(float(l)-0.5-ftldmean)**2
        enddo

      ftldsd=sqrt(devftld/100.)

      return
      end

c---
      function gzr(r)

      gzr=alog(1-r)

      return
      end

c---
      function xinvzr (gr)

      xinvzr=1-exp(gr)

      return
      end

c---
      function temperaturezr (temp,time,t,n)

      real*4   temp(n),time(n)

      temperaturezr=temp(1)

        do i=1,n-1
          if ((t-time(i))*(t-time(i+1)).le.0.) then
          rat=(t-time(i))/(time(i+1)-time(i))
          temperaturezr=temp(i)+rat*(temp(i+1)-temp(i))
          return
          endif
        enddo

      temperaturezr=temp(n)

      return
      end


c---
      function xkzr (x)

      xkzr=0.39894228*exp(-(x**2)/2.)

      return
      end
