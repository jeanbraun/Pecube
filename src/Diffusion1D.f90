      subroutine Diffusion1D (time,temperature,ntime,tempheat,duration,nheat, &
                              agei,released,agereleased,size)

! subroutine to compute He Diffusion profiles in Apatite and apparent age
!
! in input:
! time-temperature: two arrays of length ntime defining the geological temperature history of the sample
! durationheat and tempheat: two arrays of length nheat defining duration and temperature of nheat heating steps in second
!
! in output
! agei: sample age in Myr
! relesased and agereleased: two arrays of length nheat containing the age spectrum, i.e. R/Rbulk vs fraction parent released

      implicit none

      double precision time(ntime),temperature(ntime)
      double precision tempheat(nheat),released(nheat),agereleased(nheat),duration(nheat)
      double precision, dimension(:),allocatable::diag,sup,inf,f,age,prod,parent,g,parentreleased
      integer ntime,itime,istep,nstep,i,j,n,nheat,iheat
      double precision D0a2,Ea,R,dt0,dt,alpha,dr,beta,temps,tempf,stopping
      double precision Da2now,fstep,Da2then,temp,f1,f2,agei,fact,rat,prodi,fraci,fracp,parent0
      double precision fracin,fracpn,bulk,size

! first set the diffusion parameters

!      D0a2=5.e5*3600.*24.*365.25e6
      D0a2=0.0032/(size/1.d6)**2*3600.*24.*365.25e6
      Ea=138.d3
      R=8.314d0

! stopping distance is expressed in radius units

      stopping=0.2d0

! time step for geological cooling/history

      dt0=0.1d0

! number of points used to discretize the radius of the grain

      n=100
      dr=1.d0/(n-1)

! alpha defines the time integration scheme

      alpha=0.5d0

! allocate memory

      allocate (age(n),diag(n),sup(n),inf(n),f(n),prod(n),parent(n),g(n))
      allocate (parentreleased(nheat))

      prod=1.d0

! includes stopping distance in production term

      do i=1,n
      fact=1.d0-float(i-1)/(n-1)
      if (fact.lt.stopping) prod(i)=prod(i)*fact/stopping
      enddo

! normalizes production by integrated production

      prodi=0.d0
        do i=1,n
        fact=1.d0
        if (i.eq.1.or.i.eq.n) fact=0.5d0
        prodi=prodi+prod(i)*fact*dr**3*(i-1)*(i-1)
        enddo
      prodi=3.d0*prodi

      prod=prod/prodi

! PART 1
! starts calculating age by solving the diffusion equation through geological time
! age will contain the profile in age across the grain

      age=0.d0

      do itime=1,ntime-1

      dt=dt0
      nstep=max(1,int((time(itime)-time(itime+1)+tiny(dt))/dt))
      dt=(time(itime)-time(itime+1))/nstep

! beta determines the geometry: 2 = spherical
!                               1 = cylindrical
      beta=2.d0

      temps=temperature(itime)
      tempf=temperature(itime+1)

      Da2now=D0a2*exp(-Ea/R/(temps+273.15d0))

        do istep=1,nstep

        fstep=float(istep)/(nstep)
        temp=temps+(tempf-temps)*fstep
        Da2then=Da2now
        Da2now=D0a2*exp(-Ea/R/(temp+273.15d0))
        f1=alpha*dt*Da2now/dr**2
        f2=(1.d0-alpha)*dt*Da2then/dr**2

          do i=2,n-1
          diag(i)=1.d0+2.d0*f1
          sup(i)=-f1*(1.d0+beta/(i-1)/2.d0)
          inf(i)=-f1*(1.d0-beta/(i-1)/2.d0)
          f(i)=age(i)+f2* &
               ((age(i+1)-2.d0*age(i)+age(i-1)) + &
                beta*(age(i+1)-age(i-1))/(i-1)/2.d0)+dt*prod(i)
          enddo

        diag(1)=1.d0
        sup(1)=-1.d0
        f(1)=0.d0
        diag(n)=1.d0
        inf(n)=0.d0
        f(n)=0.d0

        call tridiagonal (inf,diag,sup,f,age,n)

        enddo

      enddo

! integrate daughter product to compute age

      agei=0.d0
        do i=1,n
        fact=1.d0
        if (i.eq.1.or.i.eq.n) fact=0.5d0
        agei=agei+age(i)*fact*dr**3*(i-1)*(i-1)
        enddo
      agei=3.d0*agei

! PART 2
! simulates step heating experiment by solving the diffusion equation
! using a minute time stepping

      dt0=60.d0

! released is released 4He
! age is 4He in sample
! parent is 3He in sample

      released=0.d0
      parent=1.d0

! goes through heating steps given by durationheat

      do iheat=1,nheat

      dt=dt0
      nstep=max(1,int((duration(iheat)*3600.d0+tiny(dt))/dt))
      dt=duration(iheat)/nstep/(24.d0*365.25d6)

! beta determines the geometry: 2 = spherical
!                               1 = cylindrical
      beta=2.d0

      temp=tempheat(iheat)

      Da2now=D0a2*exp(-Ea/R/(temp+273.15d0))

! fraci is integrated 4He in sample at start of time step

      fraci=0.d0
        do i=1,n
        fact=1.d0
        if (i.eq.1.or.i.eq.n) fact=0.5d0
        fraci=fraci+age(i)*fact*dr**3*(i-1)*(i-1)
        enddo
      fraci=3.d0*fraci

! fracp is integrated 3He in sample at start of time step

      fracp=0.d0
        do i=1,n
        fact=1.d0
        if (i.eq.1.or.i.eq.n) fact=0.5d0
        fracp=fracp+parent(i)*fact*dr**3*(i-1)*(i-1)
        enddo
      fracp=3.d0*fracp

! time steping by 60 seconds time steps

        do istep=1,nstep

        fstep=float(istep)/(nstep)
        f1=alpha*dt*Da2now/dr**2
        f2=(1.d0-alpha)*dt*Da2now/dr**2

          do i=2,n-1
          diag(i)=1.d0+2.d0*f1
          sup(i)=-f1*(1.d0+beta/(i-1)/2.d0)
          inf(i)=-f1*(1.d0-beta/(i-1)/2.d0)
          f(i)=age(i)+f2* &
               ((age(i+1)-2.d0*age(i)+age(i-1)) + &
                beta*(age(i+1)-age(i-1))/(i-1)/2.d0)
          g(i)=parent(i)+f2* &
               ((parent(i+1)-2.d0*parent(i)+parent(i-1)) + &
                beta*(parent(i+1)-parent(i-1))/(i-1)/2.d0)
          enddo

        diag(1)=1.d0
        sup(1)=-1.d0
        f(1)=0.d0
        g(1)=0.d0
        diag(n)=1.d0
        inf(n)=0.d0
        f(n)=0.d0
        g(n)=0.d0

        call tridiagonal (inf,diag,sup,f,age,n)
        call tridiagonal (inf,diag,sup,g,parent,n)

        enddo

! fracin is integrated 4He in sample at end of time step

      fracin=0.d0
        do i=1,n
        fact=1.d0
        if (i.eq.1.or.i.eq.n) fact=0.5d0
        fracin=fracin+age(i)*fact*dr**3*(i-1)*(i-1)
        enddo
      fracin=3.d0*fracin

! fracpn is integrated 3He in sample at end of time step

      fracpn=0.d0
        do i=1,n
        fact=1.d0
        if (i.eq.1.or.i.eq.n) fact=0.5d0
        fracpn=fracpn+parent(i)*fact*dr**3*(i-1)*(i-1)
        enddo
      fracpn=3.d0*fracpn

! released is total amount of 3He released
! age released is apparent age of gas released during time step (delta 4He/delta 3He)
! parentreleased is 3He released over time step

      released(iheat)=1.d0-fracpn
      agereleased(iheat)=(fraci-fracin)/(fracp-fracpn)
      parentreleased(iheat)=fracp-fracpn

      enddo

      bulk=sum(agereleased*parentreleased)/sum(parentreleased)
      agereleased=agereleased/bulk

      return
      end

!!!!!

! utility routine

      SUBROUTINE tridiagonal(a,b,c,r,u,n)
      INTEGER n
      double precision a(n),b(n),c(n),r(n),u(n)
      INTEGER j
      double precision bet
	double precision,dimension(:),allocatable::gam
	allocate (gam(n))
      if(b(1).eq.0.) stop 'in tridag'
      bet=b(1)
      u(1)=r(1)/bet
      do 11 j=2,n
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j)*gam(j)
        if(bet.eq.0.) then
        print*,'tridag failed'
        stop
        endif
        u(j)=(r(j)-a(j)*u(j-1))/bet
11    continue
      do 12 j=n-1,1,-1
        u(j)=u(j)-gam(j+1)*u(j+1)
12    continue
      deallocate (gam)
      return
      END
