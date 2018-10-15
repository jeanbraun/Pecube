program Cascade2Pecube

implicit none

real time,latmean,lapse_rate,tmsl,dum
character a5*5,c5*5,fnme*256
integer nt,i,k,nnode,nnodep,istep,nc5,nfnme
integer,dimension(:,:),allocatable::icon
real,dimension(:),allocatable::x,y,h,e

  do i=1,256
  fnme(i:i)=' '
  enddo
print*,'Enter directory name (in data directory) where cascade output is to be saved >'
read (*,'(a)') fnme
  do i=1,256
  if (fnme(i:i).ne.' ') nfnme=i
  enddo

open (7,file='cascade_output/connectivity',status='old')
read (7,*) a5,time,nt
print*,nt,'triangles'
allocate (icon(3,nt))
  do i=1,nt
  read (7,*) k,icon(:,i)
  enddo
close (7)

open (7,file='cascade_output/geometry',status='old')
read (7,*) a5,time,nnode
print*,nnode,'nodes'
allocate (x(nnode),y(nnode))
  do i=1,nnode
  read (7,*) k,x(i),y(i)
  enddo
close (7)

latmean=sum(y)/nnode/111111.
open (8,file='data/'//fnme(1:nfnme)//'/geometry',status='unknown') 
  do i=1,nnode
  write (8,*) x(i)/111111./cos(latmean/180.*3.141592654),y(i)/111111.
  enddo
  do i=1,nt
  write (8,*) icon(:,i)
  enddo
close (8)

allocate (h(nnode),e(nnode))
istep=-1

print*,'Enter temperature at mean sea level (in °C) >'
read*,tmsl
print*,'Enter lapse rate (in°C/km) >'
read*,lapse_rate
e=0.

open (7,file='cascade_output/topography',status='old')
open (9,file='cascade_output/uplift',status='old')
2 read (7,*,end=1) a5,time,nnodep
if (istep.ne.-1) read (9,*) a5,time,nnodep
if (nnodep.ne.nnode) stop 'this version of Pecube cannot handle cascade problems with time-varying nodal geometry...'
  do i=1,nnode
  read (7,*) k,h(i)
  enddo
  if (istep.ne.-1) then
    do i=1,nnode
    read (9,*) k,e(i)
    enddo
  endif
istep=istep+1
print*,'doing step',istep,'at time',time
  if (istep.lt.10) then 
  write (c5(1:1),'(i1)') istep
  nc5=1
  elseif (istep.lt.100) then
  write (c5(1:2),'(i2)') istep
  nc5=2
  elseif (istep.lt.1000) then
  write (c5(1:3),'(i3)') istep
  nc5=3
  elseif (istep.lt.10000) then
  write (c5(1:4),'(i4)') istep
  nc5=4
  else
  write (c5(1:5),'(i5)') istep
  nc5=5
  endif
open (8,file='data/'//fnme(1:nfnme)//'/topo'//c5(1:nc5),status='unknown')
  do i=1,nnode
  write (8,*) h(i)
  enddo
close (8)
open (8,file='data/'//fnme(1:nfnme)//'/uplift'//c5(1:nc5),status='unknown')
  do i=1,nnode
  write (8,*) e(i)*1.e3
  enddo
close (8)
open (8,file='data/'//fnme(1:nfnme)//'/temp'//c5(1:nc5),status='unknown')
  do i=1,nnode
  write (8,*) tmsl-h(i)*lapse_rate/1.e3
  enddo
close (8)
goto 2
1 close (7)
close (9)
deallocate (h,e)

end program Cascade2Pecube
