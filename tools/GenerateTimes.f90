program GenerateTimes

character c1*1,c2*2,c3*3,c4*4
integer :: ntime,i
double precision :: tend

print*,'Enter number of times'
read*,ntime

print*,'Enter starting time (in Myr)'
read*,tend

do i=1,ntime
if (i.lt.10) then
write (c1,'(i1)') i
write (*,'(a,f12.6)') 'time_topo'//c1//' = ',tend*float(ntime-i+1)/ntime
elseif (i.lt.100) then
write (c2,'(i2)') i
write (*,'(a,f12.6)') 'time_topo'//c2//' = ',tend*float(ntime-i+1)/ntime
elseif (i.lt.1000) then
write (c3,'(i3)') i
write (*,'(a,f12.6)') 'time_topo'//c3//' = ',tend*float(ntime-i+1)/ntime
elseif (i.lt.10000) then
write (c4,'(i4)') i
write (*,'(a,f12.6)') 'time_topo'//c4//' = ',tend*float(ntime-i+1)/ntime
endif
enddo

end program GenerateTimes
