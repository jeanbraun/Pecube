!------------------------------------------------------------------------------|

subroutine iscanfile (fnme, text, res, res_desc, ires, vocal, nd, range, par)

!------------------------------------------------------------------------------|
! subroutine to read the value of an integer parameter whose name is stored in 
! text  from file fnme the result is stored in res and the flag ires is set 
! to 1 if all went well. if vocal is set yo 1 the routine echoes its find to
! standard output

character*(*) fnme,text,res_desc
integer,intent(out)::res,ires
integer,intent(in)::vocal
integer,intent(inout)::nd
real*4,intent(inout)::range(2,1024),par(1024)

character line*256,record*256,cres0*256
integer ios

if (vocal.eq.4) then
write (*,*) text//' is '//trim(res_desc)
return
endif

ires=0
write (cres0,*) res
cres0=trim(adjustl(cres0))

open (7,file=fnme(1:len_trim(fnme)),status='old',iostat=ios) 
if (ios/=0) write(*,*) 'problem opening '//trim(fnme)//' in iscanfile'

111 continue

read (7,'(a)',end=112) line

ieq=scan(line,'=')
if (ieq/=0) then
   if (trim(adjustl(line(1:ieq-1))).eq.text(1:len_trim(text)))  then
      nrec=len_trim(adjustl(line(ieq+1:len_trim(line))))
      record(1:nrec)=trim(adjustl(line(ieq+1:len_trim(line))))
        if (record(1:1).eq.text(1:1)) then
        res=0
        read (record(nrec-1:nrec),'(i2)',err=888) ires
        ires=-ires
        goto 889
888     read (record(nrec:nrec),'(i1)',err=999) ires
        ires=-ires
889     continue
        else
        icolon=scan(record(1:nrec),':')
          if (icolon.eq.0) then
          read (record(1:nrec),*,err=999) res
          else
          nd=nd+1
          read (record(1:icolon-1),*,err=999) range(1,nd)
          read (record(icolon+1:nrec),*,err=999) range(2,nd)
          res=par(nd)
          endif
        ires=1
        endif
      if (vocal.eq.1) then
      write (*,*) text//' - Value read is: '//record(1:nrec)
      elseif (vocal.eq.2) then
      write (*,*) text//' - Value read is: '//record(1:nrec)// &
      ' to replace default value : '//trim(cres0)
      elseif (vocal.eq.3) then
      write (*,*) text//' is '//trim(res_desc)//' - Value read is: '//record(1:nrec)// &
      ' to replace default value : '//trim(cres0)
      endif
   end if
end if

goto 111

112 close (7)

return

999 print*,''
    print*,'Error in input file'
    print*,'You have entered the following value : (',record(1:nrec),')'
    print*,'while parameter ('//text//') needs to be specified as an integer'
    print*,'Parameter Description: '//trim(res_desc)
    stop ' Terminating job'
end subroutine iscanfile



!---------------------------------------------------------------------------

subroutine dscanfile (fnme, text, res, res_desc, ires, vocal, nd, range, par)

!------------------------------------------------------------------------------|
! subroutine to read the value of a double precision  parameter whose name 
! is stored in  text  from file fnme the result is stored in res and the flag 
! ires is set  to 1 if all went well. if vocal is set yo 1 the routine echoes its find to
! standard output

character*(*) fnme,text,res_desc
double precision res
integer ires,vocal
integer,intent(inout)::nd
real*4,intent(inout)::range(2,1024),par(1024)

character line*256,record*256,cres0*256
integer ios,icolon,istar,itime,i

if (vocal.eq.4) then
write (*,*) text//' is '//trim(res_desc)
return
endif

ires=0
write (cres0,*) res
cres0=trim(adjustl(cres0))

open (7,file=fnme(1:len_trim(fnme)),status='old',iostat=ios) 
if (ios/=0) write(*,*) 'problem opening '//trim(fnme)//' in iscanfile'

111     continue

read (7,'(a)',end=112) line

ieq=scan(line,'=')
if (ieq/=0) then
   if (trim(adjustl(line(1:ieq-1))).eq.text(1:len_trim(text))) then
      nrec=len_trim(adjustl(line(ieq+1:len_trim(line))))
      record(1:nrec)=trim(adjustl(line(ieq+1:len_trim(line))))
      itime=scan(record(1:nrec),'time_topo')
        if (itime.ne.0) then
        read(record(itime+9:itime+10),'(i2)',err=991) i
991     read(record(itime+9:itime+9),'(i1)',err=999) i
        res=-i
        goto 113
        endif
        if (record(1:1).eq.text(1:1)) then
        res=0.d0
        read (record(nrec-1:nrec),'(i2)',err=888) ires
        ires=-ires
        goto 889
888     read (record(nrec:nrec),'(i1)',err=999) ires
        ires=-ires
889     continue
        else
        icolon=scan(record(1:nrec),':')
        istar=scan(record(1:nrec),'*')
          if (icolon.eq.0 .and. istar.eq.0 .and. itime.eq.0) then
          read (record(1:nrec),*,err=999) res
          ires=1
          elseif (icolon.ne.0) then
          nd=nd+1
          read (record(1:icolon-1),*,err=999) range(1,nd)
          read (record(icolon+1:nrec),*,err=999) range(2,nd)
          res=par(nd)
          ires=1
          else
          nrec=1
          record(1:1)='*'
          ires=999
          endif
        endif
113   continue
      if (vocal.eq.1) then
      write (*,*) text//' - Value read is: '//record(1:nrec)
      elseif (vocal.eq.2) then
      write (*,*) text//' - Value read is: '//record(1:nrec)// &
      ' to replace default value : '//trim(cres0)
      elseif (vocal.eq.3) then
      write (*,*) text//' is '//trim(res_desc)//' - Value read is: '//record(1:nrec)// &
      ' to replace default value : '//trim(cres0)
      endif
   endif
end if

goto 111

112   close (7)

return

999 print*,''
    print*,'Error in input file'
    print*,'You have entered the following value : (',record(1:nrec),')'
    print*,'while parameter ('//text//') needs to be specified as a floating point number'
    print*,'Parameter Description: '//trim(res_desc)
    stop ' Terminating job'
end subroutine dscanfile


!------------------------------------------------------------------------------|

subroutine cscanfile (fnme, text, res, res_desc, ires, vocal, nd, range, par)

!------------------------------------------------------------------------------|
! subroutine to read the value of a character string parameter whose name 
! is stored in  text  from file fnme the result is stored in res and the flag 
! ires is set  to 1 if all went well. if vocal is set yo 1 the routine echoes its find to
! standard output
!------------------------------------------------------------------------------|
      
character*(*) fnme,text,res,res_desc
integer ires,vocal
integer,intent(inout)::nd
real*4,intent(inout)::range(2,1024),par(1024)

character line*256,record*256,cres0*256
integer ios

if (vocal.eq.4) then
write (*,*) text//' is '//trim(res_desc)
return
endif

ires=0
cres0=res
cres0=trim(adjustl(cres0))

open (7,file=fnme(1:len_trim(fnme)),status='old',iostat=ios) 
if (ios/=0) write(*,*) 'problem opening '//trim(fnme)//' in iscanfile'

111     continue

read (7,'(a)',end=112) line

ieq=scan(line,'=')
if (ieq/=0) then
   if (trim(adjustl(line(1:ieq-1))).eq.text(1:len_trim(text))) then
      nrec=len_trim(adjustl(line(ieq+1:len_trim(line))))
      record(1:nrec)=trim(adjustl(line(ieq+1:len_trim(line))))
      res=record(1:nrec)
        if (vocal.eq.1) then
        write (*,*) text//' - Value read is: '//record(1:nrec)
        elseif (vocal.eq.2) then
        write (*,*) text//' - Value read is: '//record(1:nrec)// &
        ' to replace default value : '//trim(cres0)
        elseif (vocal.eq.3) then
        write (*,*) text//' is '//trim(res_desc)//' - Value read is: '//record(1:nrec)// &
        ' to replace default value : '//trim(cres0)
        endif
      ires=1
   endif
end if

goto 111

112   close (7)

return
end subroutine cscanfile


!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
