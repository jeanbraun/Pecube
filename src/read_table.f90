subroutine read_table(fnme, tag0, ntags0, sample_names, sample_values, nsamples, nsample_max)

! this routine reads from a table from a file

! it transfers into an array "sample_values" the values of the fields (in the table) corresponding to the
! tags given in tag(ntags) in the same order as they are given in the tag array

! sample_names returns the name of the samples
! nsamples returns the number of samples

! the table in the file should start with a column labeled SAMPLE

! extra field are neglected
! if a requested field value is missing (or an entire field) it will appear as -9999.


character*1024 line
character*128 fnme, tag0(ntags0), sample_names(nsample_max)
double precision sample_values(ntags0, nsample_max)
integer nsamples

integer, dimension(:), allocatable :: ind
integer ntags
character*128, dimension(:), allocatable :: tag

type :: record
character :: name*128
double precision :: x(10000)
end type record

type(record) :: rec

ntags = ntags0 + 1

allocate(ind(ntags), tag(ntags))
tag(1) = 'SAMPLE'
tag(2:ntags) = tag0

open (777, file=trim(fnme), status='old')

ind=0

    do while (ind(1).eq.0)
    read (777, '(a)') line
    itab = 1
    istart = 1
    i = index(line(istart:), ',')
        do k=1,ntags
        if (adjustl(trim(tag(k))).eq.adjustl(trim(line(istart:istart+i-2)))) ind(k)=itab
        enddo
        do while (i.ne.0)
        itab = itab + 1
        istart = istart + i
        i = index(line(istart:), ',')
            if (i.ne.0) then
                do k=1,ntags
                if (adjustl(trim(tag(k))).eq.adjustl(trim(line(istart:istart+i-2)))) ind(k)=itab
                enddo
            else
                do k=1,ntags
                if (adjustl(trim(tag(k))).eq.adjustl(trim(line(istart:)))) ind(k)=itab
                enddo
            endif
        enddo
    enddo

    if (ind(1).ne.1) then
    print*,'ERROR'
    print*,'The first non-empty column should be the sample name and thus labeled SAMPLE'
    goto 999
    endif

nsamples = 0
    do
    read (777, '(a)', end=999) line
        if (adjustl(trim(line)).ne.'') then
        nsamples = nsamples + 1
        rec%x = -9999.
        read (line,*, end=998) rec
998     sample_names(nsamples) = adjustl(trim(rec%name))
            do k=1,ntags0
                if (ind(k+1).ne.0) then
                sample_values(k,nsamples) = rec%x(ind(k+1)-1)
                else
                sample_values(k,nsamples) = -9999.
                endif
            enddo
        endif
    enddo

999 close (777)

deallocate (ind, tag)


end subroutine read_table
