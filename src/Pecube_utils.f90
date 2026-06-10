! This file gather several subroutines useful for Pecube
! author: Maxime Bernard
! date: 17/10/2023

!--------------------------------------
subroutine GetOperatingSystem(is_unix)
  logical, intent(out) :: is_unix
  character(256) :: os_name
  integer :: status
  character(256) :: temp_os_name

  ! Get the value of the "OS" environment variable
  call get_environment_variable('OS', os_name, status)
  temp_os_name = os_name
  
  ! Check the value of the "OS" environment variable
  if (status == 0) then
    ! Convert the OS name to uppercase for case-insensitive comparison
    do i=1, len_trim(temp_os_name)
        temp_os_name(i:i) = char(ichar(temp_os_name(i:i), kind=1), kind=1)
    enddo
    ! Check if it's Windows or Unix
    if (index(temp_os_name, 'WINDOWS') > 0) then
      is_unix = .false.
      !print *, "OS system: Windows."
    else
      is_unix = .true.
      !print *, "OS system: Unix."
    end if
  else
    !print *, "Unable to determine the operating system. Assume Windows"
    is_unix = .false.
  end if
end subroutine GetOperatingSystem


!-----------------------------
subroutine CheckFileExistence(filename)
  character(len=*), intent(in) :: filename
  logical :: file_exists
  character(50) :: command

  inquire(file=filename, exist=file_exists)
  if (file_exists) then
    ! Windows system
    command = "del /Q "// filename
    status = system(command)
 endif
end subroutine CheckFileExistence

!----------------------------
subroutine RemoveDirectory(directory_path)
  character(*), intent(in) :: directory_path
  logical :: success
  character(256) :: cmd
  integer :: exit_status

  ! Initialize the success variable to false
  success = .false.

  ! Construct the system command to remove the directory
  ! The exact command may vary depending on the operating system
  cmd = "rmdir /s /q " // trim(directory_path) // " 2>nul" ! Windows

  ! Execute the system command to remove the directory
  exit_status = system(cmd)

  ! Check the exit status to determine success
  if (exit_status == 0) then
    success = .true.
    !print*, 'Remove directory.'
  else
    success = .false.
  end if
end subroutine RemoveDirectory

! ---------------------------------------------------------------------------------

subroutine sort_index(arr, indices, nstep, order)
        double precision :: arr(nstep)
        double precision :: sorted_arr(nstep)
        integer, intent(out) :: indices(nstep)
        double precision :: temp_val
        integer :: i, j, temp_idx, n, order

        ! Initiate variables
        n = size(arr)
        sorted_arr = arr
        indices = [(i, i=1, n)]

        ! sort array
        if (order.eq.0) then !'ascendant order'
            do i = 1, n-1
            do j = i+1, n
                if (sorted_arr(i) < sorted_arr(j)) then
                    ! sort values
                    temp_val = sorted_arr(i)
                    sorted_arr(i) = sorted_arr(j)
                    sorted_arr(j) = temp_val

                    ! sort corresponding indices
                    temp_idx = indices(i)
                    indices(i) = indices(j)
                    indices(j) = temp_idx
                end if
            end do
        end do
        else !Descendant order
            do i = 1, n-1
                do j = i+1, n
                    if (sorted_arr(i) > sorted_arr(j)) then
                        ! sort values
                        temp_val = sorted_arr(i)
                        sorted_arr(i) = sorted_arr(j)
                        sorted_arr(j) = temp_val

                        ! sort corresponding indices
                        temp_idx = indices(i)
                        indices(i) = indices(j)
                        indices(j) = temp_idx
                    end if
                end do
            end do
        endif
        
        arr = sorted_arr

    end subroutine sort_index
    
    
! ---------------------------------------------------------------------------------
subroutine isinarray(array,val2find, nstep,res)

  double precision :: array(nstep),val2find
  integer :: nstep, res, i
  
  res = 0
  do i=1, nstep
    if ((array(i)+1e-3 - val2find) * (array(i)-1e-3 - val2find).lt.0.d0) then
        res = 1
    endif
  enddo
  
end subroutine isinarray

!----------------------------------------------------------------------------------------------------
subroutine interp_1D(t1,t2,nt,z1,z2,nz)
   
   ! t1 = initial time array
   ! t2 = time array interpolated
   ! z1 = former array
   ! z2 = array to interpolate
   
    implicit none
    integer nt,nz
    double precision t1(nt),t2(nz),z1(nt),z2(nz)
    integer ti, zi
    double precision weight, tol

    tol = 1e-4
    
    ! Loop array to be interpolated
    ! t1(1) = 0
    do zi=1,nz
      ! loop through bigger array (minute time step)
        ti = 0
      do ti=2,nt
        if (t1(ti).ge.(t2(zi)-tol)) then
          weight = (z1(ti) - z1(ti-1)) / (t1(ti) - t1(ti-1))
          z2(zi) = z1(ti-1) + weight*(t2(zi) - t1(ti-1))
          goto 1
        endif
      enddo
      ! If issue print variables
      print*,'problem in interp_1D'
      print*, nt, nz, ti, zi
      print*, t2
      print*, t1
      print*, z1
      print*, z2
      print*, 't2 = ', t2(zi)
      print*, 't1 = ', t1(ti)
      print*, 't1 i-1 = ', t1(ti-1)
      print*, 'z1 i = ', z1(ti)
      print*, 'z1 i-1 = ', z1(ti-1)
      print*, 'zi = ', zi
      print*, 'ti = ', ti
      stop
      ! else continue
1      CONTINUE
    enddo
    z2(nz) = z2(nz-1)+1e-10
    
end subroutine interp_1D


!----------------------------------------------------------------------------------------------------
subroutine interpolate_1D(t1,t2,nt,z1,z2,nz)
   
   ! t1 = initial time array
   ! t2 = time array interpolated
   ! z1 = former array
   ! z2 = array to interpolate
   
    implicit none
    integer nt,nz
    double precision t1(nt),t2(nz),z1(nt),z2(nz)
    integer ti, zi
    double precision weight, tol

    tol = 1e-4
    
    ! Loop array to be interpolated
    ! t1(1) = 0
    do zi=1,nz
      ! loop through bigger array (minute time step)
        ti = 0
      do ti=2,nt
        if (t1(ti).ge.(t2(zi)-tol)) then
          weight = (z1(ti) - z1(ti-1)) / (t1(ti) - t1(ti-1))
          z2(zi) = z1(ti-1) + weight*(t2(zi) - t1(ti-1))
          goto 1
        endif
      enddo
      print*,'problem in interpolate_1D'
        print*, nt, nz, ti, zi
      print*, t2
      print*, t1
        print*, z1
        print*, z2
      print*, 't2 = ', t2(zi)
      print*, 't1 = ', t1(ti)
      print*, 't1 i-1 = ', t1(ti-1)
      print*, 'z1 i = ', z1(ti)
      print*, 'z1 i-1 = ', z1(ti-1)
      print*, 'zi = ', zi
      print*, 'ti = ', ti
      stop
1      CONTINUE
    enddo
    
end subroutine interpolate_1D


!----------------------------------------------------------------------------------------------------
subroutine interpolate2_1D(t1,t2,nt,z1,z2,nz)
   
  ! Interpolate form higher resolution array
  ! to lower resolution array
   ! t1 = initial time array (high res)
   ! t2 = time array (low res)
   ! z1 = former array (high res)
   ! z2 = array to interpolate (low res)
   
    implicit none
    integer nt,nz
    double precision t1(nt),t2(nz),z1(nt)
    double precision z2(nz)
    integer ti, zi
    double precision weight, tol

    tol = 1e-4
    
    do zi=1,nz
      ! loop through low res array
        ti = 0
      do ti=1,nt ! loop through high res array
        if (t2(zi).ge.(t1(ti)-tol)) then
          if (ti.eq.1) then
            z2(zi) = z1(1)
          elseif (ti.eq.nt) then
            z2(zi) = z1(nt)
          else
            weight = (z1(ti-1) - z1(ti)) / (t1(ti-1) - t1(ti))
            z2(zi) = z1(ti-1) + weight*(t1(ti-1)-t2(zi))
          endif
          goto 1
        endif
      enddo
      print*,'problem in interpolate2_1D'
        print*, nt, nz, ti, zi
      print*, 't2 = ', t2(zi)
      print*, 't1 = ', t1(ti)
      print*, 't1 i-1 = ', t1(ti-1)
      print*, 'z1 i = ', z1(ti)
      print*, 'z1 i-1 = ', z1(ti-1)
      print*, 'zi = ', zi
      print*, 'ti = ', ti
      stop
1      CONTINUE
    enddo
    
end subroutine interpolate2_1D

!###########################################################

SUBROUTINE InterpolateTTPathKet(numTTDefs, time_in, temp_in, numTTNodes, time_out, temp_out,numTTNodes_init,&
            Temp_crit, status)
  ! Code translated from the c routine of Ketcham et al. (1999) by Maxime Bernard
  ! Date: 12/09/2025

  ! InterpolateTTPathKet
  ! Takes the time-temperature path specification and subdivides it for
  ! calculation in isothermal intervals.
  ! Does it based on model of Ketcham et al., in review.
  ! It is calibrated to facilitate 0.5% accuracy for end-member F-apatite by
  ! having a maximum temperature step of 3.5 degrees C when the model
  ! temperature is within 10 C of the total annealing temperature.  Before this
  ! cutoff the maximum temperature step required is 8 C.  If the overall model
  ! time steps are too large, these more distant requirements may not be met.
    
  IMPLICIT NONE
  ! Inputs
  integer numTTDefs,numTTNodes  !numTTNodes = MAX_NUM_TIME_STEP
  double precision time_in(numTTDefs),temp_in(numTTDefs),time_out(numTTNodes),temp_out(numTTNodes)
  real*4 pctPerTimeStep
  real*4 NEAR_ANNEAL_CUTOFF_KET 
  real*4 MAX_TEMP_STEP_NEAR_TA_KET 
  real*4 MAX_TEMP_STEP_KET 
  real*4 KELVINS_AT_0C 
  real*4 SECS_PER_MA
  real*4 Temp_crit ! Critical temperature to decrease time step (depends on thermochronometer)

  ! Outputs
  integer numTTNodes_init 
  integer status  ! return 1 = success, 0 = failure

  ! Locals
  integer dN, n
  real*8 rate, absRate
  real*8 maxTMult, maxTemp
  real*8 nearAnnealTemp, timeStep
  real*8 defTimeStep, tempPerTimeStep
  real*8 currDefTimeStep, altTimeStep
  real*8 endTemp, min_dt

  ! Initialize
  numTTNodes_init = 1
  status = 1
  altTimeStep = 0.0
  pctPerTimeStep = 0.2 ! Needed for He rhov - Minimum time step (Myr)
  NEAR_ANNEAL_CUTOFF_KET = 10.0
  MAX_TEMP_STEP_NEAR_TA_KET = 3.5
  MAX_TEMP_STEP_KET = 8.0
  KELVINS_AT_0C = 273.15
  SECS_PER_MA = 3.15576E13

  temp_out(1) = temp_in(numTTDefs) + KELVINS_AT_0C
  time_out(1) = time_in(numTTDefs) 
  defTimeStep = time_in(numTTDefs) * pctPerTimeStep / 100.0
  min_dt = 0.2 ! Needed for He rhov - Minimum time step (Myr)

  ! Loop through segments (backward)
  do dN = numTTDefs, 2, -1
    rate = (temp_in(dN) - temp_in(dN-1)) / &
          (time_in(dN) - time_in(dN-1) + 1.0E-4)
    absRate = ABS(rate)
    tempPerTimeStep = absRate * defTimeStep

    if (tempPerTimeStep <= MAX_TEMP_STEP_KET) then
      currDefTimeStep = defTimeStep
    elseif (temp_in(dN-1) < Temp_crit) then ! if we passed below critical temperature
      currDefTimeStep = min_dt
    else
      currDefTimeStep = MAX_TEMP_STEP_KET / absRate
    endif

    if (rate > 0.0) then
      maxTMult = 0.0
    else
      maxTMult = -1.0
    endif
    endTemp = temp_in(dN-1) + KELVINS_AT_0C

    ! Near anneal temp
    if (absRate < 0.1) then
      nearAnnealTemp = 1000.0
    else
      nearAnnealTemp = 3.7767 * absRate**0.019837 - NEAR_ANNEAL_CUTOFF_KET
      altTimeStep = MAX_TEMP_STEP_NEAR_TA_KET / absRate
    endif

    ! While loop: subdivide segment
    do while (time_out(numTTNodes_init) > time_in(dN-1))
      if (numTTNodes_init + 1 > numTTNodes) then
        status = 0
        temp_out = temp_out - KELVINS_AT_0C
        return
      endif

      maxTemp = temp_out(numTTNodes_init) + defTimeStep * rate * maxTMult
      if ((rate < 0.0) .and. (maxTemp > endTemp)) maxTemp = endTemp

      timeStep = currDefTimeStep
      if (maxTemp > nearAnnealTemp) then
        if (altTimeStep < defTimeStep) timeStep = altTimeStep
      endif

      ! Final step check
      if (timeStep + 0.001 > time_out(numTTNodes_init) - time_in(dN-1)) then
        time_out(numTTNodes_init+1) = time_in(dN-1)
        temp_out(numTTNodes_init+1) = endTemp
      else
        time_out(numTTNodes_init+1) = time_out(numTTNodes_init) - timeStep
        temp_out(numTTNodes_init+1) = temp_out(numTTNodes_init) - rate * timeStep
      endif

      numTTNodes_init = numTTNodes_init + 1
    enddo
  enddo

  ! Convert Ma → seconds
  do n = 1, numTTNodes_init
    time_out(n) = time_out(n) !* SECS_PER_MA
    temp_out(n) = temp_out(n) - KELVINS_AT_0C
  enddo

end subroutine InterpolateTTPathKet

