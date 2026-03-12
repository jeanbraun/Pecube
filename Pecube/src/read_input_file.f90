!--------------------------------------------------------------------------------------------

subroutine read_input_file (fnme, vocal, p, nd, range, par)

! This routine reads all input parameters from an input file (fnme)
! They are placed in the p (param) structure for use in Pecube
! All parameters have a default value
!
! Each time a pair of parameter values separated by a column are found in the input file,
! the nd parameter is incremented by one and the resulting range is stored in the range array
! and the value returned is taken from par(nd)

! If vocal is 1, each parameter value read from the input file is echoed to the standard output
! If vocal is 2, each parameter value read from the input file and the default value are echoed to the standard output

! This is where a couple of lines must be added per new parameter

use Pecube

implicit none

character*(*) :: fnme
integer :: vocal, ij, num, k, istep1, istep2
type (parameters) p
integer nd, itime
real*4 range(2,1024), par(1024)

integer :: res, i, j, npoint_max, nstep_max, j2
character*10 ci, cj
character*100 text, text_temp
integer :: tecto_indices(100)
double precision :: tecto_values(100),new_par_val


call scanfile (fnme, "echo_input_file", p%echo_input_file, p%echo_input_file_desc, res, vocal, nd, range, par)

call scanfile (fnme, "topo_file_name", p%topo_file_name, p%topo_file_name_desc, res, vocal, nd, range, par)

call scanfile (fnme, "nx", p%nx, p%nx_desc, res, vocal, nd, range, par)

call scanfile (fnme, "ny", p%ny, p%ny_desc, res, vocal, nd, range, par)

call scanfile (fnme, "nz", p%nz, p%nz_desc, res, vocal, nd, range, par)

call scanfile (fnme, "lon0", p%lon0, p%lon0_desc, res, vocal, nd, range, par)

call scanfile (fnme, "lat0", p%lat0, p%lat0_desc, res, vocal, nd, range, par)

call scanfile (fnme, "dlon", p%dlon, p%dlon_desc, res, vocal, nd, range, par)

call scanfile (fnme, "dlat", p%dlat, p%dlat_desc, res, vocal, nd, range, par)

call scanfile (fnme, "nskip", p%nskip, p%nskip_desc, res, vocal, nd, range, par)

call scanfile (fnme, "thickness", p%thickness, p%thickness_desc, res, vocal, nd, range, par)

call scanfile (fnme, "ntime", p%ntime, p%ntime_desc, res, vocal, nd, range, par)

call scanfile (fnme, "erosional_time_scale", p%erosional_time_scale, p%erosional_time_scale_desc, res, vocal, nd, range, par)

call scanfile (fnme, "topo_ref", p%topo_ref, p%topo_ref_desc, res, vocal, nd, range, par) !Maxime

call scanfile (fnme, "topo_ref_custom", p%topo_ref_custom, p%topo_ref_custom_desc, res, vocal, nd, range, par) !Maxime

call scanfile (fnme, "topo_wavelength", p%topo_wavelength, p%topo_wavelength_desc, res, vocal, nd, range, par) !Maxime

call scanfile (fnme, "topo_offset", p%topo_offset, p%topo_offset_desc, res, vocal, nd, range, par) !Maxime

call scanfile (fnme, "topo_amp", p%topo_amp, p%topo_amp_desc, res, vocal, nd, range, par) !Maxime

call scanfile (fnme, "topo_phase", p%topo_phase, p%topo_phase_desc, res, vocal, nd, range, par) !Maxime

allocate (p%time_topo(p%ntime+1), p%amplification(p%ntime+1), p%offset(p%ntime+1), p%output(p%ntime+1),&
            p%phase(p%ntime+1))
p%time_topo(:) = 0.d0
p%amplification(:) = 1.d0
p%offset(:) = 0.d0
p%output(:) = 1
p%phase(:) = 0.d0

  do i = 1, p%ntime + 1
    write (ci,'(i10)') i
    j = i
    call scanfile (fnme, "time_topo"//trim(adjustl(ci)), p%time_topo(j), p%time_topo_desc, res, vocal, nd, range, par)
    if (res.eq.999) then
        p%time_topo(j)=p%time_topo(i-1)!Maxime
    else
     if (j.gt.1) then ! Check if previous time is older than current time
      if (p%time_topo(j) .gt. p%time_topo(j-1)) then
        range(2,nd) = range(2,nd-1)
        new_par_val = 0
        call rescale_parameter (p%time_topo(j-1),range,nd, new_par_val,0)
        par(nd) = new_par_val ! make sure that the new inverted parameter is stored in param to write the value in the NA_int_results.csv file
      endif
     endif
    endif
  enddo
  
  do i = 1, p%ntime + 1
      write (ci,'(i10)') i
      j = i
      call scanfile (fnme, "amplification"//trim(adjustl(ci)), p%amplification(j), p%amplification_desc, res, vocal, nd, range, par)
      if (res.eq.999) then
          p%amplification(j)=p%amplification(j-1) !Maxime
      elseif (res.eq.996) then ! Take as maximum value
        range(2,nd) = range(2,nd-1)
          if (p%amplification(j).gt.p%amplification(j-1)) then
            new_par_val = 0
            call rescale_parameter (p%amplification(j-1),range,nd, new_par_val,0)
            par(nd) = new_par_val ! make sure that the new inverted parameter is stored in param to write the value in the NA_int_results.csv file
        endif
      elseif (res.eq.997) then ! Take as minimum value
        range(1,nd) = range(1,nd-1)
        if (p%amplification(j).lt.p%amplification(j-1)) then
            new_par_val = 0
            call rescale_parameter (p%amplification(j-1),range,nd, new_par_val,1)
            par(nd) = new_par_val ! make sure that the new inverted parameter is stored in param to write the value in the NA_int_results.csv file
        endif
      endif
  enddo
  
   do i = 1, p%ntime + 1
      write (ci,'(i10)') i
      j = i
    call scanfile (fnme, "offset"//trim(adjustl(ci)), p%offset(j), p%offset_desc, res, vocal, nd, range, par)
    if (res.eq.999) then
        p%offset(j)=p%offset(j-1)!Maxime
    elseif (res.eq.996) then
        range(2,nd) = range(2,nd-1)
        if (p%offset(j).gt.p%offset(j-1)) then
            new_par_val = 0
            call rescale_parameter (p%offset(j-1),range,nd, new_par_val,0)
            par(nd) = new_par_val ! make sure that the new inverted parameter is stored in param to write the value in the NA_int_results.csv file
        endif
    endif
    enddo
    
  do i = 1, p%ntime + 1
      write (ci,'(i10)') i
      j = i
      call scanfile (fnme, "output"//trim(adjustl(ci)), p%output(j), p%output_desc, res, vocal, nd, range, par)
  enddo
  
  do i = 1, p%ntime + 1
      write (ci,'(i10)') i
      j = i
      call scanfile (fnme, "phase"//trim(adjustl(ci)), p%phase(j), p%phase_desc, res, vocal, nd, range, par)
      if (res.eq.999) then
        p%phase(j)=p%phase(j-1) !Maxime
      elseif (res.eq.996) then
        range(2,nd) = range(2,nd-1)
        if (p%phase(j).gt.p%phase(j-1)) then
            new_par_val = 0
            call rescale_parameter (p%phase(j-1),range,nd, new_par_val,0)
            par(nd) = new_par_val ! make sure that the new inverted parameter is stored in param to write the value in the NA_int_results.csv file
        endif
    endif
  enddo

!do i = 1, p%ntime
!if (p%time_topo(i).eq.0.d0.and.vocal.ne.4) stop 'There must be as many time_topo as the number of time steps'
!enddo

call scanfile (fnme, "do_headward", p%do_Headward, p%do_Headward_desc, res, vocal, nd, range, par) !Maxime

call scanfile (fnme, "tauH", p%tauH, p%tauH_desc, res, vocal, nd, range, par) !Maxime

call scanfile (fnme, "depth_incision", p%depth_incision, p%depth_incision_desc, res, vocal, nd, range, par) !Maxime

call scanfile (fnme, "time_incision_start", p%time_incision_start, p%time_incision_start_desc, res, vocal, nd, range, par) !Maxime

call scanfile (fnme, "time_incision_stop", p%time_incision_stop, p%time_incision_stop_desc, res, vocal, nd, range, par) !Maxime

call scanfile (fnme, "isostasy", p%isostasy, p%isostasy_desc, res, vocal, nd, range, par)

call scanfile (fnme, "rho_crust", p%rho_crust, p%rho_crust_desc, res, vocal, nd, range, par)

call scanfile (fnme, "rho_asthenosphere", p%rho_asthenosphere, p%rho_asthenosphere_desc, res, vocal, nd, range, par)

call scanfile (fnme, "young_modulus", p%young_modulus, p%young_modulus_desc, res, vocal, nd, range, par)

call scanfile (fnme, "poisson_ratio", p%poisson_ratio, p%poisson_ratio_desc, res, vocal, nd, range, par)

call scanfile (fnme, "EET", p%EET, p%EET_desc, res, vocal, nd, range, par)

call scanfile (fnme, "nx_isostasy", p%nx_isostasy, p%nx_isostasy_desc, res, vocal, nd, range, par)

call scanfile (fnme, "ny_isostasy", p%ny_isostasy, p%ny_isostasy_desc, res, vocal, nd, range, par)

call scanfile (fnme, "thermal_diffusivity", p%thermal_diffusivity, p%thermal_diffusivity_desc, res, vocal, nd, range, par)

call scanfile (fnme, "basal_temperature", p%basal_temperature, p%basal_temperature_desc, res, vocal, nd, range, par)

call scanfile (fnme, "min_dt", p%min_dt, p%min_dt_desc, res, vocal, nd, range, par) !Maxime

call scanfile (fnme, "D_heat", p%D_heat, p%D_heat_desc, res, vocal, nd, range, par)

call scanfile (fnme, "sea_level_temperature", p%sea_level_temperature, p%sea_level_temperature_desc, res, vocal, nd, range, par)

call scanfile (fnme, "lapse_rate", p%lapse_rate, p%lapse_rate_desc, res, vocal, nd, range, par)

call scanfile (fnme, "heat_production", p%heat_production, p%heat_production_desc, res, vocal, nd, range, par)

call scanfile (fnme, "heatprod_flag", p%heatprod_flag, p%heatprod_flag_desc, res, vocal, nd, range, par)

call scanfile (fnme, "data_folder", p%data_folder, p%data_folder_desc, res, vocal, nd, range, par)

call scanfile (fnme, "default_age", p%default_age, p%default_age_desc, res, vocal, nd, range, par)

call scanfile (fnme, "FT_code_flag", p%FT_code_flag, p%FT_code_flag_desc, res, vocal, nd, range, par)

call scanfile (fnme, "misfit_slope", p%misfit_slope, p%misfit_slope_desc, res, vocal, nd, range, par)

call scanfile (fnme, "misfit_function", p%misfit_function, p%misfit_function_desc, res, vocal, nd, range, par)

call scanfile (fnme, "fault_advect_flag", p%fault_advect_flag, p%fault_advect_flag_desc, res, vocal, nd, range, par)

call scanfile (fnme, "shear_heating", p%shear_heating, p%shear_heating_desc, res, vocal, nd, range, par)

call scanfile (fnme, "age_TL_flag", p%age_TL_flag, p%age_TL_flag_desc, res, vocal, nd, range, par)

call scanfile (fnme, "age_OSL_flag", p%age_OSL_flag, p%age_OSL_flag_desc, res, vocal, nd, range, par)

call scanfile (fnme, "age_ESR_flag", p%age_ESR_flag, p%age_ESR_flag_desc, res, vocal, nd, range, par)   

call scanfile (fnme, "ESR_misfit_target", p%ESR_misfit_Target, p%ESR_misfit_Target_desc, res, vocal, nd, range, par)    

call scanfile (fnme, "OSL_misfit_target", p%OSL_misfit_Target, p%OSL_misfit_Target_desc, res, vocal, nd, range, par)
                                                                                                 
call scanfile (fnme, "age_AHe_flag", p%age_AHe_flag, p%age_AHe_flag_desc, res, vocal, nd, range, par)

call scanfile (fnme, "age_ZHe_flag", p%age_ZHe_flag, p%age_ZHe_flag_desc, res, vocal, nd, range, par)

call scanfile (fnme, "age_AFT_flag", p%age_AFT_flag, p%age_AFT_flag_desc, res, vocal, nd, range, par)

call scanfile (fnme, "Init_FTL_Model", p%Init_FTL_Model, p%Init_FTL_Model_desc, res, vocal, nd, range, par)

call scanfile (fnme, "Init_FTL_value", p%Init_FTL_value, p%Init_FTL_value_desc, res, vocal, nd, range, par)

call scanfile (fnme, "Kinetic_FTL_Parameter", p%Kinetic_FTL_Parameter, p%Kinetic_FTL_Parameter_desc, res, vocal, nd, range, par)

call scanfile (fnme, "KINAFT", p%Kinetic_FTL_Parameter_value_AFT, p%Kinetic_FTL_Parameter_value_AFT_desc,&
                res, vocal, nd, range, par)

call scanfile (fnme, "KINAHE", p%Kinetic_FTL_Parameter_value_AHe, p%Kinetic_FTL_Parameter_value_AHe_desc,&
                res, vocal, nd, range, par)

call scanfile (fnme, "MFTL_error", p%MFTL_error, p%MFTL_error_desc,res, vocal, nd, range, par)

call scanfile (fnme, "MFTL_std_error", p%MFTL_std_error, p%MFTL_std_error_desc,res, vocal, nd, range, par)

call scanfile (fnme, "age_ZFT_flag", p%age_ZFT_flag, p%age_ZFT_flag_desc, res, vocal, nd, range, par)

call scanfile (fnme, "age_FTL_flag", p%age_FTL_flag, p%age_FTL_flag_desc, res, vocal, nd, range, par)
if (p%age_FTL_flag.eq.1) p%age_AFT_flag = 1

call scanfile (fnme, "age_KAr_flag", p%age_KAr_flag, p%age_KAr_flag_desc, res, vocal, nd, range, par)

call scanfile (fnme, "age_BAr_flag", p%age_BAr_flag, p%age_BAr_flag_desc, res, vocal, nd, range, par)

call scanfile (fnme, "age_MAr_flag", p%age_MAr_flag, p%age_MAr_flag_desc, res, vocal, nd, range, par)

call scanfile (fnme, "age_HAr_flag", p%age_HAr_flag, p%age_HAr_flag_desc, res, vocal, nd, range, par)

call scanfile (fnme, "He43_flag", p%He43_flag, p%He43_flag_desc, res, vocal, nd, range, par)

call scanfile (fnme, "He43_inv_flag", p%He43_inv_flag, p%He43_inv_flag_desc, res, vocal, nd, range, par)

!!! For sample-specific parameters - Maxime
call scanfile (fnme, "rhoST", p%rhoST, p%rhoST_desc, res, vocal, nd, range, par)

call scanfile (fnme, "age_computation", p%age_computation, p%age_computation_desc, res, vocal, nd, range, par)

call scanfile (fnme, "grainradius", p%grainradius, p%grainradius_desc, res, vocal, nd, range, par)

call scanfile (fnme, "ASIZE", p%ASize, p%ASize_desc, res, vocal, nd, range, par)

call scanfile (fnme, "AUPPM", p%AUPPM, p%AUppm_desc, res, vocal, nd, range, par)

call scanfile (fnme, "AThPPM", p%AThPPM, p%AThppm_desc, res, vocal, nd, range, par)

call scanfile (fnme, "RDmodel", p%RDmodel, p%RDmodel_desc, res, vocal, nd, range, par)

call scanfile (fnme, "D0_AHe", p%D0_AHe, p%D0_AHe_desc, res, vocal, nd, range, par)

call scanfile (fnme, "Ea_AHe", p%Ea_AHe, p%Ea_AHe_desc, res, vocal, nd, range, par)

call scanfile (fnme, "D0_ZHe", p%D0_ZHe, p%D0_ZHe_desc, res, vocal, nd, range, par)

call scanfile (fnme, "Ea_ZHe", p%Ea_ZHe, p%Ea_ZHe_desc, res, vocal, nd, range, par)

call scanfile (fnme, "RDmodelz", p%RDmodelz, p%RDmodelz_desc, res, vocal, nd, range, par)

call scanfile (fnme, "ZSize", p%ZSize, p%ZSize_desc, res, vocal, nd, range, par)

call scanfile (fnme, "ZUPPM", p%ZUppm, p%ZUppm_desc, res, vocal, nd, range, par)

call scanfile (fnme, "ZThPPM", p%ZThppm, p%ZThppm_desc, res, vocal, nd, range, par)

call scanfile (fnme, "D0_KAr", p%D0_KAr, p%D0_KAr_desc, res, vocal, nd, range, par)

call scanfile (fnme, "Ea_KAr", p%Ea_KAr, p%Ea_KAr_desc, res, vocal, nd, range, par)

call scanfile (fnme, "D0_BAr", p%D0_BAr, p%D0_BAr_desc, res, vocal, nd, range, par)

call scanfile (fnme, "Ea_BAr", p%Ea_BAr, p%Ea_BAr_desc, res, vocal, nd, range, par)

call scanfile (fnme, "D0_MAr", p%D0_MAr, p%D0_MAr_desc, res, vocal, nd, range, par)

call scanfile (fnme, "Ea_MAr", p%Ea_MAr, p%Ea_MAr_desc, res, vocal, nd, range, par)

call scanfile (fnme, "D0_HAr", p%D0_HAr, p%D0_HAr_desc, res, vocal, nd, range, par)

call scanfile (fnme, "Ea_HAr", p%Ea_HAr, p%Ea_HAr_desc, res, vocal, nd, range, par)

call scanfile (fnme, "Alpha_Ejec_Flag", p%Alpha_Ejec_Flag, p%Alpha_Ejec_Flag_desc, res, vocal, nd, range, par)

call scanfile (fnme, "ESR_Model", p%ESR_Model, p%ESR_Model_desc, res, vocal, nd, range, par) !Maxime

call scanfile (fnme, "ESR_Lmax", p%ESR_Lmax, p%ESR_Lmax_desc, res, vocal, nd, range, par)

call scanfile (fnme, "OSL_Model", p%OSL_Model, p%OSL_Model_desc, res, vocal, nd, range, par) !Maxime

call scanfile (fnme, "OSL_a", p%OSL_a, p%OSL_a_desc, res, vocal, nd, range, par) !Maxime

call scanfile (fnme, "OSL_Lmax", p%OSL_Lmax, p%OSL_Lmax_desc, res, vocal, nd, range, par)

call scanfile (fnme, "TL_Model", p%TL_Model, p%TL_Model_desc, res, vocal, nd, range, par) !Maxime

!!! end Maxime

call scanfile (fnme, "nfault", p%nfault, p%nfault_desc, res, vocal, nd, range, par)
if (vocal.eq.4) p%nfault = 1

call scanfile (fnme, "lon1", p%x1, p%x1_desc, res, vocal, nd, range, par)

call scanfile (fnme, "lat1", p%y1, p%y1_desc, res, vocal, nd, range, par)

call scanfile (fnme, "lon2", p%x2, p%x2_desc, res, vocal, nd, range, par)

call scanfile (fnme, "lat2", p%y2, p%y2_desc, res, vocal, nd, range, par)

!write (*,*) p%npoint
allocate (p%npoint(p%nfault), p%nstep(p%nfault), p%static(p%nfault))
p%npoint = 0
p%nstep = 0
p%static = 0            
  do i = 1, p%nfault
  write (ci,'(i10)') i
  call scanfile (fnme, "npoint"//trim(adjustl(ci)), p%npoint(i), p%npoint_desc, res, vocal, nd, range, par)
  call scanfile (fnme, "nstep"//trim(adjustl(ci)), p%nstep(i), p%nstep_desc, res, vocal, nd, range, par)
  call scanfile (fnme, "static"//trim(adjustl(ci)), p%static(i), p%static_desc, res, vocal, nd, range, par)
  enddo
if (vocal.eq.4) p%npoint(1)=2
if (vocal.eq.4) p%nstep(1)=2

call scanfile (fnme, "bottom_left", p%bottom_left, p%bottom_left_desc, res, vocal, nd, range, par)

call scanfile (fnme, "bottom_right", p%bottom_right, p%bottom_right_desc, res, vocal, nd, range, par)

call scanfile (fnme, "top_right", p%upper_right, p%upper_right_desc, res, vocal, nd, range, par)

call scanfile (fnme, "top_left", p%upper_left, p%upper_left_desc, res, vocal, nd, range, par)

npoint_max = maxval(p%npoint)
if (npoint_max.lt.0) npoint_max = 4
nstep_max = maxval(p%nstep)

allocate (p%r(npoint_max, p%nfault), p%s(npoint_max, p%nfault))
p%r = 0.d0
p%s = 0.d0
  do i = 1, p%nfault
  write (ci,'(i10)') i
    if (p%npoint(i).lt.0) then
! these lines were commented out to allow for combining two velocity fields
! one defined by a fault and one defined by a uniform velocity field (vertical only)
! These lines were a legacy from Pecube v3                                                                                                                                            
!      do j = 1, 4
!      write (cj,'(i10)') j
!      call scanfile (fnme, "r"//trim(adjustl(ci))//"_"//trim(adjustl(cj)), p%r(j,i), p%r_desc, res, vocal, nd, range, par)
!      if (res.eq.999) p%r(j,i)=p%r(j-1,i)
!      enddo
    else
      do j = 1, p%npoint(i)
      write (cj,'(i10)') j
      call scanfile (fnme, "r"//trim(adjustl(ci))//"_"//trim(adjustl(cj)), p%r(j,i), p%r_desc, res, vocal, nd, range, par)
      if (res.eq.999) p%r(j,i)=p%r(j-1,i)
      call scanfile (fnme, "s"//trim(adjustl(ci))//"_"//trim(adjustl(cj)), p%s(j,i), p%s_desc, res, vocal, nd, range, par)
      if (res.eq.999) p%s(j,i)=p%s(j-1,i)
      enddo
    endif
  enddo

allocate (p%time_start(nstep_max, p%nfault), p%time_end(nstep_max, p%nfault), p%velo(nstep_max, p%nfault))
p%time_start = 0.d0
p%time_end = 0.d0
p%velo = 0.d0
j2 = 1
ij = 0
tecto_indices(:) = 0
  do i = 1, p%nfault
  write (ci,'(i10)') i
    do j = 1, p%nstep(i)
    ! Record indices to recover tectonic timing
    ij = ij+1
    text=""
    write(text_temp,'(I1)') i
    text = trim(text) //trim(text_temp)
    write(text_temp,'(I1)') j
    text = trim(text) // trim(text_temp)
    read(text, '(I2)'), num
    tecto_indices(ij) = num
    write (cj,'(i10)') j
    
    
    call scanfile (fnme, "time_start"//trim(adjustl(ci))//"_"//trim(adjustl(cj)), p%time_start(j,i), p%time_start_desc, &
                   res, vocal, nd, range, par)
    ! if * time start from end previous step (same fault)
    if (res.eq.999) p%time_start(j,i)=p%time_end(j-1,i) 
  ! if # time start from end of first fault
    if (res.eq.998) p%time_start(j,i)=p%time_end(j2,i-1) 
    if (j.gt.1) then
        if (p%time_start(j,i).gt.p%time_end(j-1,i)) then ! Handle if time start >= time end step before
            new_par_val = 0
            call rescale_parameter (p%time_end(j-1,i)-1e-6,range,nd, new_par_val,0)
            par(nd) = new_par_val ! make sure that the new inverted parameter is stored in param to write the value in the NA_int_results.csv file
        endif
    endif
    if (p%time_start(j,i).lt.0.d0) then
         text = ""
         p%time_start(j,i)=p%time_topo(int(-p%time_start(j,i)))
         write(text_temp, '(I2)') j
         text = trim('time_start') // trim(text_temp) //'_'
         write(text_temp,'(I2)') i
         text = trim(text) // trim(text_temp) // ' = '
         write(text_temp,'(F6.2)') p%time_start(j,i)
         text = trim(text) // trim(text_temp)
         if (vocal.gt.1) print *,trim(text)
    endif
    
    call scanfile (fnme, "time_end"//trim(adjustl(ci))//"_"//trim(adjustl(cj)), p%time_end(j,i), p%time_end_desc, &
                   res, vocal, nd, range, par)
    if (j.gt.1) then
        if (p%time_end(j,i).gt.p%time_start(j,i)) then
            !print *, 'Apply change in time_end'
            new_par_val = 0
            call rescale_parameter (p%time_start(j,i)-1e-6,range,nd, new_par_val,0)
            par(nd) = new_par_val ! make sure that the new inverted parameter is stored in param to write the value in the NA_int_results.csv file
        endif
    endif

    if (p%time_end(j,i).lt.0.d0) then
        p%time_end(j,i)=p%time_topo(int(-p%time_end(j,i)))
        text = ""
        write(text_temp, '(I2)') j
        text = trim('time_start') // trim(text_temp) // '_'
        write(text_temp,'(I2)') i
        text = trim(text) // trim(text_temp) // ' = '
        write(text_temp,'(F6.2)') p%time_end(j,i)
        text = trim(text) // trim(text_temp)
        if (vocal.gt.1) print *,trim(text)
    endif
  enddo
 enddo
    
  ! Check again for time topo values
    do i = 1, p%ntime + 1
        j = i
        if (p%time_topo(j).lt.0.d0) then
            do k = 1, 30
                if(abs(p%time_topo(j)).eq.tecto_indices(k)) then
                    p%time_topo(j) = tecto_values(k)
                    goto 113
                endif
            enddo
            text = ""
            write(text_temp, '(I2)') j
            text = trim('time_topo') // trim(text_temp)
            write(text_temp,'(F6.2)') p%time_topo(j)
            text = trim(text) // trim(text_temp)
            if (vocal.gt.1) print *,trim(text)
        endif
113     continue
    enddo
    
    do i = 1, p%nfault
        write (ci,'(i10)') i
        do j = 1, p%nstep(i)
            write (cj,'(i10)') j
            ! Record indices to recover tectonic timing
            ij = ij+1
            call scanfile (fnme, "velo"//trim(adjustl(ci))//"_"//trim(adjustl(cj)), p%velo(j,i), p%velo_desc, &
                      res, vocal, nd, range, par)
            if (res.eq.996) then
                range(2,nd) = range(2,nd-1)
                if (p%velo(j,i).gt.p%velo(j-1,i)) then ! force velocity to decrease relative to value from previous step
                    new_par_val = 0
                    call rescale_parameter (p%velo(j-1,i),range,nd, new_par_val,0)
                    par(nd) = new_par_val ! make sure that the new inverted parameter is stored in param to write the value in the NA_int_results.csv file
                endif
            elseif (res.eq.997) then ! Take as minimum value
                range(1,nd) = range(1,nd-1)
                if (p%velo(j,i).lt.p%velo(j-1,i)) then
                    new_par_val = 0
                    call rescale_parameter (p%velo(j-1,i),range,nd, new_par_val,1)
                    par(nd) = new_par_val ! make sure that the new inverted parameter is stored in param to write the value in the NA_int_results.csv file
                endif
            endif
        enddo
    enddo

call scanfile (fnme, "logarithmic_velocity", p%logarithmic_velocity, p%logarithmic_velocity_desc, res, vocal, nd, range, par)

call scanfile (fnme, "debug", p%debug, p%debug_desc, res, vocal, nd, range, par)

call scanfile (fnme, "save_PTT_paths", p%save_PTT_paths, p%save_PTT_paths_desc, res, vocal, nd, range, par)

call scanfile (fnme, "save_cooling_rates", p%save_cooling_rates, p%save_cooling_rates_desc, res, vocal, nd, range, par) !Maxime

call scanfile (fnme, "save_eroded_volume", p%save_eroded_volume, p%save_eroded_volume_desc, res, vocal, nd, range, par)

call scanfile (fnme, "error_predictions", p%error_predictions, p%error_predictions_desc, res, vocal, nd, range, par)

call scanfile (fnme, "misfit_weight_AHE", p%misfit_weight_AHE, p%misfit_weight_AHE_desc, res, vocal, nd, range, par)

call scanfile (fnme, "misfit_weight_ZHE", p%misfit_weight_ZHE, p%misfit_weight_ZHE_desc, res, vocal, nd, range, par)

call scanfile (fnme, "misfit_weight_AFT", p%misfit_weight_AFT, p%misfit_weight_AFT_desc, res, vocal, nd, range, par)

call scanfile (fnme, "misfit_weight_ZFT", p%misfit_weight_ZFT, p%misfit_weight_ZFT_desc, res, vocal, nd, range, par)

call scanfile (fnme, "misfit_weight_BAR", p%misfit_weight_BAR, p%misfit_weight_BAR_desc, res, vocal, nd, range, par)

call scanfile (fnme, "misfit_weight_MAR", p%misfit_weight_MAR, p%misfit_weight_MAR_desc, res, vocal, nd, range, par)

call scanfile (fnme, "misfit_weight_HAR", p%misfit_weight_HAR, p%misfit_weight_HAR_desc, res, vocal, nd, range, par)

call scanfile (fnme, "misfit_weight_KAR", p%misfit_weight_KAR, p%misfit_weight_KAR_desc, res, vocal, nd, range, par)

call scanfile (fnme, "misfit_weight_FTLD", p%misfit_weight_FTLD, p%misfit_weight_FTLD_desc, res, vocal, nd, range, par)

call scanfile (fnme, "misfit_weight_TH", p%misfit_weight_TH, p%misfit_weight_TH_desc, res, vocal, nd, range, par)

call scanfile (fnme, "misfit_weight_43He", p%misfit_weight_43He, p%misfit_weight_43He_desc, res, vocal, nd, range, par)

call scanfile (fnme, "misfit_weight_TL", p%misfit_weight_TL, p%misfit_weight_TL_desc, res, vocal, nd, range, par)

call scanfile (fnme, "misfit_weight_OSL", p%misfit_weight_OSL, p%misfit_weight_OSL_desc, res, vocal, nd, range, par)

call scanfile (fnme, "score_with_data", p%score_with_data, p%score_with_data_desc, res, vocal, nd, range, par)

call scanfile (fnme, "inversion_mode", p%inversion_mode, p%inversion_mode_desc, res, vocal, nd, range, par)

call scanfile (fnme, "misfit_weight_ESR", p%misfit_weight_ESR, p%misfit_weight_ESR_desc,&
                res, vocal, nd, range, par)                                                                                                                 
call scanfile (fnme, "maximum_number_of_iterations", p%maximum_number_of_iterations, &
                      p%maximum_number_of_iterations_desc, res, vocal, nd, range, par)

call scanfile (fnme, "sample_size_for_first_iteration", p%sample_size_for_first_iteration, &
                      p%sample_size_for_first_iteration_desc, res, vocal, nd, range, par)

call scanfile (fnme, "sample_size_for_all_other_iterations", p%sample_size_for_all_other_iterations, &
                      p%sample_size_for_all_other_iterations_desc, res, vocal, nd, range, par)

call scanfile (fnme, "number_of_cells_to_resample", p%number_of_cells_to_resample, &
                      p%number_of_cells_to_resample_desc, res, vocal, nd, range, par)
call scanfile (fnme, "save_ages_inversion", p%save_ages_inversion, &
                      p%save_ages_inversion_desc, res, vocal, nd, range, par)
call scanfile (fnme, "TL_doser", p%TL_doser, p%TL_doser_desc, res, vocal, nd, range, par)

call scanfile (fnme, "TL_D0", p%TL_D0, p%TL_D0_desc, res, vocal, nd, range, par)

call scanfile (fnme, "TL_a", p%TL_a, p%TL_a_desc, res, vocal, nd, range, par)

call scanfile (fnme, "TL_b", p%TL_b, p%TL_b_desc, res, vocal, nd, range, par)

call scanfile (fnme, "TL_Et", p%TL_Et, p%TL_Et_desc, res, vocal, nd, range, par)

call scanfile (fnme, "TL_logs", p%TL_logs, p%TL_logs_desc, res, vocal, nd, range, par)

call scanfile (fnme, "TL_logrho", p%TL_logrho, p%TL_logrho_desc, res, vocal, nd, range, par)

call scanfile (fnme, "OSL_doser", p%OSL_doser, p%OSL_doser_desc, res, vocal, nd, range, par)

call scanfile (fnme, "OSL_D0", p%OSL_D0, p%OSL_D0_desc, res, vocal, nd, range, par)

call scanfile (fnme, "OSL_Eu", p%OSL_Eu, p%OSL_Eu_desc, res, vocal, nd, range, par)

call scanfile (fnme, "OSL_Et", p%OSL_Et, p%OSL_Et_desc, res, vocal, nd, range, par)

call scanfile (fnme, "OSL_logs", p%OSL_logs, p%OSL_logs_desc, res, vocal, nd, range, par)

call scanfile (fnme, "OSL_logrho", p%OSL_logrho, p%OSL_logrho_desc, res, vocal, nd, range, par)

call scanfile (fnme, "ESR_doser", p%ESR_doser, p%ESR_doser_desc, res, vocal, nd, range, par)

call scanfile (fnme, "ESR_D0", p%ESR_D0, p%ESR_D0_desc, res, vocal, nd, range, par)

call scanfile (fnme, "ESR_logs", p%ESR_logs, p%ESR_logs_desc, res, vocal, nd, range, par)

call scanfile (fnme, "ESR_Et", p%ESR_Et, p%ESR_Et_desc, res, vocal, nd, range, par)

call scanfile (fnme, "ESR_sigmaEt", p%ESR_sigmaEt, p%ESR_sigmaEt_desc, res, vocal, nd, range, par)

call scanfile (fnme, "ESR_a", p%ESR_a, p%ESR_a_desc, res, vocal, nd, range, par)            

call scanfile (fnme, "ESR_b", p%ESR_b, p%ESR_b_desc, res, vocal, nd, range, par)            
return
end

!--------------------------------------------------------------------------------------------

subroutine rescale_parameter (max_val,range,nd, new_par_val,ref)

! This routine give a new value for the parameter if its value is 
! above another parameter value set by the user (for inversion)
! ref is whether the user wish the maximum value (ref = 0) or
! the minum value (ref = 1)
integer nd, ref
double precision new_par_val,max_val,u
real*4,intent(inout)::range(2,1024)

call RANDOM_NUMBER(u)

if (ref.eq.0) then
    new_par_val = range(1,nd) + (max_val-range(1,nd))*u
elseif (ref.eq.1) then
    new_par_val = max_val + (range(2,nd) - max_val)*u
endif
! print *, 'u = ', u, max_val, range(1,nd), range(2,nd)

end subroutine rescale_parameter


