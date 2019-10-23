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
integer :: vocal
type (parameters) p
integer nd
real*4 range(2,nd), par(nd)

integer :: res, i, j, npoint_max, nstep_max
character*10 ci, cj

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

allocate (p%time_topo(p%ntime+1), p%amplification(p%ntime+1), p%offset(p%ntime+1), p%output(p%ntime+1))
p%time_topo(:) = 0.d0
p%amplification(:) = 1.d0
p%offset(:) = 0.d0
p%output(:) = 1

  do i = 1, p%ntime + 1
  write (ci,'(i10)') i
  j = i
  call scanfile (fnme, "time_topo"//trim(adjustl(ci)), p%time_topo(j), p%time_topo_desc, res, vocal, nd, range, par)
  call scanfile (fnme, "amplification"//trim(adjustl(ci)), p%amplification(j), p%amplification_desc, res, vocal, nd, range, par)
  if (res.eq.999) p%amplification(j)=p%amplification(i)
  call scanfile (fnme, "offset"//trim(adjustl(ci)), p%offset(j), p%offset_desc, res, vocal, nd, range, par)
  if (res.eq.999) p%offset(j)=p%offset(i)
  call scanfile (fnme, "output"//trim(adjustl(ci)), p%output(j), p%output_desc, res, vocal, nd, range, par)
  enddo

!do i = 1, p%ntime
!if (p%time_topo(i).eq.0.d0.and.vocal.ne.4) stop 'There must be as many time_topo as the number of time steps'
!enddo

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

call scanfile (fnme, "sea_level_temperature", p%sea_level_temperature, p%sea_level_temperature_desc, res, vocal, nd, range, par)

call scanfile (fnme, "lapse_rate", p%lapse_rate, p%lapse_rate_desc, res, vocal, nd, range, par)

call scanfile (fnme, "heat_production", p%heat_production, p%heat_production_desc, res, vocal, nd, range, par)

call scanfile (fnme, "data_folder", p%data_folder, p%data_folder_desc, res, vocal, nd, range, par)

call scanfile (fnme, "default_age", p%default_age, p%default_age_desc, res, vocal, nd, range, par)

call scanfile (fnme, "FT_code_flag", p%FT_code_flag, p%FT_code_flag_desc, res, vocal, nd, range, par)

call scanfile (fnme, "misfit_slope", p%misfit_slope, p%misfit_slope_desc, res, vocal, nd, range, par)

call scanfile (fnme, "misfit_corrected", p%misfit_corrected, p%misfit_corrected_desc, res, vocal, nd, range, par)

call scanfile (fnme, "fault_advect_flag", p%fault_advect_flag, p%fault_advect_flag_desc, res, vocal, nd, range, par)

call scanfile (fnme, "shear_heating", p%shear_heating, p%shear_heating_desc, res, vocal, nd, range, par)

call scanfile (fnme, "age_TL_flag", p%age_TL_flag, p%age_TL_flag_desc, res, vocal, nd, range, par)

call scanfile (fnme, "age_TL_flag", p%age_TL_flag, p%age_TL_flag_desc, res, vocal, nd, range, par)

call scanfile (fnme, "age_AHe_flag", p%age_AHe_flag, p%age_AHe_flag_desc, res, vocal, nd, range, par)

call scanfile (fnme, "age_ZHe_flag", p%age_ZHe_flag, p%age_ZHe_flag_desc, res, vocal, nd, range, par)

call scanfile (fnme, "age_AFT_flag", p%age_AFT_flag, p%age_AFT_flag_desc, res, vocal, nd, range, par)

call scanfile (fnme, "age_ZFT_flag", p%age_ZFT_flag, p%age_ZFT_flag_desc, res, vocal, nd, range, par)

call scanfile (fnme, "age_FTL_flag", p%age_FTL_flag, p%age_FTL_flag_desc, res, vocal, nd, range, par)
if (p%age_FTL_flag.eq.1) p%age_AFT_flag = 1

call scanfile (fnme, "age_KAr_flag", p%age_KAr_flag, p%age_KAr_flag_desc, res, vocal, nd, range, par)

call scanfile (fnme, "age_BAr_flag", p%age_BAr_flag, p%age_BAr_flag_desc, res, vocal, nd, range, par)

call scanfile (fnme, "age_MAr_flag", p%age_MAr_flag, p%age_MAr_flag_desc, res, vocal, nd, range, par)

call scanfile (fnme, "age_HAr_flag", p%age_HAr_flag, p%age_HAr_flag_desc, res, vocal, nd, range, par)


call scanfile (fnme, "nfault", p%nfault, p%nfault_desc, res, vocal, nd, range, par)
if (vocal.eq.4) p%nfault = 1

call scanfile (fnme, "lon1", p%x1, p%x1_desc, res, vocal, nd, range, par)

call scanfile (fnme, "lat1", p%y1, p%y1_desc, res, vocal, nd, range, par)

call scanfile (fnme, "lon2", p%x2, p%x2_desc, res, vocal, nd, range, par)

call scanfile (fnme, "lat2", p%y2, p%y2_desc, res, vocal, nd, range, par)

allocate (p%npoint(p%nfault), p%nstep(p%nfault))
p%npoint = 0
p%nstep = 0
  do i = 1, p%nfault
  write (ci,'(i10)') i
  call scanfile (fnme, "npoint"//trim(adjustl(ci)), p%npoint(i), p%npoint_desc, res, vocal, nd, range, par)
  call scanfile (fnme, "nstep"//trim(adjustl(ci)), p%nstep(i), p%nstep_desc, res, vocal, nd, range, par)
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
      do j = 1, 4
      write (cj,'(i10)') j
      call scanfile (fnme, "r"//trim(adjustl(ci))//"_"//trim(adjustl(cj)), p%r(j,i), p%r_desc, res, vocal, nd, range, par)
      if (res.eq.999) p%r(j,i)=p%r(j-1,i)
      enddo
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
  do i = 1, p%nfault
  write (ci,'(i10)') i
    do j = 1, p%nstep(i)
    write (cj,'(i10)') j
    call scanfile (fnme, "time_start"//trim(adjustl(ci))//"_"//trim(adjustl(cj)), p%time_start(j,i), p%time_start_desc, &
                   res, vocal, nd, range, par)
    if (res.eq.999) p%time_start(j,i)=p%time_end(j-1,i)
    if (p%time_start(j,i).lt.0.d0) p%time_start(j,i)=p%time_topo(int(-p%time_start(j,i)))
    if (p%time_start(j,i).lt.0.d0) print*,int(-p%time_start(j,i))
    call scanfile (fnme, "time_end"//trim(adjustl(ci))//"_"//trim(adjustl(cj)), p%time_end(j,i), p%time_end_desc, &
                   res, vocal, nd, range, par)
    if (p%time_end(j,i).lt.0.d0) p%time_end(j,i)=p%time_topo(int(-p%time_end(j,i)))
    if (p%time_end(j,i).lt.0.d0) print*,int(-p%time_end(j,i))
    call scanfile (fnme, "velo"//trim(adjustl(ci))//"_"//trim(adjustl(cj)), p%velo(j,i), p%velo_desc, &
                  res, vocal, nd, range, par)
    enddo
  enddo

call scanfile (fnme, "logarithmic_velocity", p%logarithmic_velocity, p%logarithmic_velocity_desc, res, vocal, nd, range, par)

call scanfile (fnme, "debug", p%debug, p%debug_desc, res, vocal, nd, range, par)

call scanfile (fnme, "save_PTT_paths", p%save_PTT_paths, p%save_PTT_paths_desc, res, vocal, nd, range, par)

call scanfile (fnme, "save_eroded_volume", p%save_eroded_volume, p%save_eroded_volume_desc, res, vocal, nd, range, par)

call scanfile (fnme, "misfit_weight_age", p%misfit_weight_age, p%misfit_weight_age_desc, res, vocal, nd, range, par)

call scanfile (fnme, "misfit_weight_FTLD", p%misfit_weight_FTLD, p%misfit_weight_FTLD_desc, res, vocal, nd, range, par)

call scanfile (fnme, "misfit_weight_TH", p%misfit_weight_TH, p%misfit_weight_TH_desc, res, vocal, nd, range, par)

call scanfile (fnme, "misfit_weight_43He", p%misfit_weight_43He, p%misfit_weight_43He_desc, res, vocal, nd, range, par)

call scanfile (fnme, "misfit_weight_TL", p%misfit_weight_TL, p%misfit_weight_TL_desc, res, vocal, nd, range, par)

call scanfile (fnme, "maximum_number_of_iterations", p%maximum_number_of_iterations, &
                      p%maximum_number_of_iterations_desc, res, vocal, nd, range, par)

call scanfile (fnme, "sample_size_for_first_iteration", p%sample_size_for_first_iteration, &
                      p%sample_size_for_first_iteration_desc, res, vocal, nd, range, par)

call scanfile (fnme, "sample_size_for_all_other_iterations", p%sample_size_for_all_other_iterations, &
                      p%sample_size_for_all_other_iterations_desc, res, vocal, nd, range, par)

call scanfile (fnme, "number_of_cells_to_resample", p%number_of_cells_to_resample, &
                      p%misfit_weight_TL_desc, res, vocal, nd, range, par)

return
end

!--------------------------------------------------------------------------------------------

