double precision function thermal_diffusivity (x0,y0,z0,temp,time)

! this file is used to define the thermal diffusivity in (km^2/Myr) as a function of initial position (x0 and y0 in km),
! depth (z0 in km), current temperature (in °C) and time (in Myr) in case the diffusivity in the input file is negative

implicit none

double precision x0,y0,z0,temp,time

thermal_diffusivity=25.d0

return
end
