double precision function heat_production (x0,y0,z0,temp,time)

! this file is used to define the heat production in (°C/Myr) as a funciton of initial position (x0 and y0 in km),
! initial depth (z0 in km), current temperature (in °C) and time (in Myr) in case the heat production in the input file 
! (topo_parameters.txt) is negative

implicit none

double precision x0,y0,z0,temp,time

heat_production=3.d0

return
end
