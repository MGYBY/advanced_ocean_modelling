!=======================================
! Exercise 2:  The bottom Ekman layer 
!=======================================
! Author: J. Kaempf, August 2009

PROGRAM ekman
USE param
USE sub

INTEGER :: n,ntot ! time level and output parameters
REAL :: locdep ! vertical location

CALL init

 mode = 1 ! constant eddy viscosity

! runtime parameters
ntot = 5.*24.*3600./dt ! runtime is 5 days
time = 0.0

! open files for output
OPEN(10,file ='uvprof1.dat',form='formatted')

!---------------------------
! simulation loop
!---------------------------

DO n = 1,ntot

time = time + dt

write(6,*)"Time (days) = ", time/(24.*3600)

! calculate eddy viscosity
CALL eddy

! predictor step for currents
CALL dyn

END DO

! data output

DO i = 1,nz
  locdep = -REAL(i)*DZ+0.5*DZ
  WRITE(10,*)locdep,u(i)+ugeo,v(i)+vgeo
END DO

END PROGRAM ekman
