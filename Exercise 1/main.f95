!=======================================
! Exercise 1:  The surface Ekman layer 
!=======================================
! Author: J. Kaempf, August 2009

PROGRAM ekman
USE param
USE sub

INTEGER :: n,ntot ! time level and output parameters
REAL :: locdep ! distance from sea surface

CALL init

! various eddy viscosity settings
 mode = 1 ! constant
! mode = 2 ! intermediate minimum
! mode = 3 ! Prandtl mixing-length closure

! runtime parameters
ntot = 5.*24.*3600./dt ! runtime is 5 days
time = 0.0

! open files for output
IF(mode==1)OPEN(10,file ='uvprof1.dat',form='formatted')
IF(mode==2)OPEN(10,file ='uvprof2.dat',form='formatted')
IF(mode==3)OPEN(10,file ='uvprof3.dat',form='formatted')

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

! data output (steady state)

DO i = 1,nz
  locdep = -REAL(i)*DZ+0.5*DZ
  WRITE(10,*)locdep,u(i),v(i)
END DO

END PROGRAM ekman
