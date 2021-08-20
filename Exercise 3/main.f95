!=========================================
! Exercise 3: Short Surface Gravity Waves
!=========================================
! Author: J. Kaempf, August 2009 

PROGRAM slice

USE param
USE sub

! local parameters
INTEGER :: n, ntot, nout
REAL :: wl, ps

period = 8.0 ! forcing period in seconds
amplitude = 1.0 ! forcing amplitude

wl = G*period*period/(2.*PI)
write(6,*)"deep-water wavelength (m) is ",wl
ps = wl/period
write(6,*)"deep-water phase speed (m/s) is ",ps

pause

!**********
CALL INIT  ! initialisation
!**********

! runtime parameters
ntot = INT(100./dt)
time = 0.0

! output parameter
nout = INT(1./dt)

! open files for data output
OPEN(10,file ='dp.dat',form='formatted')
OPEN(20,file ='u.dat',form='formatted')
OPEN(30,file ='w.dat',form='formatted')
OPEN(40,file ='eta.dat',form='formatted')

!---------------------------
! simulation loop
!---------------------------
DO n = 1,ntot

time = time + dt
write(6,*)"time (hours)", time/(3600.)

! variation of forcing
ad = amplitude*SIN(2.*PI*time/period)

! prognostic equations
CALL dyn

! data output
IF(MOD(n,nout)==0)THEN
  DO i = 1,nz
    WRITE(10,'(101F12.6)')(dp(i,k)/(RHO*G),k=1,nx)
    WRITE(20,'(101F12.6)')(u(i,k),k=1,nx)
    WRITE(30,'(101F12.6)')(w(i,k),k=1,nx)
  END DO
  WRITE(40,'(101F12.6)')(dp(0,k)/(RHO*G),k=1,nx)
  WRITE(6,*)"Data output at time = ",time/(24.*3600.)

ENDIF

END DO ! end of iteration loop

END PROGRAM slice
