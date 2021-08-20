!=========================================
! Exercise 21: Eddy formation in a Strait
!=========================================
! Author: J. Kaempf, August 2009 

PROGRAM full3D

USE param
USE sub

! local parameters
INTEGER :: ntot, nout

!**********
CALL INIT  ! initialisation
!**********

! runtime parameters 
ntot = INT(10.*24.*60.*60./dt)
time = 0.0

! output frequency
nout = INT(3.*60.*60./dt)

! open files and output of initial distributions
OPEN(10,file ='rhoS.dat',form='formatted')
OPEN(20,file ='uS.dat',form='formatted')
OPEN(30,file ='vS.dat',form='formatted')
OPEN(40,file ='eta.dat',form='formatted')
OPEN(50,file ='cS.dat',form='formatted')

DO j = 1,ny
  WRITE(10,'(100F12.6)')(rho(1,j,k),k=1,nx)
  WRITE(20,'(100F12.6)')(u(1,j,k),k=1,nx)
  WRITE(30,'(100F12.6)')(v(1,j,k),k=1,nx)
  WRITE(50,'(100F12.6)')(c(1,j,k),k=1,nx)
  WRITE(40,'(100F12.6)')(q(0,j,k)/(RHOREF*G),k=1,nx)
END DO

!---------------------------
! simulation loop
!---------------------------
DO n = 1,ntot

time = time + dt
write(6,*)"time (days)", time/(24.*3600.)
ad = MIN(time/(1.*24.0*3600.0),1.0)

! prognostic equations
CALL dyn

! data output
IF(MOD(n,nout)==0)THEN
DO j = 1,ny
  WRITE(10,'(100F12.6)')(rho(1,j,k),k=1,nx)
  WRITE(20,'(100F12.6)')(u(1,j,k),k=1,nx)
  WRITE(30,'(100F12.6)')(v(1,j,k),k=1,nx)
  WRITE(40,'(100F12.6)')(q(0,j,k)/(RHOREF*G),k=1,nx)
  WRITE(50,'(100F12.6)')(c(1,j,k),k=1,nx)
END DO
  WRITE(6,*)"Data output at time = ",time/(24.*3600.)
ENDIF

END DO ! end of iteration loop

END PROGRAM full3D
