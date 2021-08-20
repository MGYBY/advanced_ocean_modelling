!===========================================
! Exercise 13: Strfatified flows on a slope
!===========================================
! Author: J. Kaempf, August 2009 

PROGRAM slice

USE param
USE sub
USE random

! local parameters
INTEGER :: ntot, nout

!**********
CALL INIT  ! initialisation
!**********

! runtime parameters
ntot = INT(120.*60./dt)
time = 0.0

! output frequency
nout = INT(60./dt)

OPEN(10,file ='q.dat',form='formatted')
OPEN(20,file ='u.dat',form='formatted')
OPEN(30,file ='w.dat',form='formatted')
OPEN(40,file ='eta.dat',form='formatted')
OPEN(50,file ='rho.dat',form='formatted')

!initialisation of random function
ist = -1
randm = ran3(ist)

DO k = 0,nx+1
DO i = 41,nz+1
  rho(i,k) = 0.2 + ran3(ist)*0.00001
END DO
END DO

! bottom inclination (5 degrees converted to radians)
alpha = 5.0*PI/180.

!---------------------------
! simulation loop
!---------------------------
DO n = 1,ntot

time = time + dt
write(6,*)"time (hours)", time/(3600.)

! prognostic equations
CALL dyn

! data output
IF(MOD(n,nout)==0)THEN
  DO i = 1,nz
    WRITE(10,'(101F12.6)')(q(i,k)/(RHOREF*G),k=1,nx)
    WRITE(20,'(101F12.6)')(u(i,k),k=1,nx)
    WRITE(30,'(101F12.6)')(w(i,k),k=1,nx)
    WRITE(50,'(101F12.6)')(rho(i,k),k=1,nx)
  END DO
  WRITE(40,'(101F12.6)')(q(0,k)/(RHOREF*G),k=1,nx)
  WRITE(6,*)"Data output at time = ",time/(24.*3600.)
ENDIF
!ENDIF

END DO ! end of iteration loop

END PROGRAM slice
