!=============================
! Exercise 19: Ekman pumping
!=============================
! Author: J. Kaempf, August 2009 

PROGRAM slice

USE param
USE sub

! local parameters
INTEGER :: ntot, nout
REAL :: ttt

!**********
CALL INIT  ! initialisation
!**********

! runtime parameters
ntot = INT(50.*24.*60.*60./dt)
time = 0.0

! output frequency
nout = INT(12.*60.*60./dt)

! open files for output
OPEN(10,file ='q.dat',form='formatted')
OPEN(20,file ='u.dat',form='formatted')
OPEN(25,file ='v.dat',form='formatted')
OPEN(30,file ='w.dat',form='formatted')
OPEN(40,file ='eta.dat',form='formatted')
OPEN(50,file ='rho.dat',form='formatted')
OPEN(60,file ='tauy.dat',form='formatted')

! output of initial distributions
DO i = 1,nz
  WRITE(10,'(101F12.6)')(q(i,k)/(RHOREF*G),k=1,nx)
  WRITE(20,'(101F12.6)')(u(i,k),k=1,nx)
  WRITE(25,'(101F12.6)')(v(i,k),k=1,nx)
  WRITE(30,'(101F12.6)')(w(i,k),k=1,nx)
  WRITE(50,'(101F12.6)')(rho(i,k),k=1,nx)
END DO

WRITE(40,'(101F12.6)')(q(0,k)/(RHOREF*G),k=1,nx)

! wind-stress forcing
DO k = 1,nx
  tauy(k) = -0.1*COS(REAL(k-1)/REAL(nx-1)*PI)
END DO

WRITE(60,'(101F12.6)')(tauy(k),k=1,nx)

!---------------------------
! simulation loop
!---------------------------
DO n = 1,ntot

time = time + dt
write(6,*)"time (days)", time/(24.*3600.)
ad = MIN(time/(2.*24.*3600.0),1.0)

! adjustment of wind-stress forcing
DO k = 1,nx
  tauy(k) = -ad*0.1*COS(REAL(k-1)/REAL(nx-1)*PI)
END DO

! prognostic equations
CALL dyn

! data output
IF(MOD(n,nout)==0)THEN
  DO i = 1,nz
    WRITE(10,'(101F12.6)')(q(i,k)/(RHOREF*G),k=1,nx)
    WRITE(20,'(101F12.6)')(u(i,k),k=1,nx)
    WRITE(25,'(101F12.6)')(v(i,k),k=1,nx)
    WRITE(30,'(101F12.6)')(w(i,k),k=1,nx)
    WRITE(50,'(101F12.6)')(rho(i,k),k=1,nx)
  END DO
  WRITE(40,'(101F12.6)')(q(0,k)/(RHOREF*G),k=1,nx)
  WRITE(6,*)"Data output at time = ",time/(24.*3600.)
ENDIF

END DO ! end of iteration loop

END PROGRAM slice
