!=============================================
! Exercise 22: Exchange flow through a strait
!=============================================
! Author: J. Kaempf, August 2009 

PROGRAM full3D

USE param
USE sub

! local parameters
INTEGER :: ntot, nout
REAL :: ttt, xx, yy, rad

!**********
CALL INIT  ! initialisation
!**********

! runtime parameters
ntot = INT(20.*24.*60.*60./dt)
time = 0.0

! output frequency
nout = INT(6.*60.*60./dt)

! open files for output
OPEN(10,file ='rhoS.dat',form='formatted')
OPEN(20,file ='uS.dat',form='formatted')
OPEN(30,file ='vS.dat',form='formatted')
OPEN(40,file ='eta.dat',form='formatted')
OPEN(50,file ='cS.dat',form='formatted')

OPEN(11,file ='rhoV.dat',form='formatted')
OPEN(21,file ='uV.dat',form='formatted')
OPEN(31,file ='vV.dat',form='formatted')
OPEN(51,file ='cV.dat',form='formatted')

OPEN(12,file ='rhoB.dat',form='formatted')
OPEN(22,file ='uB.dat',form='formatted')
OPEN(32,file ='vB.dat',form='formatted')
OPEN(52,file ='cB.dat',form='formatted')

! output of initial distributions
! surface distributions
DO j = 1,ny
  WRITE(10,'(101F12.6)')(rho(1,j,k),k=1,nx)
  WRITE(20,'(101F12.6)')(u(1,j,k),k=1,nx)
  WRITE(30,'(101F12.6)')(v(1,j,k),k=1,nx)
  WRITE(50,'(101F12.6)')(c(1,j,k),k=1,nx)
  WRITE(40,'(101F12.6)')(q(0,j,k)/(RHOREF*G),k=1,nx)
END DO

! vertical transects
DO i = 1,nz
  WRITE(11,'(51F12.6)')(rho(i,j,51),j=1,ny)
  WRITE(21,'(51F12.6)')(u(i,j,51),j=1,ny)
  WRITE(31,'(51F12.6)')(v(i,j,51),j=1,ny)
  WRITE(51,'(51F12.6)')(c(i,j,51),j=1,ny)
END DO

! bottom distributions
DO j = 1,ny
  WRITE(12,'(101F12.6)')(rho(ib(j,k),j,k),k=1,nx)
  WRITE(22,'(101F12.6)')(u(ib(j,k),j,k),k=1,nx)
  WRITE(32,'(101F12.6)')(v(ib(j,k),j,k),k=1,nx)
  WRITE(52,'(101F12.6)')(c(ib(j,k),j,k),k=1,nx)
END DO

!---------------------------
! simulation loop
!---------------------------
DO n = 1,ntot

time = time + dt
write(6,*)"time (days)", time/(24.*3600.)
ad = MIN(time/(2.*24.0*3600.0),1.0) ! initial adjustment

IF(ad<1.0)THEN

  DO k = 1,51
  DO j = 1,ny
  DO i = 1,nz
   IF(wet(i,j,k)) rho(i,j,k) = rho(i,j,k) + dt*2.0/(2.*24*3600.)
  END DO
  END DO
  END DO

ENDIF

! prognostic equations
CALL dyn

! data output
IF(MOD(n,nout)==0)THEN
DO j = 1,ny
    WRITE(10,'(101F12.6)')(rho(1,j,k),k=1,nx)
    WRITE(20,'(101F12.6)')(u(1,j,k),k=1,nx)
    WRITE(30,'(101F12.6)')(v(1,j,k),k=1,nx)
    WRITE(40,'(101F12.6)')(q(0,j,k)/(RHOREF*G),k=1,nx)
    WRITE(50,'(101F12.6)')(c(1,j,k),k=1,nx)
END DO

DO j = 1,ny
    WRITE(12,'(101F12.6)')(rho(ib(j,k),j,k),k=1,nx)
    WRITE(22,'(101F12.6)')(u(ib(j,k),j,k),k=1,nx)
    WRITE(32,'(101F12.6)')(v(ib(j,k),j,k),k=1,nx)
    WRITE(52,'(101F12.6)')(c(ib(j,k),j,k),k=1,nx)
END DO

DO i = 1,nz
    WRITE(11,'(51F12.6)')(rho(i,j,51),j=1,ny)
    WRITE(21,'(51F12.6)')(u(i,j,51),j=1,ny)
    WRITE(31,'(51F12.6)')(v(i,j,51),j=1,ny)
    WRITE(51,'(51F12.6)')(c(i,j,51),j=1,ny)
END DO
  WRITE(6,*)"Data output at time = ",time/(24.*3600.)
ENDIF

END DO ! end of iteration loop

END PROGRAM full3D
