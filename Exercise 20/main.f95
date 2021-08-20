!=============================================
! Exercise 20: Geostrophic adjustment in 3d
!=============================================
! Author: J. Kaempf, November 2009 

PROGRAM full3D

USE param
USE sub

! local parameters
INTEGER :: ntot, nout
REAL :: xx, yy, rad

!**********
CALL INIT  ! initialisation
!**********

! runtime parameters
ntot = INT(60.*60.*60./dt)
time = 0.0

! output frequency
nout = INT(60.*60./dt)

! open data output files
OPEN(10,file ='rhoS.dat',form='formatted')
OPEN(20,file ='uS.dat',form='formatted')
OPEN(30,file ='vS.dat',form='formatted')
OPEN(40,file ='eta.dat',form='formatted')

OPEN(11,file ='rhoV.dat',form='formatted')
OPEN(21,file ='uV.dat',form='formatted')
OPEN(31,file ='vV.dat',form='formatted')

! output of initial fields
DO j = 1,ny
  WRITE(10,'(25F12.6)')(rho(1,j,k),k=1,nx)
  WRITE(20,'(25F12.6)')(u(1,j,k),k=1,nx)
  WRITE(30,'(25F12.6)')(v(1,j,k),k=1,nx)
  WRITE(40,'(25F12.6)')(q(0,j,k)/(RHOREF*G),k=1,nx)
END DO

DO i = 1,nz
  WRITE(11,'(25F12.6)')(rho(i,ny/2,k),k=1,nx)
  WRITE(21,'(25F12.6)')(u(i,ny/2,k),k=1,nx)
  WRITE(31,'(25F12.6)')(v(i,ny/2,k),k=1,nx)
END DO

!---------------------------
! simulation loop
!---------------------------
DO n = 1,ntot

time = time + dt
write(6,*)"time (days)", time/(24.*3600.)
ad = MIN(time/(6.0*3600.0),1.0) ! adjustment

! prescription of low-density patch 
IF(ad<1.0)THEN

  DO k = 1,nx
  DO j = 1,ny

  xx = REAL(k-13)*dx
  yy = REAL(j-13)*dy
  rad = SQRT(xx*xx+yy*yy)

  IF(rad < 5000.0)THEN
     DO i = 1,10
       rho(i,j,k) = rho(i,j,k) - dt*0.1/(6*3600.)
     END DO
  END IF
  END DO
  END DO

ENDIF

! prognostic equations
CALL dyn

! data output
IF(MOD(n,nout)==0)THEN
DO j = 1,ny
    WRITE(10,'(25F12.6)')(rho(1,j,k),k=1,nx)
    WRITE(20,'(25F12.6)')(u(1,j,k),k=1,nx)
    WRITE(30,'(25F12.6)')(v(1,j,k),k=1,nx)
    WRITE(40,'(25F12.6)')(q(0,j,k)/(RHOREF*G),k=1,nx)
END DO

DO i = 1,nz
    WRITE(11,'(25F12.6)')(rho(i,ny/2,k),k=1,nx)
    WRITE(21,'(25F12.6)')(u(i,ny/2,k),k=1,nx)
    WRITE(31,'(25F12.6)')(v(i,ny/2,k),k=1,nx)
END DO
  WRITE(6,*)"Data output at time = ",time/(24.*3600.)
ENDIF

END DO ! end of iteration loop

END PROGRAM full3D
