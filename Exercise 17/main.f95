!==================================
! Exercise 17: Tidal-mixing fronts
!==================================
! Author: J. Kaempf, August 2009 

PROGRAM slice

USE param
USE sub

! local parameters
INTEGER :: ntot, nout
REAL :: azamb(0:nx+1) ! eddy

!**********
CALL INIT  ! initialisation
!**********

DO k = 0,60
 azamb(k) = 0.0
END DO
DO k = 61,70
 azamb(k) = 0.05*REAL(k-60)/10.
END DO
DO k = 71,nx+1 
 azamb(k) = azamb(70)
END DO

! runtime parameters
ntot = INT(24.*60*60./dt)
time = 0.0

! output parameter
nout = INT(15.*60./dt)

OPEN(10,file ='q.dat',form='formatted')
OPEN(20,file ='u.dat',form='formatted')
OPEN(25,file ='v.dat',form='formatted')
OPEN(30,file ='w.dat',form='formatted')
OPEN(40,file ='eta.dat',form='formatted')
OPEN(50,file ='rho.dat',form='formatted')
OPEN(60,file ='c1.dat',form='formatted')
OPEN(61,file ='c2.dat',form='formatted')
OPEN(62,file ='c3.dat',form='formatted')


! data output of initial distributions
DO i = 1,nz
  WRITE(10,'(101F12.6)')(q(i,k)/(RHOREF*G),k=1,nx)
  WRITE(20,'(101F12.6)')(u(i,k),k=1,nx)
  WRITE(25,'(101F12.6)')(v(i,k),k=1,nx)
  WRITE(30,'(101F12.6)')(w(i,k),k=1,nx)
  WRITE(50,'(101F12.6)')(rho(i,k),k=1,nx)
  WRITE(60,'(101F12.6)')(c1(i,k),k=1,nx)
  WRITE(61,'(101F12.6)')(c2(i,k),k=1,nx)
  WRITE(62,'(101F12.6)')(c3(i,k),k=1,nx)
END DO
WRITE(40,'(101F12.6)')(q(0,k)/(RHOREF*G),k=1,nx)

!---------------------------
! simulation loop
!---------------------------
DO n = 1,ntot

time = time + dt
write(6,*)"time (hours)", time/(3600.)

! temporal variation of turbulence fields
ad = ABS(SIN(time/(12.*3600.0)*2.*PI))

DO i = 0,nz+1
DO k = 1,nx+1
 az(i,k) = azamb(k)*ad + azmin
 kz(i,k) = az(i,k)
END DO
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
    WRITE(60,'(101F12.6)')(c1(i,k),k=1,nx)
    WRITE(61,'(101F12.6)')(c2(i,k),k=1,nx)
    WRITE(62,'(101F12.6)')(c3(i,k),k=1,nx)
  END DO
  WRITE(40,'(101F12.6)')(q(0,k)/(RHOREF*G),k=1,nx)
  WRITE(6,*)"Data output at time = ",time/(24.*3600.)
ENDIF

END DO ! end of iteration loop

END PROGRAM slice
