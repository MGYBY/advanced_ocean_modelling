!==========================================
! Exercise 6: Kelvin-Helmholtz instability
!==========================================
! Author: J. Kaempf, August 2009 

PROGRAM slice

USE param
USE sub

! local parameters
INTEGER :: ntot, nout
REAL :: wl, ps

!**********
CALL INIT  ! initialisation
!**********

! runtime parameters
ntot = INT(100.*60./dt)
time = 0.0

! output frequency
nout = INT(30./dt)

! open files for output
OPEN(10,file ='q.dat',form='formatted')
OPEN(20,file ='u.dat',form='formatted')
OPEN(30,file ='w.dat',form='formatted')
OPEN(40,file ='eta.dat',form='formatted')
OPEN(50,file ='rho.dat',form='formatted')

!---------------------------
! simulation loop
!---------------------------
DO n = 1,ntot

time = time + dt
write(6,*)"time (hours)", time/(3600.)

! prognostic equations
CALL dyn

! data output (start after 20 minutes)
IF(n > ntot/5)then
IF(MOD(n,nout)==0)THEN
  DO i = 1,nz
    WRITE(10,'(101F12.6)')(q(i,k)/(RHOREF*G),k=1,nx)
    WRITE(20,'(101F12.6)')(u(i,k),k=1,nx)
    WRITE(30,'(101F12.6)')(w(i,k),k=1,nx)
    WRITE(50,'(101F12.6)')(rho(i,k)-RHOREF,k=1,nx)
  END DO
  WRITE(40,'(101F12.6)')(q(0,k)/(RHOREF*G),k=1,nx)
  WRITE(6,*)"Data output at time = ",time/(24.*3600.)
ENDIF
ENDIF

END DO ! end of iteration loop

END PROGRAM slice
