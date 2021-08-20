!==============================================
! Exercise 10: Slope Convection near the Shore
!==============================================
! Author: J. Kaempf, August 2009 

PROGRAM slice

USE param
USE sub

! local parameters
INTEGER :: ntot, nout

!**********
CALL INIT  ! initialisation
!**********

! runtime parameters
ntot = INT(12.*3600./dt)
time = 0.0

! output frequency
nout = INT(6.*60./dt)

! open files for data output
OPEN(10,file ='q.dat',form='formatted')
OPEN(20,file ='u.dat',form='formatted')
OPEN(30,file ='w.dat',form='formatted')
OPEN(40,file ='eta.dat',form='formatted')
OPEN(50,file ='rho.dat',form='formatted')
OPEN(60,file ='c.dat',form='formatted')

OPEN(70,file ='TRx.dat',form='formatted',recl = 10000000)
OPEN(80,file ='TRz.dat',form='formatted',recl = 10000000)

!---------------------------
! simulation loop
!---------------------------
DO n = 1,ntot

time = time + dt
write(6,*)"time (hours)", time/(3600.)

! prognostic equations
CALL dyn

! float prediction scheme
DO ii = 1,ntra

! locate grid cell of tracer
ipos = INT(tra(ii,1)/dz)+1
kpos = INT(tra(ii,2)/dx+0.5)+1
uu = 0.5*(u(ipos,kpos)+u(ipos,kpos-1))
ww = 0.5*(w(ipos,kpos)+w(ipos+1,kpos))

! avoid stranding
IF( .not.wet(ipos,kpos).AND.uu>0 ) uu = 0.0
IF( .not.wet(ipos,kpos-1).AND.uu<0 ) uu = 0.0
IF( .not.wet(ipos+1,kpos).AND.ww<0 ) ww = 0.0
IF( ipos==1.AND.ww>0 ) ww = 0.0
IF( ipos==nz.AND.ww<0 ) ww = 0.0

! change of location
tra(ii,1) = tra(ii,1)-dt*ww
tra(ii,2) = tra(ii,2)+dt*uu

tra(ii,1) = MAX(tra(ii,1),0.0)

END DO

! data output
IF(MOD(n,nout)==0)THEN
  DO i = 1,nz
    WRITE(10,'(201F12.6)')(q(i,k)/(RHOREF*G),k=1,nx)
    WRITE(20,'(201F12.6)')(u(i,k),k=1,nx)
    WRITE(30,'(201F12.6)')(w(i,k),k=1,nx)
    WRITE(50,'(201F12.6)')(rho(i,k),k=1,nx)
    WRITE(60,'(201F12.6)')(c(i,k),k=1,nx)
  END DO
    WRITE(70,'(3000F12.6)')(tra(ii,2),ii=1,ntra)
    WRITE(80,'(3000F12.6)')(tra(ii,1),ii=1,ntra)

  WRITE(40,'(201F12.6)')(q(0,k)/(RHOREF*G),k=1,nx)
  WRITE(6,*)"Data output at time = ",time/(24.*3600.)
ENDIF


END DO ! end of iteration loop

END PROGRAM slice
