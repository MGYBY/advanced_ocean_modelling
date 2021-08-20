!=============================================
! Exercise 25: Simulation of an El-Nino Event
!=============================================
!* Author: J. Kaempf, August 2009 

PROGRAM full3D

USE param
USE sub
USE random

! local parameters
INTEGER :: ntot, nout, ist, count
REAL :: ttt, xx, yy, rad, dist

!**********
CALL INIT  ! initialisation
!**********

! wind-stress forcing
DO j = 0,ny+1
DO k = 0,nx+1
  taux(j,k) = 0.0
  tauy(j,k) = 0.0
END DO
END DO

DO j = 51-10,51+10
  dist = REAL(j-51)
  DO k = 0,10
    taux(j,k) = 0.1*COS(2.*PI*dist/40.)
    tauy(j,k) = 0.0
  END DO
END DO

! output of wind-stress forcing
OPEN(10,file ='taux.dat',form='formatted')
DO j = 1,ny
  WRITE(10,'(101F12.6)')(taux(j,k),k=1,nx)
END DO
CLOSE(10)

! initialise Lagrangian floats
ist = -1
randm = ran3(ist)

xlen = REAL(nx)*dx
ylen = REAL(ny)*dy
zlen = REAL(nz)*dz

count = 0

DO ii = 1,100*ntra
  xpos = ran3(ist)*xlen
  ypos = ran3(ist)*ylen
  zpos = (0.1 + 0.8*ran3(ist))*zlen
  ipos = INT(zpos/dz+1)
  jpos = INT(ypos/dy+0.5)+1
  kpos = INT(xpos/dx+0.5)+1
  IF(wet(ipos,jpos,kpos))THEN
   count = count + 1
   IF(count <= ntra)THEN
       tra(count,1) = zpos
       tra(count,2) = ypos
       tra(count,3) = xpos
   ELSE
    WRITE(6,*) "Float initialisation finished"
    GOTO 99
   END IF
  END IF
END DO

WRITE(6,*) "WARNING: Some floats were left over"

99 CONTINUE

PAUSE

OPEN(70,file ='TRx.dat',form='formatted',recl = 10000000)
OPEN(80,file ='TRy.dat',form='formatted',recl = 10000000)
OPEN(90,file ='TRz.dat',form='formatted',recl = 10000000)

! runtime parameters
ntot = INT(60.*24.*60.*60./dt)
time = 0.0

! output frequency
nout = INT(12.*60.*60./dt)

! open files for output
OPEN(10,file ='rhoS.dat',form='formatted')
OPEN(20,file ='uS.dat',form='formatted')
OPEN(30,file ='vS.dat',form='formatted')
OPEN(40,file ='eta.dat',form='formatted')
OPEN(50,file ='cS.dat',form='formatted')

OPEN(11,file ='rho2.dat',form='formatted')
OPEN(21,file ='u2.dat',form='formatted')
OPEN(31,file ='v2.dat',form='formatted')
OPEN(51,file ='c2.dat',form='formatted')

OPEN(14,file ='rho3.dat',form='formatted')
OPEN(24,file ='u3.dat',form='formatted')
OPEN(34,file ='v3.dat',form='formatted')
OPEN(54,file ='c3.dat',form='formatted')

OPEN(15,file ='rho4.dat',form='formatted')
OPEN(25,file ='u4.dat',form='formatted')
OPEN(35,file ='v4.dat',form='formatted')
OPEN(55,file ='c4.dat',form='formatted')

OPEN(12,file ='rhoB.dat',form='formatted')
OPEN(22,file ='uB.dat',form='formatted')
OPEN(32,file ='vB.dat',form='formatted')
OPEN(52,file ='cB.dat',form='formatted')

OPEN(13,file ='rhoT.dat',form='formatted')
OPEN(23,file ='uT.dat',form='formatted')
OPEN(33,file ='vT.dat',form='formatted')
OPEN(53,file ='cT.dat',form='formatted')

! output of initial distributions
DO j = 1,ny
  WRITE(10,'(101F12.6)')(rho(1,j,k),k=1,nx)
  WRITE(20,'(101F12.6)')(u(1,j,k),k=1,nx)
  WRITE(30,'(101F12.6)')(v(1,j,k),k=1,nx)
  WRITE(50,'(101F12.6)')(c(1,j,k),k=1,nx)
  WRITE(40,'(101F12.6)')(q(0,j,k)/(RHOREF*G),k=1,nx)
END DO

DO j = 1,ny
  WRITE(11,'(101F12.6)')(rho(2,j,k),k=1,nx)
  WRITE(21,'(101F12.6)')(u(2,j,k),k=1,nx)
  WRITE(31,'(101F12.6)')(v(2,j,k),k=1,nx)
  WRITE(51,'(101F12.6)')(c(2,j,k),k=1,nx)
END DO

DO j = 1,ny
  WRITE(14,'(101F12.6)')(rho(3,j,k),k=1,nx)
  WRITE(24,'(101F12.6)')(u(3,j,k),k=1,nx)
  WRITE(34,'(101F12.6)')(v(3,j,k),k=1,nx)
  WRITE(54,'(101F12.6)')(c(3,j,k),k=1,nx)
END DO

DO j = 1,ny
  WRITE(15,'(101F12.6)')(rho(4,j,k),k=1,nx)
  WRITE(25,'(101F12.6)')(u(4,j,k),k=1,nx)
  WRITE(35,'(101F12.6)')(v(4,j,k),k=1,nx)
  WRITE(55,'(101F12.6)')(c(4,j,k),k=1,nx)
END DO

DO j = 1,ny
  WRITE(12,'(101F12.6)')(rho(ib(j,k),j,k),k=1,nx)
  WRITE(22,'(101F12.6)')(u(ibu(j,k),j,k),k=1,nx)
  WRITE(32,'(101F12.6)')(v(ibv(j,k),j,k),k=1,nx)
  WRITE(52,'(101F12.6)')(c(ib(j,k),j,k),k=1,nx)
END DO

DO i = 1,nz
  WRITE(13,'(101F12.6)')(rho(i,j,50),j=1,ny)
  WRITE(23,'(101F12.6)')(u(i,j,50),j=1,ny)
  WRITE(33,'(101F12.6)')(v(i,j,50),j=1,ny)
  WRITE(53,'(101F12.6)')(c(i,j,50),j=1,ny)
END DO

!---------------------------
! simulation loop
!---------------------------
DO n = 1,ntot

time = time + dt
write(6,*)"time (days)", time/(24.*3600.)
ad = MIN(time/(5.*24.0*3600.0),1.0) ! adjustment parameter

IF(time > 5.*24.0*3600.0) ad = 0.0

! prognostic equations
CALL dyn

DO ii = 1,ntra

! locate grid cell of tracer

ipos = INT(tra(ii,1)/dz+1)
jpos = INT(tra(ii,2)/dy+0.5)+1
kpos = INT(tra(ii,3)/dx+0.5)+1

uu = 0.5*(u(ipos,jpos,kpos)+u(ipos,jpos,kpos-1))
vv = 0.5*(v(ipos,jpos,kpos)+v(ipos,jpos-1,kpos))
ww = 0.5*(w(ipos,jpos,kpos)+w(ipos+1,jpos,kpos))

! avoid stranding
IF( dry(ipos,jpos,kpos).AND.uu>0 ) uu = 0.0
IF( dry(ipos,jpos,kpos-1).AND.uu<0 ) uu = 0.0
IF( dry(ipos,jpos,kpos).AND.vv>0 ) vv = 0.0
IF( dry(ipos,jpos-1,kpos).AND.vv<0 ) vv = 0.0
IF( dry(ipos+1,jpos,kpos).AND.ww<0 ) ww = 0.0
IF( ipos==1.AND.ww>0 ) ww = 0.0
IF( ipos==nz.AND.ww<0 ) ww = 0.0

! change of location
tra(ii,1) = tra(ii,1)-dt*ww
IF(tra(ii,1) < 0.0) tra(ii,1) = - tra(ii,1)
IF(tra(ii,1) > zlen) tra(ii,1) = 2.0*zlen - tra(ii,1)
tra(ii,2) = tra(ii,2)+dt*vv
IF(tra(ii,2) < 0.0) tra(ii,2) = ylen + tra(ii,2)
IF(tra(ii,2) > ylen) tra(ii,2) = tra(ii,2) - ylen
tra(ii,3) = tra(ii,3)+dt*uu
IF(tra(ii,3) < 0.0) tra(ii,3) = - tra(ii,3)
IF(tra(ii,3) > xlen) tra(ii,3) = 2.0*xlen - tra(ii,3)

END DO

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
  WRITE(11,'(101F12.6)')(rho(2,j,k),k=1,nx)
  WRITE(21,'(101F12.6)')(u(2,j,k),k=1,nx)
  WRITE(31,'(101F12.6)')(v(2,j,k),k=1,nx)
  WRITE(51,'(101F12.6)')(c(2,j,k),k=1,nx)
END DO

DO j = 1,ny
  WRITE(12,'(101F12.6)')(rho(ib(j,k),j,k),k=1,nx)
  WRITE(22,'(101F12.6)')(u(ibu(j,k),j,k),k=1,nx)
  WRITE(32,'(101F12.6)')(v(ibv(j,k),j,k),k=1,nx)
  WRITE(52,'(101F12.6)')(c(ib(j,k),j,k),k=1,nx)
END DO

DO i = 1,nz
  WRITE(13,'(101F12.6)')(rho(i,j,50),j=1,ny)
  WRITE(23,'(101F12.6)')(u(i,j,50),j=1,ny)
  WRITE(33,'(101F12.6)')(v(i,j,50),j=1,ny)
  WRITE(53,'(101F12.6)')(c(i,j,50),j=1,ny)
END DO

DO j = 1,ny
  WRITE(14,'(101F12.6)')(rho(3,j,k),k=1,nx)
  WRITE(24,'(101F12.6)')(u(3,j,k),k=1,nx)
  WRITE(34,'(101F12.6)')(v(3,j,k),k=1,nx)
  WRITE(54,'(101F12.6)')(c(3,j,k),k=1,nx)
END DO

DO j = 1,ny
  WRITE(15,'(101F12.6)')(rho(4,j,k),k=1,nx)
  WRITE(25,'(101F12.6)')(u(4,j,k),k=1,nx)
  WRITE(35,'(101F12.6)')(v(4,j,k),k=1,nx)
  WRITE(55,'(101F12.6)')(c(4,j,k),k=1,nx)
END DO

  WRITE(70,'(5000F12.6)')(tra(ii,3)/1000.0,ii=1,ntra)
  WRITE(80,'(5000F12.6)')(tra(ii,2)/1000.0,ii=1,ntra)
  WRITE(90,'(5000F12.6)')(tra(ii,1),ii=1,ntra)


  WRITE(6,*)"Data output at time = ",time/(24.*3600.)
ENDIF

END DO ! end of iteration loop

END PROGRAM full3D
