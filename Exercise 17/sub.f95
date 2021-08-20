MODULE sub
USE param
USE random

CONTAINS

!******************
SUBROUTINE init
!local parameters
REAL :: ccc, cmc, depth(0:nx+1),uu
INTEGER :: nb
REAL :: N2

azmin = 0.0005 ! background viscosity/diffusivity

! set initial arrays
DO i = 0,nz+1
DO k = 0,nx+1
 p(i,k) = 0.0
 q(i,k) = 0.0
 dq(i,k) = 0.0
 rho(i,k) = 0.0
 u(i,k) = 0.0
 un(i,k) = 0.0
 c1(i,k) = 0.0
 c1n(i,k) = 0.0
 c2(i,k) = 0.0
 c2n(i,k) = 0.0
 c3(i,k) = 0.0
 c3n(i,k) = 0.0
 ustar(i,k) = 0.0
 v(i,k) = 0.0
 vn(i,k) = 0.0
 vstar(i,k) = 0.0
 w(i,k) = 0.0
 wn(i,k) = 0.0
 wstar(i,k) = 0.0
 wet(i,k) = .true.
 az(i,k) = azmin
 kz(i,k) = azmin
END DO
END DO

DO k = 0,nx+1
 usum(k) = 0.0
END DO

! grid parameters
dx = 100.0
dz = 2.0
dt = 1.5
drdx = 1.0/(RHOREF*dx)
drdz = 1.0/(RHOREF*dz)
f = 1.e-4
alpha = dt*f

kh = 1.e-4
ah = 1.e-4
r = 1.e-3

N2 = 1.e-6 ! stability frequency squared

!initialisation of random function
ist = -1
randm = ran3(ist)

depth(0) = 0.0
depth(nx+1) = 0.0

DO k = 1,nx
 depth(k) = 100.0 - 50.0*REAL(k-1)/REAL(nx-1)
END DO

OPEN(50,file ='h.dat',form='formatted')
  WRITE(50,'(201F12.6)')(depth(k),k=1,nx)
CLOSE(50)

DO k = 0,nx+1
  nb = INT(depth(k)/dz)
  nb = MIN(nb,nz)
  DO i = nb+1,nz+1
    wet(i,k) = .false.
  END DO
! wet(0,k) = .false.
END DO

DO k = 0,nx+1
DO i = 0,nz+1
  dry(i,k)=.NOT.wet(i,k)
END DO
END DO

DO k = 0,nx+1
DO i = 1,10
  IF(wet(i,k))rho(i,k) = 0.0
END DO
DO i = 11,nz+1
  IF(wet(i,k))rho(i,k) = 1.0
END DO
END DO

! coefficients for SOR
omega = 1.4
peps = 1.e-3

DO k = 1,nx
DO i = 1,nz
  at(i,k) = dx/dz
  ab(i,k) = dx/dz
  ae(i,k) = dz/dx
  aw(i,k) = dz/dx
  IF(dry(i,k+1)) ae(i,k) = 0.0
  IF(dry(i,k-1)) aw(i,k) = 0.0
  IF(dry(i+1,k)) ab(i,k) = 0.0
  IF(dry(i-1,k)) at(i,k) = 0.0
  atot(i,k) = ab(i,k)+at(i,k)+ae(i,k)+aw(i,k)
END DO
END DO

! bottom pointer
DO i = 0,nz
DO k = 0,nx+1
 ib(k) = 0
END DO
END DO

DO i = 0,nz
DO k = 0,nx+1
  IF(wet(i,k).AND.dry(i+1,k))ib(k) = i
END DO
END DO

! initial distributions of concentration fields
DO k = 0,60
  DO i = 0,10
    IF(wet(i,k))c1(i,k) = 1.0
  END DO
  IF(ib(k)>9)THEN
    DO i = ib(k)-9,ib(k)
      IF(wet(i,k))c2(i,k) = 1.0
    END DO
  END IF
END DO

DO k = 71,nx+1
DO i = 0,ib(k)
  IF(wet(i,k))c3(i,k) = 1.0
END DO
END DO 

END SUBROUTINE init

!*******************
SUBROUTINE dyn

! local parameters
REAL :: pressx, pressz
REAL :: perr, q1, q2, term1
INTEGER :: nsor, nstop
REAL :: advx(0:nz+1,0:nx+1), advy(0:nz+1,0:nx+1), advz(0:nz+1,0:nx+1)
REAL :: div, div1, div2, um, vm, unn
REAL :: dif1, dif2, difh, difv, diff, diffu, diffv, diffw
REAL :: drho, speed, aztop, azbot

! density prediction
DO k = 0,nx+1
DO i = 0,nz+1
  CuP(i,k) = 0.5*(u(i,k)+abs(u(i,k)))*dt/dx
  CuN(i,k) = 0.5*(u(i,k)-abs(u(i,k)))*dt/dx
  CwP(i,k) = 0.5*(w(i,k)+abs(w(i,k)))*dt/dz
  CwN(i,k) = 0.5*(w(i,k)-abs(w(i,k)))*dt/dz
END DO
END DO

DO i = 0,nz+1
DO k = 0,nx+1
  B(i,k) = rho(i,k)
END DO
END DO

CALL advect

DO k = 1,nx
DO i = 1,nz
 div1 = (u(i,k)-u(i,k-1))/dx
 div2 = (w(i,k)-w(i+1,k))/dz
 div = dt*B(i,k)*(div1+div2)
 dif1 = 0.0
 IF(WET(i,k+1)) dif1 = (rho(i,k+1)-rho(i,k))/dx
 dif2 = 0.0
 IF(WET(i,k-1)) dif2 = (rho(i,k)-rho(i,k-1))/dx
 difh = kh*(dif1-dif2)/dx
 aztop = 0.5*(kz(i,k)+kz(i-1,k))
 dif1 = aztop*(rho(i-1,k)-rho(i,k))/dz
 IF(dry(i-1,k)) dif1 = 0.0
 azbot = 0.5*(kz(i,k)+kz(i+1,k))
 dif2 = azbot*(rho(i,k)-rho(i+1,k))/dz
 IF(dry(i+1,k)) dif2 = 0.0
 difz = (dif1-dif2)/dz
 diff = dt*(difh+difz)  
 rhon(i,k)= rho(i,k)+BN(i,k)+div+diff
END DO
END DO

! boundary conditions
DO i = 0,nz+1
 rhon(i,0) = rhon(i,1)
 rhon(i,nx+1) = rhon(i,nx)
END DO

DO k = 0,nx+1
 rhon(0,k) = rhon(1,k)
 rhon(nz+1,k) = rhon(nz,k)
END DO

! updating
DO k = 0,nx+1
DO i = 0,nz+1
  rho(i,k) = rhon(i,k)
END DO
END DO

! calculate hydrostatic pressure

DO k = 0,nx+1
  P(0,k) = 0.0
  DO i = 1,nz+1
    P(i,k) = P(i-1,k) + 0.5*(rho(i-1,k)+rho(i,k))*G*dz
  END DO
END DO

! tracer prediction
DO i = 0,nz+1
DO k = 0,nx+1
  B(i,k) = c1(i,k)
END DO
END DO

CALL advect

DO i = 1,nz
DO k = 1,nx
 div1 = (u(i,k)-u(i,k-1))/dx
 div2 = (w(i,k)-w(i+1,k))/dz
 div = dt*B(i,k)*(div1+div2)
 dif1 = 0.0
 IF(WET(i,k+1)) dif1 = (c1(i,k+1)-c1(i,k))/dx
 dif2 = 0.0
 IF(WET(i,k-1)) dif2 = (c1(i,k)-c1(i,k-1))/dx
 difh = kh*(dif1-dif2)/dx
 aztop = 0.5*(kz(i,k)+kz(i-1,k))
 dif1 = aztop*(c1(i-1,k)-c1(i,k))/dz
 IF(dry(i-1,k)) dif1 = 0.0 
 azbot = 0.5*(kz(i,k)+kz(i+1,k))
 dif2 = azbot*(c1(i,k)-c1(i+1,k))/dz
 IF(dry(i+1,k)) dif2 = 0.0 
 difz = (dif1-dif2)/dz
 diff = dt*(difh+difz)  
 c1n(i,k)= c1(i,k)+BN(i,k)+div+diff
END DO
END DO

! boundary conditions
DO i = 0,nz+1
 c1n(i,0) = c1n(i,1)
 c1n(i,nx+1) = c1n(i,nx)
END DO

DO k = 0,nx+1
 c1n(0,k) = c1n(1,k)
 c1n(nz+1,k) = c1n(nz,k)
END DO

! updating
DO k = 0,nx+1
DO i = 0,nz+1
 c1(i,k) = c1n(i,k)
END DO
END DO

DO i = 0,nz+1
DO k = 0,nx+1
  B(i,k) = c2(i,k)
END DO
END DO

CALL advect

DO i = 1,nz
DO k = 1,nx
 div1 = (u(i,k)-u(i,k-1))/dx
 div2 = (w(i,k)-w(i+1,k))/dz
 div = dt*B(i,k)*(div1+div2)
 dif1 = 0.0
 IF(WET(i,k+1)) dif1 = (c2(i,k+1)-c2(i,k))/dx
 dif2 = 0.0
 IF(WET(i,k-1)) dif2 = (c2(i,k)-c2(i,k-1))/dx
 difh = kh*(dif1-dif2)/dx
 aztop = 0.5*(kz(i,k)+kz(i-1,k))
 dif1 = aztop*(c2(i-1,k)-c2(i,k))/dz
 IF(dry(i-1,k)) dif1 = 0.0 
 azbot = 0.5*(kz(i,k)+kz(i+1,k))
 dif2 = azbot*(c2(i,k)-c2(i+1,k))/dz
 IF(dry(i+1,k)) dif2 = 0.0 
 difz = (dif1-dif2)/dz
 diff = dt*(difh+difz)  
 c2n(i,k)= c2(i,k)+BN(i,k)+div+diff
END DO
END DO

! boundary conditions
DO i = 0,nz+1
 c2n(i,0) = c2n(i,1)
 c2n(i,nx+1) = c2n(i,nx)
END DO

DO k = 0,nx+1
 c2n(0,k) = c2n(1,k)
 c2n(nz+1,k) = c2n(nz,k)
END DO

! updating
DO k = 0,nx+1
DO i = 0,nz+1
 c2(i,k) = c2n(i,k)
END DO
END DO

DO i = 0,nz+1
DO k = 0,nx+1
  B(i,k) = c3(i,k)
END DO
END DO

CALL advect

DO i = 1,nz
DO k = 1,nx
 div1 = (u(i,k)-u(i,k-1))/dx
 div2 = (w(i,k)-w(i+1,k))/dz
 div = dt*B(i,k)*(div1+div2)
 dif1 = 0.0
 IF(WET(i,k+1)) dif1 = (c3(i,k+1)-c3(i,k))/dx
 dif2 = 0.0
 IF(WET(i,k-1)) dif2 = (c3(i,k)-c3(i,k-1))/dx
 difh = kh*(dif1-dif2)/dx
 aztop = 0.5*(kz(i,k)+kz(i-1,k))
 dif1 = aztop*(c3(i-1,k)-c3(i,k))/dz
 IF(dry(i-1,k)) dif1 = 0.0 
 azbot = 0.5*(kz(i,k)+kz(i+1,k))
 dif2 = azbot*(c3(i,k)-c3(i+1,k))/dz
 IF(dry(i+1,k)) dif2 = 0.0 
 difz = (dif1-dif2)/dz
 diff = dt*(difh+difz)  
 c3n(i,k)= c3(i,k)+BN(i,k)+div+diff
END DO
END DO

! boundary conditions
DO i = 0,nz+1
 c3n(i,0) = c3n(i,1)
 c3n(i,nx+1) = c3n(i,nx)
END DO

DO k = 0,nx+1
 c3n(0,k) = c3n(1,k)
 c3n(nz+1,k) = c3n(nz,k)
END DO

! updating
DO k = 0,nx+1
DO i = 0,nz+1
 c3(i,k) = c3n(i,k)
END DO
END DO

! calculate the nonlinear terms for v-momentum equation
DO i = 0,nz+1
DO k = 0,nx+1
  B(i,k) = v(i,k)
END DO
END DO

CALL advect

DO i = 1,nz
DO k = 1,nx
  div1 = (u(i,k)-u(i,k-1))/dx
  div2 = (w(i,k)-w(i+1,k))/dz
  div = dt*B(i,k)*(div1+div2)  
  advy(i,k)= BN(i,k)+div
END DO
END DO

! calculate the nonlinear terms for u-momentum equation
DO i = 0,nz+1
DO k = 0,nx
  CuP(i,k) = 0.25*(u(i,k)+u(i,k+1)+abs(u(i,k))+abs(u(i,k+1)))*dt/dx
  CuN(i,k) = 0.25*(u(i,k)+u(i,k+1)-abs(u(i,k))-abs(u(i,k+1)))*dt/dx
  CwP(i,k) = 0.25*(w(i,k)+w(i,k+1)+abs(w(i,k))+abs(w(i,k+1)))*dt/dz
  CwN(i,k) = 0.25*(w(i,k)+w(i,k+1)-abs(w(i,k))-abs(w(i,k+1)))*dt/dz
END DO
END DO

DO i = 0,nz+1
DO k = 0,nx+1
  B(i,k) = u(i,k)
END DO
END DO

CALL advect

DO i = 1,nz
DO k = 1,nx
  div1 = 0.5*(u(i,k+1)-u(i,k-1))/dx
  div2 = 0.5*(w(i,k)+w(i,k+1)-w(i+1,k)-w(i+1,k+1))/dz
  div = dt*B(i,k)*(div1+div2)  
  advx(i,k)= BN(i,k)+div
END DO
END DO

! calculate the nonlinear terms for w-momentum equation
DO i = 0,nz
DO k = 0,nx+1
  CuP(i,k) = 0.25*(u(i,k)+u(i-1,k)+abs(u(i,k))+abs(u(i-1,k)))*dt/dx
  CuN(i,k) = 0.25*(u(i,k)+u(i-1,k)-abs(u(i,k))-abs(u(i-1,k)))*dt/dx
  CwP(i,k) = 0.25*(w(i,k)+w(i-1,k)+abs(w(i,k))+abs(w(i-1,k)))*dt/dz
  CwN(i,k) = 0.25*(w(i,k)+w(i-1,k)-abs(w(i,k))-abs(w(i-1,k)))*dt/dz
END DO
END DO

DO i = 0,nz+1
DO k = 0,nx+1
  B(i,k) = w(i,k)
END DO
END DO

CALL advect

DO i = 1,nz
DO k = 1,nx
  div1 = 0.5*(u(i,k)+u(i-1,k)-u(i,k-1)-u(i-1,k-1))/dx
  div2 = 0.5*(w(i-1,k)-w(i+1,k))/dz
  div = dt*B(i,k)*(div1+div2)  
  advz(i,k)= BN(i,k) + div
END DO
END DO

! calculate bottom shear-stress components
DO k = 1,nx
  IF(ib(k)>0)THEN
    um = u(ib(k),k)
    vm = 0.5*(v(ib(k),k)+v(ib(k),k+1))
    speed = SQRT(um*um+vm*vm)
    tbu(k) = r*um*speed
    vm = v(ib(k),k)
    um = 0.5*(u(ib(k),k)+u(ib(k),k-1))
    speed = SQRT(um*um+vm*vm)
    tbv(k) = r*vm*speed
  END IF
END DO

! calculate ustar, vstar and wstar 
DO i = 1,nz
DO k = 1,nx

! diffusion of u
 dif1 = 0.0
 IF(WET(i,k+1)) dif1 = (u(i,k+1)-u(i,k))/dx
 dif2 = 0.0
 IF(WET(i,k-1)) dif2 = (u(i,k)-u(i,k-1))/dx
 difh = ah*(dif1-dif2)/dx
 aztop = 0.25*(az(i,k)+az(i-1,k)+az(i,k+1)+az(i-1,k+1))
 dif1 = aztop*(u(i-1,k)-u(i,k))/dz
 IF(dry(i-1,k)) dif1 = 0.0 
 azbot = 0.25*(az(i,k)+az(i+1,k)+az(i,k+1)+az(i+1,k+1))
 dif2 = azbot*(u(i,k)-u(i+1,k))/dz
 IF(i == ib(k) ) dif2 = tbu(k)
 difz = (dif1-dif2)/dz
 diffu = dt*(difh+difz)  

! diffusion of v
 dif1 = 0.0
 IF(WET(i,k+1)) dif1 = (v(i,k+1)-v(i,k))/dx
 dif2 = 0.0
 IF(WET(i,k-1)) dif2 = (v(i,k)-v(i,k-1))/dx
 difh = ah*(dif1-dif2)/dx
 aztop = 0.5*(az(i,k)+az(i-1,k))
 dif1 = aztop*(v(i-1,k)-v(i,k))/dz
 IF(dry(i-1,k)) dif1 = 0.0 
 azbot = 0.5*(az(i,k)+az(i+1,k))
 dif2 = azbot*(v(i,k)-v(i+1,k))/dz
 IF(i == ib(k) ) dif2 = tbv(k)

 difz = (dif1-dif2)/dz
 diffv = dt*(difh+difz)  

! diffusion of w
 dif1 = 0.0
 IF(WET(i,k+1)) dif1 = (w(i,k+1)-w(i,k))/dx
 dif2 = 0.0
 IF(WET(i,k-1)) dif2 = (w(i,k)-w(i,k-1))/dx
 difh = ah*(dif1-dif2)/dx
 aztop = az(i-1,k)
 dif1 = aztot*(w(i-1,k)-w(i,k))/dz
 IF(dry(i-1,k)) dif1 = 0.0
 azbot = az(i,k)
 dif2 = azbot*(w(i,k)-w(i+1,k))/dz
 IF(dry(i+1,k)) dif2 = 0.0 
 difz = (dif1-dif2)/dz
 diffw = dt*(difh+difz)  

  IF(wet(i,k))THEN
    um = 0.5*(u(i,k)+u(i,k-1))
    vm = 0.5*(v(i,k)+v(i,k+1))
    pressx = -drdx*(q(i,k+1)-q(i,k))-drdx*(p(i,k+1)-p(i,k))
    IF(wet(i,k+1)) ustar(i,k) = cos(alpha)*u(i,k) + sin(alpha)*vm + dt*pressx + advx(i,k) + diffu
    vstar(i,k) = cos(alpha)*v(i,k) - sin(alpha)*um + advy(i,k) + diffv
    pressz = -drdz*(q(i-1,k)-q(i,k))
    IF(wet(i-1,k)) wstar(i,k) = w(i,k) + dt*pressz + advz(i,k) + diffw
  END IF

END DO
END DO

! boundary conditions
DO i = 1,nz
 ustar(i,0) = 0.0
 ustar(i,nx) = 0.0
 ustar(i,nx+1) = 0.0 
 vstar(i,0) = 0.0
 vstar(i,nx+1) = 0.0
 wstar(i,0) = 0.0
 wstar(i,nx+1) = 0.0
END DO

! calculate right-hand side of Poisson equation
DO i = 1,nz
DO k = 1,nx
 qstar(i,k) = -1.*RHOREF/dt*(  &
 &  (ustar(i,k)-ustar(i,k-1))*dz + &
 &  (wstar(i,k)-wstar(i+1,k))*dx )
END DO
END DO

! STEP 4: S.O.R. ITERATION

nstop = 8000

!*****************
DO nsor = 1,8000
!*****************

perr = 0.0

! STEP 1: predict pressure correction
DO i = 1,nz
DO k = 1,nx

 IF(wet(i,k))THEN
 
 q1 = dq(i,k)
 term1 = qstar(i,k) + & 
  &      at(i,k)*dq(i-1,k) + ab(i,k)*dq(i+1,k) + & 
  &      aw(i,k)*dq(i,k-1) + ae(i,k)*dq(i,k+1)
 q2 = (1.0-omega)*q1 + omega*term1/atot(i,k) 
 dq(i,k) = q2
 IF(k == 1) dq(i,nx+1) = q2 ! boundary condition
 IF(k == nx) dq(i,0) = q2 ! boundary condition
 perr = MAX(ABS(q2-q1),perr)

 END IF

END DO
END DO

! STEP 2: predict new velocities 
DO i = 1,nz
DO k = 1,nx
  IF(wet(i,k))THEN
    pressx = -drdx*(dq(i,k+1)-dq(i,k))
    IF(wet(i,k+1))un(i,k) = ustar(i,k) + dt*pressx
    vn(i,k) = vstar(i,k)
    pressz = -drdz*(dq(i-1,k)-dq(i,k))
    IF(wet(i-1,k)) wn(i,k) = wstar(i,k) + dt*pressz
  END IF
END DO
END DO

! STEP 3a: predict depth-integrated flow
DO k = 1,nx
  usum(k) = 0.
  DO i = 1,nz
    usum(k) = usum(k) + dz*un(i,k)
  END DO
END DO

! lateral boundary conditions
 usum(0) = 0.0
 usum(nx) = 0.0
 usum(nx+1) = 0.0

! STEP 3b: predict surface pressure field
DO k = 1,nx
 dq(0,k) = -dt*RHOREF*G*(usum(k)-usum(k-1))/dx
END DO

IF(perr <= peps)THEN
  nstop = nsor
  GOTO 33
END IF

!********************
END DO
!********************

    GOTO 34
 33 WRITE(*,*) "No. of Interactions =>", nstop
 34 CONTINUE

! updating for next time step
DO i = 1,nz
DO k = 1,nx
  q(i,k) = q(i,k)+dq(i,k)
  u(i,k) = un(i,k)
  v(i,k) = vn(i,k)
  w(i,k) = wn(i,k)
END DO
END DO

DO k = 1,nx
  q(0,k) = q(0,k)+dq(0,k)
END DO

! lateral boundary conditions
 q(0,0) = q(0,1)
 q(0,nx+1) = q(0,nx)

DO i = 1,nz
 u(i,0) = 0.0
 u(i,nx) = 0.0
 u(i,nx+1) = 0.0
 v(i,0) = 0.0
 v(i,nx) = 0.0
 w(i,0) = 0.0
 w(i,nx+1) = 0.0
 q(i,0) = q(i,1)
 q(i,nx+1) = q(i,nx)
END DO

RETURN
END SUBROUTINE dyn

SUBROUTINE advect
! local parameters
REAL :: RxP(0:nz+1,0:nx+1), RxN(0:nz+1,0:nx+1)
REAL :: RzP(0:nz+1,0:nx+1), RzN(0:nz+1,0:nx+1)
REAL :: dB, term1, term2, term3, term4
REAL :: BwP, BwN, BeP, BeN, BbP, BbN, BtP, BtN 

DO i = 0,nz+1
DO k = 0,nx+1
  RxP(i,k) = 0.0
  RxN(i,k) = 0.0
  RzP(i,k) = 0.0
  RzN(i,k) = 0.0
END DO
END DO

DO i = 1,nz
DO k = 1,nx
  dB =  B(i,k+1)-B(i,k)
  IF(ABS(dB) > 0.0) RxP(i,k) = (B(i,k)-B(i,k-1))/dB
  dB =  B(i-1,k)-B(i,k)
  IF(ABS(dB) > 0.0) RzP(i,k) = (B(i,k)-B(i+1,k))/dB
END DO
END DO

DO i = 1,nz
DO k = 0,nx-1
  dB =  B(i,k+1)-B(i,k)
  IF(ABS(dB) > 0.0) RxN(i,k) = (B(i,k+2)-B(i,k+1))/dB
END DO
END DO

DO i = 2,nz+1
DO k = 1,nx
  dB =  B(i-1,k)-B(i,k)
  IF(ABS(dB) > 0.0) RzN(i,k) = (B(i-2,k)-B(i-1,k))/dB
END DO
END DO   

DO i = 1,nz
DO k = 1,nx

term1 = (1.0-CuP(i,k-1))*(B(i,k)-B(i,k-1))
BwP = B(i,k-1)+0.5*PSI(RxP(i,k-1))*term1

term1 = (1.0+CuN(i,k-1))*(B(i,k)-B(i,k-1))
BwN = B(i,k)-0.5*PSI(RxN(i,k-1))*term1

term1 = (1.0-CuP(i,k))*(B(i,k+1)-B(i,k))
BeP = B(i,k)+0.5*PSI(RxP(i,k))*term1

term1 = (1.0+CuN(i,k))*(B(i,k+1)-B(i,k))  
BeN = B(i,k+1)-0.5*PSI(RxN(i,k))*term1

term1 = (1.0-CwP(i+1,k))*(B(i,k)-B(i+1,k))
BbP = B(i+1,k)+0.5*PSI(RzP(i+1,k))*term1

term1 = (1.0+CwN(i+1,k))*(B(i,k)-B(i+1,k))  
BbN = B(i,k)-0.5*PSI(RzN(i+1,k))*term1

term1 = (1.0-CwP(i,k))*(B(i-1,k)-B(i,k)) 
BtP = B(i,k)+0.5*PSI(RzP(i,k))*term1

term1 = (1.0+CwN(i,k))*(B(i-1,k)-B(i,k)) 
BtN = B(i-1,k)-0.5*PSI(RzN(i,k))*term1

term1 = CuP(i,k-1)*BwP+CuN(i,k-1)*BwN
term2 = CuP(i,k)*BeP+CuN(i,k)*BeN
term3 = CwP(i+1,k)*BbP+CwN(i+1,k)*BbN
term4 = CwP(i,k)*BtP+CwN(i,k)*BtN

BN(i,k) = term1-term2+term3-term4

END DO
END DO

RETURN

END SUBROUTINE advect

REAL FUNCTION psi(r)

! input parameters
REAL, INTENT(IN) :: r  

! local parameters 
REAL :: term1, term2, term3

  term1 = MIN(2.0*r,1.0)
  term2 = MIN(r,2.0)
  term3 = MAX(term1,term2)
  psi = MAX(term3,0.0)

RETURN

END FUNCTION psi 

END MODULE sub