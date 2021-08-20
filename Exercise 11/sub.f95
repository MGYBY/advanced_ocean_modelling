MODULE sub
USE param
USE random

CONTAINS

!******************
SUBROUTINE init
!local parameters
REAL :: ccc, cmc, depth(nx), delsdelt
INTEGER :: nb
REAL :: N2

! set initial arrays
DO i = 0,nz+1
DO k = 0,nx+1
 p(i,k) = 0.0
 t(i,k) = 0.0
 tn(i,k) = 0.0
 s(i,k) = 0.0
 sn(i,k) = 0.0
 q(i,k) = 0.0
 dq(i,k) = 0.0
 rho(i,k) = 0.0
 u(i,k) = 0.0
 un(i,k) = 0.0
 ustar(i,k) = 0.0
 w(i,k) = 0.0
 wn(i,k) = 0.0
 wstar(i,k) = 0.0
 wet(i,k) = .true.
 c(i,k) = 0.0
 cn(i,k) = 0.0
END DO
END DO

delsdelt = 2.5/8.0
ist = -1
randm = ran3(ist)

! initial temperature and salinity fields
DO i = 1,10
DO k = 0,nx+1
 c(i,k) = 1.0
 cn(i,k) = 1.0
 t(i,k) = 10.0+ran3(ist)*0.0001
 tn(i,k) = t(i,k)
 s(i,k) = 10.0*delsdelt
 sn(i,k) = s(i,k)
 rho(i,k) = density(t(i,k),s(i,k))
END DO
END DO

DO k = 0,nx+1
 usum(k) = 0.0
END DO

! grid parameters
dx = 1.0
dz = 1.0
dt = 1.0
drdx = 1.0/(RHOREF*dx)
drdz = 1.0/(RHOREF*dz)
! diffusivity for heat (artificially increase by factor of 1000)
kht = 1000.*1.e-7
kzt = 1000.*1.e-7
! diffusivity for salt (artificially increase by factor of 1000)
khs = 1000.*1.e-9
kzs = 1000.*1.e-9
! ambient viscosity
ah = 1.e-4
az = 1.e-4
r = 1.e-3

! random initialisation of floats

!initialisation of random function
ist = -1
randm = ran3(ist)
xlen = REAL(nx)*dx
zlen = REAL(nz-1)*dz

DO ii = 1,ntra
  xpos = ran3(ist)*xlen
  zpos = 0.5*ran3(ist)*zlen
  tra(ii,1) = zpos
  tra(ii,2) = xpos
END DO

DO k = 1,nx
 depth(k) = 20.0
END DO

OPEN(50,file ='h.dat',form='formatted')
  WRITE(50,'(201F12.6)')(depth(k),k=1,nx)
CLOSE(50)

DO k = 1,nx
  nb = INT(depth(k)/dz)
  nb = MIN(nb,nz+1)
  DO i = nb,nz+1
    wet(i,k) = .false.
  END DO
  wet(0,k) = .false.
END DO

! coefficients for SOR
omega = 1.4
peps = 1.e-3

DO i = 1,nz
DO k = 1,nx
  at(i,k) = dx/dz
  ab(i,k) = dx/dz
  ae(i,k) = dz/dx
  aw(i,k) = dz/dx
  IF(.not.wet(i+1,k)) ab(i,k) = 0.0
  IF(.not.wet(i-1,k)) at(i,k) = 0.0
  atot(i,k) = ab(i,k)+at(i,k)+ae(i,k)+aw(i,k)
END DO
END DO

END SUBROUTINE init

!*******************
SUBROUTINE dyn

! local parameters
REAL :: pressx, pressz
REAL :: perr, q1, q2, term1
INTEGER :: nsor, nstop
REAL :: advx(0:nz+1,0:nx+1), advz(0:nz+1,0:nx+1)
REAL :: div, div1, div2
REAL :: dif1, dif2, difh, difv, diff, diffu, diffw
REAL :: drho, speed

! tracer prediction
DO i = 0,nz+1
DO k = 0,nx+1
  CuP(i,k) = 0.5*(u(i,k)+abs(u(i,k)))*dt/dx
  CuN(i,k) = 0.5*(u(i,k)-abs(u(i,k)))*dt/dx
  CwP(i,k) = 0.5*(w(i,k)+abs(w(i,k)))*dt/dz
  CwN(i,k) = 0.5*(w(i,k)-abs(w(i,k)))*dt/dz
END DO
END DO

DO i = 0,nz+1
DO k = 0,nx+1
  B(i,k) = c(i,k)
END DO
END DO

CALL advect

DO i = 1,nz
DO k = 1,nx
 div1 = (u(i,k)-u(i,k-1))/dx
 div2 = (w(i,k)-w(i+1,k))/dz
 div = dt*B(i,k)*(div1+div2)
 dif1 = 0.0
 IF(WET(i,k+1)) dif1 = (c(i,k+1)-c(i,k))/dx
 dif2 = 0.0
 IF(WET(i,k-1)) dif2 = (c(i,k)-c(i,k-1))/dx
 difh = kht*(dif1-dif2)/dx
 dif1 = 0.0
 IF(WET(i-1,k)) dif1 = (c(i-1,k)-c(i,k))/dz
 dif2 = 0.0
 IF(WET(i+1,k)) dif2 = (c(i,k)-c(i+1,k))/dz
 difz = kzt*(dif1-dif2)/dz
 diff = dt*(difh+difz)  
 cn(i,k)= c(i,k)+BN(i,k)+div+diff
END DO
END DO

! boundary conditions
DO i = 0,nz+1
 cn(i,0) = cn(i,nx)
 cn(i,nx+1) = cn(i,1)
END DO

DO k = 0,nx+1
 cn(0,k) = cn(1,k)
 cn(nz+1,k) = cn(nz,k)
END DO

! updating
DO k = 0,nx+1
DO i = 0,nz+1
 c(i,k) = cn(i,k)
END DO
END DO

! temperature prediction
DO i = 0,nz+1
DO k = 0,nx+1
  B(i,k) = t(i,k)
END DO
END DO

CALL advect

DO i = 1,nz
DO k = 1,nx
 div1 = (u(i,k)-u(i,k-1))/dx
 div2 = (w(i,k)-w(i+1,k))/dz
 div = dt*B(i,k)*(div1+div2)
 dif1 = 0.0
 IF(WET(i,k+1)) dif1 = (t(i,k+1)-t(i,k))/dx
 dif2 = 0.0
 IF(WET(i,k-1)) dif2 = (t(i,k)-t(i,k-1))/dx
 difh = kht*(dif1-dif2)/dx
 dif1 = 0.0
 IF(WET(i-1,k)) dif1 = (t(i-1,k)-t(i,k))/dz
 dif2 = 0.0
 IF(WET(i+1,k)) dif2 = (t(i,k)-t(i+1,k))/dz
 difz = kzt*(dif1-dif2)/dz
 diff = dt*(difh+difz)  
 tn(i,k)= t(i,k)+BN(i,k)+div+diff
END DO
END DO

! boundary conditions
DO i = 0,nz+1
 tn(i,0) = tn(i,nx)
 tn(i,nx+1) = tn(i,1)
END DO

DO k = 0,nx+1
 tn(0,k) = tn(1,k)
 tn(nz+1,k) = tn(nz,k)
END DO

! updating

DO k = 0,nx+1
DO i = 0,nz+1
 t(i,k) = tn(i,k)
END DO
END DO

! salinity prediction
DO i = 0,nz+1
DO k = 0,nx+1
  B(i,k) = s(i,k)
END DO
END DO

CALL advect

DO i = 1,nz
DO k = 1,nx
 div1 = (u(i,k)-u(i,k-1))/dx
 div2 = (w(i,k)-w(i+1,k))/dz
 div = dt*B(i,k)*(div1+div2)
 dif1 = 0.0
 IF(WET(i,k+1)) dif1 = (s(i,k+1)-s(i,k))/dx
 dif2 = 0.0
 IF(WET(i,k-1)) dif2 = (s(i,k)-s(i,k-1))/dx
 difh = khs*(dif1-dif2)/dx
 dif1 = 0.0
 IF(WET(i-1,k)) dif1 = (s(i-1,k)-s(i,k))/dz
 dif2 = 0.0
 IF(WET(i+1,k)) dif2 = (s(i,k)-s(i+1,k))/dz
 difz = kzs*(dif1-dif2)/dz
 diff = dt*(difh+difz)  
 sn(i,k)= s(i,k)+BN(i,k)+div+diff
END DO
END DO

! boundary conditions
DO i = 0,nz+1
 sn(i,0) = sn(i,nx)
 sn(i,nx+1) = sn(i,1)
END DO

DO k = 0,nx+1
 sn(0,k) = sn(1,k)
 sn(nz+1,k) = sn(nz,k)
END DO

! updating
DO k = 0,nx+1
DO i = 0,nz+1
 s(i,k) = sn(i,k)
END DO
END DO

! density calculation calling equation of state
DO k = 0,nx+1
DO i = 0,nz+1
 rho(i,k) = density(t(i,k),s(i,k))
END DO
END DO

! calculate hydrostatic pressure
DO k = 0,nx+1
  P(0,k) = 0.0
  DO i = 1,nz+1
    P(i,k) = P(i-1,k) + 0.5*(rho(i-1,k)+rho(i,k))*G*dz
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

! calculate ustar and vstar 
DO i = 1,nz
DO k = 1,nx

! diffusion of u
 dif1 = 0.0
 IF(WET(i,k+1)) dif1 = (u(i,k+1)-u(i,k))/dx
 dif2 = 0.0
 IF(WET(i,k-1)) dif2 = (u(i,k)-u(i,k-1))/dx
 difh = az*(dif1-dif2)/dx
 dif1 = 0.0
 IF(WET(i-1,k)) dif1 = ah*(u(i-1,k)-u(i,k))/dz
! bottom friction
 speed = ABS(u(i,k))
 dif2 = r*u(i,k)*speed
 IF(WET(i+1,k)) dif2 = az*(u(i,k)-u(i+1,k))/dz
 difz = (dif1-dif2)/dz
 diffu = dt*(difh+difz)  

! diffusion of w
 dif1 = 0.0
 IF(WET(i,k+1)) dif1 = (w(i,k+1)-w(i,k))/dx
 dif2 = 0.0
 IF(WET(i,k-1)) dif2 = (w(i,k)-w(i,k-1))/dx
 difh = az*(dif1-dif2)/dx
 dif1 = 0.0
 IF(WET(i-1,k)) dif1 = ah*(w(i-1,k)-w(i,k))/dz
 dif2 = 0.0
 IF(WET(i+1,k)) dif2 = az*(w(i,k)-w(i+1,k))/dz
 difz = (dif1-dif2)/dz
 diffw = dt*(difh+difz)  

  IF(wet(i,k))THEN
    pressx = -drdx*(q(i,k+1)-q(i,k))-drdx*(p(i,k+1)-p(i,k))
    IF(wet(i,k+1)) ustar(i,k) = u(i,k) + dt*pressx + advx(i,k) + diffu
    pressz = -drdz*(q(i-1,k)-q(i,k))
    IF(wet(i-1,k)) wstar(i,k) = w(i,k) + dt*pressz + advz(i,k) + diffw
  END IF
END DO
END DO

DO i = 1,nz
 ustar(i,0) = ustar(i,nx)
 ustar(i,nx+1) = ustar(i,1)
 wstar(i,0) = wstar(i,nx)
 wstar(i,nx+1) = wstar(i,1)
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
DO nsor = 1,nstop
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
 IF(k == 1) dq(i,nx+1) = q2
 IF(k == nx) dq(i,0) = q2
 perr = MAX(ABS(q2-q1),perr)

 END IF

END DO
END DO

! STEP 2: predict new velocities 
DO i = 1,nz
DO k = 1,nx
  IF(wet(i,k))THEN
    pressx = -drdx*(dq(i,k+1)-dq(i,k))
    IF(wet(i,k+1)) un(i,k) = ustar(i,k) + dt*pressx
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
 usum(0) = usum(nx)
 usum(nx+1) = usum(1)

! STEP 3b: predict surface pressure field
DO k = 1,nx
 dq(0,k) = -dt*RHOREF*G*(usum(k)-usum(k-1))/dx
! IF(k == 51) dq(0,k) = dq(0,k)+ad*RHOREF*G
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
  w(i,k) = wn(i,k)
END DO
END DO

DO k = 1,nx
  q(0,k) = q(0,k)+dq(0,k)
END DO

! lateral boundary conditions
 q(0,0) = q(0,nx)
 q(0,nx+1) = q(0,1)

DO i = 1,nz
 u(i,0) = u(i,nx)
 u(i,nx+1) = u(i,1)
 w(i,0) = w(i,nx)
 w(i,nx+1) = w(i,1)
 q(i,0) = q(i,nx)
 q(i,nx+1) = q(i,1)
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

REAL FUNCTION density(tt,ss)

! input parameters
REAL, INTENT(IN) :: tt, ss  

! local parameters
REAL :: alpha, beta, rhozero


alpha = 2.5e-4
beta = 8.0e-4

density = RHOREF*( beta*ss - alpha*tt)

RETURN
END FUNCTION density

END MODULE sub