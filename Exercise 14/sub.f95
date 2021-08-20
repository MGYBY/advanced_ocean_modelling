MODULE sub
USE param

CONTAINS

!******************
SUBROUTINE init
!local parameters
REAL :: ccc, cmc, depth(0:nx+1),uu
INTEGER :: nb
REAL :: N2

! set initial arrays
DO i = 0,nz+1
DO k = 0,nx+1
 p(i,k) = 0.0
 q(i,k) = 0.0
 dq(i,k) = 0.0
 rho(i,k) = 0.0
 u(i,k) = 0.0
 un(i,k) = 0.0
 c1(i,k) = 0.0 ! used for age calculation
 c1n(i,k) = 0.0
 c2(i,k) = 0.0
 c2n(i,k) = 0.0
 ustar(i,k) = 0.0
 w(i,k) = 0.0
 wn(i,k) = 0.0
 wstar(i,k) = 0.0
 wet(i,k) = .true.
 az(i,k) = 0.0
 kz(i,k) = 0.0
 uM(i,k) = 0.0
 wM(i,k) = 0.0
 rhoM(i,k) = 0.0
END DO
END DO

DO k = 0,nx+1
 usum(k) = 0.0
 etaMIN(k) = 0.0 ! used to determine tidal range
 etaMAX(k) = 0.0
END DO

! grid parameters
dx = 2000.0
dz = 1.0

! channel width
DO k = 0,60
  Width(k) = 50.0
END DO
DO k = 61,nx
  Width(k) = 50.0+REAL(k-60)/REAL(nx-60)*450.0
END DO

Width(nx+1) = Width(nx)

dt = 30.0
drdx = 1.0/(RHOREF*dx)
drdz = 1.0/(RHOREF*dz)

kh = 0.1
ah = 0.1
r = 1.e-3

! channel depth
DO k = 1,101
 depth(k) = 5.0 + 15.0*REAL(k-1)/100.0
END DO

DO k = 102,nx
 depth(k) = 20.0
END DO

depth(0) = depth(1)
depth(nx+1) = depth(nx)

OPEN(50,file ='h.dat',form='formatted')
  WRITE(50,'(201F12.6)')(depth(k),k=1,nx)
CLOSE(50)

DO k = 0,nx+1
  nb = INT(depth(k)/dz)
  nb = MIN(nb,nz)
  DO i = nb+1,nz+1
    wet(i,k) = .false.
  END DO
! enable this for the rigid-lid approximation
! wet(0,k) = .false.
END DO

DO k = 0,nx+1
DO i = 0,nz+1
  dry(i,k)=.NOT.wet(i,k)
END DO
END DO

DO i = 0,nz+1
 IF(wet(i,0))THEN
  rho(i,0) = 0.0
  c1(i,0) = 0.0
  c1(i,nx+1) = 0.0
 END IF
DO k = 1,nx+1
  IF(wet(i,k))rho(i,k) = 28.0
END DO
END DO

DO k = 0,nx+1
DO i = 0,nz+1
 rhon(i,k) = rho(i,k)
END DO
END DO

! coefficients for SOR
omega = 1.4
peps = 1.e-3

DO k = 1,nx
DO i = 1,nz
  at(i,k) = dx/dz*Width(k)
  ab(i,k) = dx/dz*Width(k)
  ae(i,k) = dz/dx*0.5*(Width(k+1)+Width(k))
  aw(i,k) = dz/dx*0.5*(Width(k-1)+Width(k))

  IF(dry(i,k+1)) ae(i,k) = 0.0
  IF(dry(i,k-1)) aw(i,k) = 0.0
  IF(dry(i+1,k)) ab(i,k) = 0.0
  IF(dry(i-1,k)) at(i,k) = 0.0
  atot(i,k) = ab(i,k)+at(i,k)+ae(i,k)+aw(i,k)
END DO
END DO

! index used for calculation of bottom friction term
DO i = 0,nz
DO k = 0,nx+1
 ibu(k) = 0
END DO
END DO

DO i = 0,nz
DO k = 0,nx+1
 IF(wet(i,k).AND.dry(i+1,k))THEN
   ibu(k) = i
 END IF
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
REAL :: div, div1, div2, uum, unn
REAL :: dif1, dif2, difh, diff, diffu, diffw
REAL :: drho, speed, aztop, azbot, ppp, dummy
REAL :: w1,w2, pgradx
REAL :: we, ww, wc

call eddy ! turbulence closure scheme

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
 rho(i,0) = 0.0 ! boundary condition for river water
 IF(u(i,nx+1)<0) rho(i,nx+1) = 28.0 ! inflow condition for seawater
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
 we = 0.5*(width(k)+width(k+1))
 ww = 0.5*(width(k)+width(k-1))
 wc = width(k)
 IF(WET(i,k+1)) dif1 = we*(rho(i,k+1)-rho(i,k))/dx
 IF(k==nx) dif1 = 0.0
 dif2 = 0.0
 IF(WET(i,k-1)) dif2 = (ww*rho(i,k)-rho(i,k-1))/dx
 if(k==1) dif2 = 0.0
 difh = kh*(dif1-dif2)/dx/wc

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
  rhon(i,0) = 0.0
  rhon(i,nx+1) = rhon(i,nx)
END DO

DO k = 0,nx+1
 rhon(0,k) = rhon(1,k)
 rhon(nz+1,k) = rhon(nz,k)
END DO

! updating
DO k = 0,nx+1
DO i = 0,nz+1
 IF(wet(i,k))THEN
   rho(i,k) = rhon(i,k)
   rhoM(i,k) = rhoM(i,k)+dt*rho(i,k) ! integrate to determine tidal average value
 END IF
END DO
END DO

DO k = 0,nx+1
  c1(0,k) = 0.0
  c1(1,k) = 0.0 
END DO

! water age prediction
DO i = 0,nz+1
DO k = 61,nx+1
  c1(i,k) = 0.0
END DO
END DO

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
 we = 0.5*(width(k)+width(k+1))
 ww = 0.5*(width(k)+width(k-1))
 wc = width(k)
 IF(WET(i,k+1)) dif1 = we*(c1(i,k+1)-c1(i,k))/dx
 dif2 = 0.0
 IF(WET(i,k-1)) dif2 = ww*(c1(i,k)-c1(i,k-1))/dx
 difh = kh*(dif1-dif2)/dx/wc

 aztop = 0.5*(kz(i,k)+kz(i-1,k))
 dif1 = aztop*(c1(i-1,k)-c1(i,k))/dz
 IF(dry(i-1,k)) dif1 = 0.0 
 azbot = 0.5*(kz(i,k)+kz(i+1,k))
 dif2 = azbot*(c1(i,k)-c1(i+1,k))/dz
 IF(dry(i+1,k)) dif2 = 0.0 
 difz = (dif1-dif2)/dz

 diff = dt*(difh+difz)
  
 c1n(i,k)= c1(i,k)+BN(i,k)+div+diff + dt
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

! calculate bottom shear-stress components
DO k = 1,nx
  uum = u(ibu(k),k)
  speed = SQRT(uum*uum)
  tbu(k) = r*uum*speed
END DO

! no tangential surface stress
DO k = 1,nx
  u(0,k) = u(1,k)
END DO

! forcing

! tidal forcing at ocean boundary
q(0,1) = -ad*(0.25*SIN(2.*PI*(time+1800.)/period)-0.5)*RHOREF*G
q(0,0) = q(0,1)
q(0,nx) = -0.25*SIN(2.*PI*time/period)*RHOREF*G
q(0,nx+1) = q(0,nx)

! calculate ustar and wstar 
DO i = 1,nz
DO k = 1,nx

! diffusion of u
 dif1 = 0.0
 we = width(k+1)
 ww = width(k)
 wc = 0.5*(width(k)+width(k+1))

 IF(WET(i,k+1)) dif1 = we*(u(i,k+1)-u(i,k))/dx
 dif2 = 0.0
 IF(WET(i,k-1)) dif2 = ww*(u(i,k)-u(i,k-1))/dx
 difh = ah*(dif1-dif2)/dx/wc

 aztop = 0.25*(az(i,k)+az(i-1,k)+az(i,k+1)+az(i-1,k+1))
 dif1 = aztop*(u(i-1,k)-u(i,k))/dz
 IF(dry(i-1,k)) dif1 = 0.0 
 azbot = 0.25*(az(i,k)+az(i+1,k)+az(i,k+1)+az(i+1,k+1))
 dif2 = azbot*(u(i,k)-u(i+1,k))/dz
 IF(i == ibu(k) ) dif2 = tbu(k)
 difz = (dif1-dif2)/dz

 diffu = dt*(difh+difz)  

! diffusion of w
 dif1 = 0.0
 we = 0.5*(width(k)+width(k+1))
 ww = 0.5*(width(k)+width(k-1))
 wc = width(k)
 IF(WET(i,k+1)) dif1 = we*(w(i,k+1)-w(i,k))/dx
 dif2 = 0.0
 IF(WET(i,k-1)) dif2 = ww*(w(i,k)-w(i,k-1))/dx
 difh = ah*(dif1-dif2)/dx/wc

 aztop = az(i-1,k)
 dif1 = aztop*(w(i-1,k)-w(i,k))/dz
 IF(dry(i-1,k)) dif1 = 0.0
 azbot = az(i,k)
 dif2 = azbot*(w(i,k)-w(i+1,k))/dz
 IF(dry(i+1,k)) dif2 = 0.0 
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

! boundary conditions
DO i = 1,nz
 ustar(i,0) = ustar(i,1)
 ustar(i,nx+1) = ustar(i,nx)
 wstar(i,0) = wstar(i,1)
 wstar(i,nx+1) = wstar(i,nx) 
END DO

! calculate right-hand side of Poisson equation
DO i = 1,nz
DO k = 1,nx
 we = 0.5*(Width(k)+Width(k+1))
 ww = 0.5*(Width(k)+Width(k-1))
 wc = Width(k)
 qstar(i,k) = -1.*RHOREF/dt*(  &
 &  (we*ustar(i,k)-ww*ustar(i,k-1))*dz + &
 &  (wstar(i,k)-wstar(i+1,k))*wc*dx )
END DO
END DO

! STEP 4: S.O.R. ITERATION

nstop = 10000

DO i = 0,nz+1
DO k = 0,nx+1
  dq(i,k) = 0.0
END DO
END DO

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
 perr = AMAX1(ABS(q2-q1),perr)

 END IF

END DO
END DO

! STEP 2: predict new velocities 
DO i = 1,nz
DO k = 1,nx
  IF(wet(i,k))THEN
    pressx = -drdx*( dq(i,k+1)-dq(i,k))
    IF(wet(i,k+1))un(i,k) = ustar(i,k) + dt*pressx
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
 usum(0) = usum(1) 
 usum(nx+1) = usum(nx)

! STEP 3b: predict surface pressure field
DO k = 1,nx
 w1 = 0.5*(Width(k)+Width(k+1))
 w2 = 0.5*(Width(k-1)+Width(k))
 dq(0,k) = -dt*RHOREF*G*(usum(k)*w1-usum(k-1)*w2)/dx
 dq(0,k) = dq(0,k)/Width(k)
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
  uM(i,k) = uM(i,k)+dt*u(i,k)
  wM(i,k) = wM(i,k)+dt*w(i,k)
END DO
END DO

DO k = 1,nx
  q(0,k) = q(0,k)+dq(0,k)
  dummy = q(0,k)/(RHOREF*G)
  IF(dummy>etaMAX(k))etaMAX(k) = dummy ! for calculation of tidal range
  IF(dummy<etaMIN(k))etaMIN(k) = dummy
END DO

! lateral boundary conditions
q(0,0) = q(0,1)
q(0,nx+1) = q(0,nx)

DO i = 1,nz
 u(i,0) = u(i,1)
 u(i,nx+1) = u(i,nx) 
 w(i,0) = w(i,1)
 w(i,nx+1) = w(i,nx)
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

SUBROUTINE eddy
! local parameters
REAL, PARAMETER :: nn = 2.0
REAL, PARAMETER :: bb = 5.0
REAL, PARAMETER :: nu0 = 5.e-3
REAL, PARAMETER :: nub = 1.e-4
REAL, PARAMETER :: kb = 1.e-5
REAL :: nu, ff
REAL :: utop, ubot, dudz
REAL :: eddyx, stab, wurz
REAL :: Ri ! Richardson number

DO k = 1,nx
DO i = 1,nz
  utop = 0.5*(u(i-1,k)+u(i-1,k-1))
  ubot = 0.5*(u(i+1,k)+u(i+1,k-1))
  dudz = (utop-ubot)/(2.0*dz)
  eddyx = dudz*dudz
  stab = -G/RHOREF*(rho(i-1,k)-rho(i+1,k))/(2.0*dz)
! Gradient Richardson number
  RI = 100.0 ! set initial value
  IF(eddyx>1.e-6) Ri = stab/eddyx
  IF(stab<0.0) Ri = 0.0 ! parameterisation of convection
  ff = 1.0+bb*Ri
  wurz = nu0/(ff**nn)+nub
!****** final value ******
  az(i,k) = wurz
  kz(i,k) = wurz/ff+kb
!*************************
END DO
END DO

DO k = 1,nx
  kz(0,k) = kz(1,k)
  kz(nz+1,k) = kz(nz,k)
  az(0,k) = az(1,k)
  az(nz+1,k) = az(nz,k)
END DO

DO i = 0,nz+1
  kz(i,0) = kz(i,1)
  kz(i,nx+1) = kz(i,nx)
  az(i,0) = az(i,1)
  az(i,nx+1) = az(i,nx)
END DO

END SUBROUTINE eddy

END MODULE sub