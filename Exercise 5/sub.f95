MODULE sub
USE param

CONTAINS

!******************
SUBROUTINE init
!local parameters
REAL :: ccc, cmc, depth(nx)
INTEGER :: nb
REAL :: N2

! set initial arrays
DO i = 0,nz+1
DO k = 0,nx+1
 p(i,k) = 0.0
 q(i,k) = 0.0
 dq(i,k) = 0.0
 rho(i,k) = RHOREF
 u(i,k) = 0.0
 un(i,k) = 0.0
 ustar(i,k) = 0.0
 w(i,k) = 0.0
 wn(i,k) = 0.0
 wstar(i,k) = 0.0
 wet(i,k) = .true.
END DO
END DO

DO k = 0,nx+1
 usum(k) = 0.0
END DO

! grid parameters
dx = 5.0
dz = 2.0
dt = 1.0
drdx = 1.0/(RHOREF*dx)
drdz = 1.0/(RHOREF*dz)
kh = 1.e-4
kz = 1.e-4

N2 = 1.e-2 ! stability frequency squared 

! ambient density stratification
rhozero(0) = RHOREF
DO i = 1,nz+1
  rhozero(i) = rhozero(i-1)+N2*RHOREF/G*DZ
END DO

DO k = 0,nx+1
DO i = 0,nz+1
  rho(i,k) = rhozero(i)
END DO
END DO

! forcing region
DO i = 0,nz+1
DO k = 51-2,51+2
 rho(i,k) = rho(i,k)+20.0
 rho(i,k) = AMIN1(rho(i,k),rhozero(nz+1))
END DO
END DO

! closed lateral boundaries
DO i = 0,nz+1
 wet(i,0) = .false.
 wet(i,nx+1) = .false.
END DO

! ambient bathymetry
DO k = 1,nx
 depth(k) = 102.0
END DO

! variable bottom topography
DO k = 1,30
  depth(k) =  60.0
END DO
DO k = 31,40
  depth(k) = 60.0+42.0*REAL(k-30)/10.0
END DO

OPEN(50,file ='h.dat',form='formatted')
  WRITE(50,'(101F12.6)')(depth(k),k=1,nx)
CLOSE(50)

! treatment of dry grid points
DO k = 1,nx
  nb = INT(depth(k)/dz)
  nbot = MIN1(nbot,nz+1)
  DO i = nb,nz+1
    wet(i,k) = .false.
    rho(i,k) = RHOREF
  END DO
  wet(0,k) = .false.
END DO

! coefficients for SOR
omega = 1.4
peps = 1.e-2

DO i = 1,nz
DO k = 1,nx
  at(i,k) = dx/dz
  ab(i,k) = dx/dz
  ae(i,k) = dz/dx
  aw(i,k) = dz/dx
  IF(.not.wet(i,k-1)) aw(i,k) = 0.0
  IF(.not.wet(i,k+1)) ae(i,k) = 0.0
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
REAL :: dif1, dif2, difh, difv, diff

! density prediction
! advection of rho
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
  B(i,k) = rho(i,k)
END DO
END DO

CALL advect

DO i = 1,nz
DO k = 1,nx
 div1 = (u(i,k)-u(i,k-1))/dx
 div2 = (w(i,k)-w(i+1,k))/dz
 div = dt*B(i,k)*(div1+div2)
! diffusion terms
 dif1 = 0.0
 IF(WET(i,k+1)) dif1 = (rho(i,k+1)-rho(i,k))/dx
 dif2 = 0.0
 IF(WET(i,k-1)) dif2 = (rho(i,k)-rho(i,k-1))/dx
 difh = kh*(dif2-dif1)/dx
 dif1 = 0.0
 IF(WET(i-1,k)) dif1 = (rho(i-1,k)-rho(i,k))/dz
 dif2 = 0.0
 IF(WET(i+1,k)) dif2 = (rho(i,k)-rho(i+1,k))/dz
 difz = kz*(dif2-dif1)/dz
 diff = dt*(difh+difz)
! prognostic equation for density  
 rhon(i,k)= rho(i,k)+BN(i,k)+div+diff
END DO
END DO

! boundary conditions
DO i = 1,nz
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
    P(i,k) = P(i-1,k) + (0.5*(rho(i-1,k)+rho(i,k))-rhoref)*G*dz
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
  IF(wet(i,k))THEN
    pressx = -drdx*(q(i,k+1)-q(i,k))-drdx*(p(i,k+1)-p(i,k))
    IF(wet(i,k+1)) ustar(i,k) = u(i,k) + dt*pressx + advx(i,k)
    pressz = -drdz*(q(i-1,k)-q(i,k))
    IF(wet(i-1,k)) wstar(i,k) = w(i,k) + dt*pressz + advz(i,k) 
  END IF
END DO
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
usum(nx) = 0.0

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
  w(i,k) = wn(i,k)
END DO
END DO

DO k = 1,nx
  q(0,k) = q(0,k)+dq(0,k)
END DO

! lateral boundary conditions
DO i = 1,nz
 u(i,nx) = 0.0
 q(i,nx+1) = q(i,nx)
 q(i,0) = q(i,1)
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