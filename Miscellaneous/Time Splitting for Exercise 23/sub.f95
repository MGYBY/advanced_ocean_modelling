MODULE sub
USE param

CONTAINS

!******************
SUBROUTINE init
!local parameters
INTEGER :: nb

! set initial arrays
DO i = 0,nz+1
DO j = 0,ny+1
DO k = 0,nx+1
 p(i,j,k) = 0.0
 q(i,j,k) = 0.0
 dq(i,j,k) = 0.0
 rho(i,j,k) = 0.0
 c(i,j,k) = 0.0
 cn(i,j,k) = 0.0
 u(i,j,k) = 0.0
 un(i,j,k) = 0.0
 ustar(i,j,k) = 0.0
 v(i,j,k) = 0.0
 vn(i,j,k) = 0.0
 vstar(i,j,k) = 0.0
 w(i,j,k) = 0.0
 wn(i,j,k) = 0.0
 wstar(i,j,k) = 0.0
 wet(i,j,k) = .true.
 ah(i,j,k) = 1.0
 kh(i,j,k) = 1.0
 az(i,j,k) = 0.0
 kz(i,j,k) = 0.0
END DO
END DO
END DO

Do j = 0,ny+1
DO k = 0,nx+1
 usum(j,k) = 0.0
 vsum(j,k) = 0.0
END DO
END DO

! grid parameters
dx = 2000.0
dy = 2000.0
dz = 10.0
dte = 30.0 ! time-step for fast propagating surface gravity waves
nfrac = 20
dti = REAL(nfrac)*dte ! time step for ocean interior

drdx = 1.0/(RHOREF*dx)
drdy = 1.0/(RHOREF*dy)
drdz = 1.0/(RHOREF*dz)
f = 1.e-4
alpha = dti*f
r = 1.e-3

OPEN(50,file ='topo.dat',form='formatted')
DO j = 0,ny+1
  READ(50,'(103F12.6)')(depth(j,k),k=0,nx+1)
END DO
CLOSE(50)

OPEN(50,file ='h.dat',form='formatted')
DO j = 1,ny
  WRITE(50,'(101F12.6)')(depth(j,k),k=1,nx)
END DO
CLOSE(50)


DO k = 0,nx+1
DO j = 0,ny+1
  nb = INT(depth(j,k)/dz)
  ib(j,k) = MIN(nb,nz)
  DO i = ib(j,k)+1,nz+1
    wet(i,j,k) = .false.
  END DO
! wet(0,j,k) = .false. ! use this for rigid-lid condition
END DO
END DO

DO k = 0,nx+1
DO j = 0,ny+1
DO i = 0,nz+1
  dry(i,j,k)=.NOT.wet(i,j,k)
END DO
END DO
END DO

DO k = 0,nx+1
DO j = 0,ny+1
DO i = 5,nz+1
  IF(wet(i,j,k))THEN
     c(i,j,k) = 1.0
     rho(i,j,k) = 1.0
  END IF
END DO
END DO
END DO


! coefficients for SOR
omega = 1.4
peps = 1.e-2

DO k = 1,nx
DO j = 1,ny
DO i = 1,nz
  at(i,j,k) = dx/dz
  ab(i,j,k) = dx/dz
  ae(i,j,k) = dz/dx
  aw(i,j,k) = dz/dx
  an(i,j,k) = dz*dx/(dy*dy)
  as(i,j,k) = dz*dx/(dy*dy)
  IF(dry(i,j,k+1).OR.k == nx) ae(i,j,k) = 0.0
  IF(dry(i,j,k-1).OR.k == 1) aw(i,j,k) = 0.0
  IF(dry(i,j+1,k).OR.j == ny) an(i,j,k) = 0.0
  IF(dry(i,j-1,k).OR.j == 1) as(i,j,k) = 0.0
  IF(dry(i+1,j,k)) ab(i,j,k) = 0.0
  IF(dry(i-1,j,k)) at(i,j,k) = 0.0
  atot(i,j,k) = ab(i,j,k)+at(i,j,k)+ae(i,j,k)+aw(i,j,k)+an(i,j,k)+as(i,j,k)
END DO
END DO
END DO

END SUBROUTINE init

!*******************
SUBROUTINE dyn

! local parameters
REAL :: pressx, pressy, pressz
REAL :: perr, q1, q2, term1
INTEGER :: nsor, nstop
REAL :: advx(0:nz+1,0:ny+1,0:nx+1), advy(0:nz+1,0:ny+1,0:nx+1), advz(0:nz+1,0:ny+1,0:nx+1)
REAL :: du(0:nz+1,0:ny+1,0:nx+1), dv(0:nz+1,0:ny+1,0:nx+1), dw(0:nz+1,0:ny+1,0:nx+1)  
REAL :: div, div1, div2, div3, um, vm
REAL :: dif1, dif2, dif3, dif4, difh, difv, diff, diffu, diffv, diffw
REAL :: drho, speed
REAL :: azt, azb, ahe, ahw, ahn, ahs
REAL :: kzt, kzb, khe, khw, khn, khs
INTEGER :: nrep

call eddy ! turbulence closure scheme

! density prediction
DO k = 0,nx+1
DO j = 0,ny+1
DO i = 0,nz+1
  CuP(i,j,k) = 0.5*(u(i,j,k)+abs(u(i,j,k)))*dti/dx
  CuN(i,j,k) = 0.5*(u(i,j,k)-abs(u(i,j,k)))*dti/dx
  CvP(i,j,k) = 0.5*(v(i,j,k)+abs(v(i,j,k)))*dti/dy
  CvN(i,j,k) = 0.5*(v(i,j,k)-abs(v(i,j,k)))*dti/dy
  CwP(i,j,k) = 0.5*(w(i,j,k)+abs(w(i,j,k)))*dti/dz
  CwN(i,j,k) = 0.5*(w(i,j,k)-abs(w(i,j,k)))*dti/dz
END DO
END DO
END DO

DO i = 0,nz+1
DO j = 0,ny+1
DO k = 0,nx+1
  B(i,j,k) = rho(i,j,k)
END DO
END DO
END DO

CALL advect

DO k = 1,nx
DO j = 1,ny
DO i = 1,nz
 div1 = (u(i,j,k)-u(i,j,k-1))/dx
 div2 = (v(i,j,k)-v(i,j-1,k))/dy
 div3 = (w(i,j,k)-w(i+1,j,k))/dz
 div = dti*B(i,j,k)*(div1+div2+div3)
 khe = 0.5*(ah(i,j,k)+ah(i,j,k+1))
 dif1 = khe*(rho(i,j,k+1)-rho(i,j,k))/dx
 IF(dry(i,j,k+1)) dif1 = 0.0
 khw = 0.5*(ah(i,j,k)+ah(i,j,k-1))
 dif2 = khw*(rho(i,j,k)-rho(i,j,k-1))/dx
 IF(dry(i,j,k-1)) dif2 = 0.0
 khn = 0.5*(ah(i,j,k)+ah(i,j+1,k))
 dif3 = khn*(rho(i,j+1,k)-rho(i,j,k))/dy
 IF(dry(i,j+1,k)) dif3 = 0.0
 khs = 0.5*(ah(i,j,k)+ah(i,j-1,k))
 dif4 = khs*(rho(i,j,k)-rho(i,j-1,k))/dy
 IF(dry(i,j-1,k)) dif4 = 0.0

 difh = (dif1-dif2)/dx + (dif3-dif4)/dy 

 azt = 0.5*(kz(i,j,k)+kz(i-1,j,k))
 dif1 = azt*(rho(i-1,j,k)-rho(i,j,k))/dz
 IF(dry(i-1,j,k)) dif1 = 0.0
 azb = 0.5*(kz(i,j,k)+kz(i+1,j,k))
 dif2 = azb*(rho(i,j,k)-rho(i+1,j,k))/dz
 IF(dry(i+1,j,k)) dif2 = 0.0
 difz = (dif1-dif2)/dz
 diff = dti*(difh+difz)  
 rhon(i,j,k)= rho(i,j,k)+BN(i,j,k)+div+diff
END DO
END DO
END DO

! boundary conditions
DO i = 0,nz+1
DO j = 0,ny+1
 rhon(i,j,0) = rhon(i,j,1)  ! keep upstream boundary value at initial value
 rhon(i,j,nx+1) = rhon(i,j,nx)
END DO
END DO

DO i = 0,nz+1
DO k = 0,nx+1
 rhon(i,0,k) = rhon(i,1,k)
 rhon(i,ny+1,k) = rhon(i,ny,k)
END DO
END DO

DO j = 0,ny+1
DO k = 0,nx+1
 rhon(0,j,k) = rhon(1,j,k)
 rhon(nz+1,j,k) = rhon(nz,j,k)
END DO
END DO

! updating

DO k = 0,nx+1
DO j = 0,ny+1
DO i = 0,nz+1
 rho(i,j,k) = rhon(i,j,k)
END DO
END DO
END DO

! calculate hydrostatic pressure

DO j = 0,ny+1
DO k = 0,nx+1
  P(0,j,k) = 0.0
  DO i = 1,nz+1
    P(i,j,k) = P(i-1,j,k) + 0.5*(rho(i-1,j,k)+rho(i,j,k))*G*dz
  END DO
END DO
END DO

! tracer prediction
DO i = 0,nz+1
DO j = 0,ny+1
DO k = 0,nx+1
  B(i,j,k) = c(i,j,k)
END DO
END DO
END DO

CALL advect

DO k = 1,nx
DO j = 1,ny
DO i = 1,nz
 div1 = (u(i,j,k)-u(i,j,k-1))/dx
 div2 = (v(i,j,k)-v(i,j-1,k))/dy
 div3 = (w(i,j,k)-w(i+1,j,k))/dz
 div = dti*B(i,j,k)*(div1+div2+div3)
 khe = 0.5*(ah(i,j,k)+ah(i,j,k+1))
 dif1 = khe*(rho(i,j,k+1)-rho(i,j,k))/dx
 IF(dry(i,j,k+1)) dif1 = 0.0
 khw = 0.5*(ah(i,j,k)+ah(i,j,k-1))
 dif2 = khw*(c(i,j,k)-c(i,j,k-1))/dx
 IF(dry(i,j,k-1)) dif2 = 0.0
 khn = 0.5*(ah(i,j,k)+ah(i,j+1,k))
 dif3 = khn*(c(i,j+1,k)-c(i,j,k))/dy
 IF(dry(i,j+1,k)) dif3 = 0.0
 khs = 0.5*(ah(i,j,k)+ah(i,j-1,k))
 dif4 = khs*(c(i,j,k)-c(i,j-1,k))/dy
 IF(dry(i,j-1,k)) dif4 = 0.0

 difh = (dif1-dif2)/dx + (dif3-dif4)/dy 

 azt = 0.5*(kz(i,j,k)+kz(i-1,j,k))
 dif1 = azt*(c(i-1,j,k)-c(i,j,k))/dz
 IF(dry(i-1,j,k)) dif1 = 0.0
 azb = 0.5*(kz(i,j,k)+kz(i+1,j,k))
 dif2 = azb*(c(i,j,k)-c(i+1,j,k))/dz
 IF(dry(i+1,j,k)) dif2 = 0.0
 difz = (dif1-dif2)/dz
 diff = dti*(difh+difz)  
 cn(i,j,k)= c(i,j,k)+BN(i,j,k)+div+diff
END DO
END DO
END DO

! boundary conditions
DO i = 0,nz+1
DO j = 0,ny+1
 cn(i,j,0) = cn(i,j,1)
 cn(i,j,nx+1) = cn(i,j,nx)
END DO
END DO

DO i = 0,nz+1
DO k = 0,nx+1
 cn(i,0,k) = cn(i,1,k)
 cn(i,ny+1,k) = cn(i,ny,k)
END DO
END DO

DO j = 0,ny+1
DO k = 0,nx+1
 cn(0,j,k) = cn(1,j,k)
 cn(nz+1,j,k) = cn(nz,j,k)
END DO
END DO

DO k = 1,nx
DO j = 1,ny
DO i = 0,nz+1

IF(wet(i,j,k))THEN

IF(dry(i,j+1,k)) cn(i,j+1,k) = cn(i,j,k)
IF(dry(i,j-1,k)) cn(i,j-1,k) = cn(i,j,k)
IF(dry(i,j,k+1)) cn(i,j,k+1) = cn(i,j,k)
IF(dry(i,j,k-1)) cn(i,j,k-1) = cn(i,j,k)

END IF

END DO
END DO
END DO

! updating
DO k = 0,nx+1
DO j = 0,ny+1
DO i = 0,nz+1
 c(i,j,k) = cn(i,j,k)
END DO
END DO
END DO

! calculate the nonlinear terms for u-momentum equation
DO i = 0,nz+1
DO j = 0,ny
DO k = 0,nx
  CuP(i,j,k) = 0.25*(u(i,j,k)+u(i,j,k+1)+abs(u(i,j,k))+abs(u(i,j,k+1)))*dti/dx
  CuN(i,j,k) = 0.25*(u(i,j,k)+u(i,j,k+1)-abs(u(i,j,k))-abs(u(i,j,k+1)))*dti/dx
  CvP(i,j,k) = 0.25*(v(i,j,k)+v(i,j,k+1)+abs(v(i,j,k))+abs(v(i,j,k+1)))*dti/dy
  CvN(i,j,k) = 0.25*(v(i,j,k)+v(i,j,k+1)-abs(v(i,j,k))-abs(v(i,j,k+1)))*dti/dy
  CwP(i,j,k) = 0.25*(w(i,j,k)+w(i,j,k+1)+abs(w(i,j,k))+abs(w(i,j,k+1)))*dti/dz
  CwN(i,j,k) = 0.25*(w(i,j,k)+w(i,j,k+1)-abs(w(i,j,k))-abs(w(i,j,k+1)))*dti/dz
END DO
END DO
END DO

DO i = 0,nz+1
DO j = 0,ny+1
DO k = 0,nx+1
  B(i,j,k) = u(i,j,k)
END DO
END DO
END DO

CALL advect

DO i = 1,nz
DO j = 1,ny
DO k = 1,nx
  div1 = 0.5*(u(i,j,k+1)-u(i,j,k-1))/dx
  div2 = 0.5*(v(i,j,k)+v(i,j,k+1)-v(i,j-1,k)-v(i,j-1,k+1))/dy
  div3 = 0.5*(w(i,j,k)+w(i,j,k+1)-w(i+1,j,k)-w(i+1,j,k+1))/dz
  div = dti*B(i,j,k)*(div1+div2+div3)  
  advx(i,j,k)= BN(i,j,k)+div
END DO
END DO
END DO

! calculate the nonlinear terms for v-momentum equation

DO i = 0,nz+1
DO j = 0,ny
DO k = 0,nx
  CuP(i,j,k) = 0.25*(u(i,j,k)+u(i,j+1,k)+abs(u(i,j,k))+abs(u(i,j+1,k)))*dti/dx
  CuN(i,j,k) = 0.25*(u(i,j,k)+u(i,j+1,k)-abs(u(i,j,k))-abs(u(i,j+1,k)))*dti/dx
  CvP(i,j,k) = 0.25*(v(i,j,k)+v(i,j+1,k)+abs(v(i,j,k))+abs(v(i,j+1,k)))*dti/dy
  CvN(i,j,k) = 0.25*(v(i,j,k)+v(i,j+1,k)-abs(v(i,j,k))-abs(v(i,j+1,k)))*dti/dy
  CwP(i,j,k) = 0.25*(w(i,j,k)+w(i,j+1,k)+abs(w(i,j,k))+abs(w(i,j+1,k)))*dti/dz
  CwN(i,j,k) = 0.25*(w(i,j,k)+w(i,j+1,k)-abs(w(i,j,k))-abs(w(i,j+1,k)))*dti/dz
END DO
END DO
END DO

DO i = 0,nz+1
DO j = 0,ny+1
DO k = 0,nx+1
  B(i,j,k) = v(i,j,k)
END DO
END DO
END DO

CALL advect

DO i = 1,nz
DO j = 1,ny
DO k = 1,nx
  div1 = 0.5*(u(i,j,k)+u(i,j+1,k)-u(i,j,k-1)-u(i,j+1,k-1))/dx
  div2 = 0.5*(v(i,j+1,k)-v(i,j-1,k))/dy
  div3 = 0.5*(w(i,j,k)+w(i,j+1,k)-w(i+1,j,k)-w(i+1,j+1,k))/dz
  div = dti*B(i,j,k)*(div1+div2+div3)  
  advy(i,j,k)= BN(i,j,k)+div
END DO
END DO
END DO

! calculate the nonlinear terms for w-momentum equation
DO i = 0,nz
DO j = 0,ny+1
DO k = 0,nx+1
  CuP(i,j,k) = 0.25*(u(i,j,k)+u(i-1,j,k)+abs(u(i,j,k))+abs(u(i-1,j,k)))*dti/dx
  CuN(i,j,k) = 0.25*(u(i,j,k)+u(i-1,j,k)-abs(u(i,j,k))-abs(u(i-1,j,k)))*dti/dx
  CvP(i,j,k) = 0.25*(v(i,j,k)+v(i-1,j,k)+abs(v(i,j,k))+abs(v(i-1,j,k)))*dti/dy
  CvN(i,j,k) = 0.25*(v(i,j,k)+v(i-1,j,k)-abs(v(i,j,k))-abs(v(i-1,j,k)))*dti/dy
  CwP(i,j,k) = 0.25*(w(i,j,k)+w(i-1,j,k)+abs(w(i,j,k))+abs(w(i-1,j,k)))*dti/dz
  CwN(i,j,k) = 0.25*(w(i,j,k)+w(i-1,j,k)-abs(w(i,j,k))-abs(w(i-1,j,k)))*dti/dz
END DO
END DO
END DO

DO i = 0,nz+1
DO j = 0,ny+1
DO k = 0,nx+1
  B(i,j,k) = w(i,j,k)
END DO
END DO
END DO

CALL advect

DO i = 1,nz
DO j = 1,ny
DO k = 1,nx
  div1 = 0.5*(u(i,j,k)+u(i-1,j,k)-u(i,j,k-1)-u(i-1,j,k-1))/dx
  div2 = 0.5*(v(i,j,k)+v(i-1,j,k)-v(i,j-1,k)-v(i-1,j-1,k))/dy
  div3 = 0.5*(w(i-1,j,k)-w(i+1,j,k))/dz
  div = dti*B(i,j,k)*(div1+div2+div3)  
  advz(i,j,k)= BN(i,j,k) + div
END DO
END DO
END DO

! calculate bottom shear-stress components

DO k = 1,nx
DO j = 1,ny
  um = u(ib(j,k),j,k)
  vm = 0.25*(v(ib(j,k),j,k)+v(ib(j,k),j,k+1)+ &
& v(ib(j,k),j-1,k)+v(ib(j,k),j-1,k+1))
  speed = SQRT(um*um+vm*vm)
  tbu(j,k) = r*um*speed
  vm = v(ib(j,k),j,k)
  um = 0.25*(u(ib(j,k),j,k)+u(ib(j,k),j,k-1)+ &
& u(ib(j,k),j+1,k)+u(ib(j,k),j+1,k-1))  
  speed = SQRT(um*um+vm*vm)
  tbv(j,k) = r*vm*speed
END DO
END DO

! wind forcing
taux = 0.1

! surface boundary conditions
DO j = 1,ny
DO k = 1,nx
  azt = 0.25*(az(1,j,k)+az(0,j,k)+az(1,j,k+1)+az(0,j,k+1))
  u(0,j,k) = u(1,j,k)+ad*taux*dz/RHOREF/azt
  v(0,j,k) = v(1,j,k)
END DO
END DO

! calculate du, dv and dw
DO i = 1,nz
DO j = 1,ny
DO k = 1,nx
  du(i,j,k) = 0.0
  dv(i,j,k) = 0.0
  dw(i,j,k) = 0.0
END DO
END DO
END DO
 
DO i = 1,nz
DO j = 1,ny
DO k = 1,nx

! diffusion of u
 ahe = ah(i,j,k+1)
 dif1 = ahe*(u(i,j,k+1)-u(i,j,k))/dx
 IF(dry(i,j,k+1)) dif1 = 0.0

 ahw = ah(i,j,k)
 dif2 = ahw*(u(i,j,k)-u(i,j,k-1))/dx
 IF(dry(i,j,k-1)) dif2 = 0.0

 ahn = 0.25*(ah(i,j,k)+ah(i,j,k+1)+ah(i,j+1,k)+ah(i,j+1,k+1))
 dif3 = ahn*(u(i,j+1,k)-u(i,j,k))/dy
 IF(dry(i,j+1,k)) dif3 = 0.0
 IF(depth(j+1,k)<=0.0) dif3 = -2.0*ahn*u(i,j,k)/dy ! no-slip condition

 ahs = 0.25*(ah(i,j,k)+ah(i,j,k+1)+ah(i,j-1,k)+ah(i,j-1,k+1))
 dif4 = ahs*(u(i,j,k)-u(i,j-1,k))/dy
 IF(dry(i,j-1,k)) dif4 = 0.0
 IF(depth(j-1,k)<=0.0) dif4 = 2.0*ahs*u(i,j,k)/dy ! no-slip condition

 azt = 0.25*(az(i,j,k)+az(i-1,j,k)+az(i,j,k+1)+az(i-1,j,k+1))
 dif1 = azt*(u(i-1,j,k)-u(i,j,k))/dz
 IF(dry(i-1,j,k)) dif1 = 0.0 
 azb = 0.25*(az(i,j,k)+az(i+1,j,k)+az(i,j,k+1)+az(i+1,j,k+1))
 dif2 = azb*(u(i,j,k)-u(i+1,j,k))/dz
 IF(i == ib(j,k) ) dif2 = tbu(j,k)
 difz = (dif1-dif2)/dz
 diffu = dti*(difh+difz)  

! diffusion of v
 ahe = 0.25*(ah(i,j,k)+ah(i,j+1,k)+ah(i,j,k+1)+ah(i,j+1,k+1))
 dif1 = ahe*(v(i,j,k+1)-v(i,j,k))/dx
 IF(dry(i,j,k+1)) dif1 = 0.0
 IF(depth(j,k+1)<=0.0) dif1 = -2.0*ahn*v(i,j,k)/dy ! no-slip condition

 ahw = 0.25*(ah(i,j,k)+ah(i,j+1,k)+ah(i,j,k-1)+ah(i,j+1,k-1))
 dif2 = ahw*(v(i,j,k)-v(i,j,k-1))/dx
 IF(dry(i,j,k-1)) dif2 = 0.0
 IF(depth(j,k-1)<=0.0) dif2 = 2.0*ahn*v(i,j,k)/dy ! no-slip condition

 ahn = ah(i,j+1,k)
 dif3 = ahn*(v(i,j+1,k)-v(i,j,k))/dy
 IF(dry(i,j+1,k)) dif3 = 0.0
 ahs = ah(i,j,k)
 dif4 = ahs*(v(i,j,k)-v(i,j-1,k))/dx
 IF(dry(i,j-1,k)) dif4 = 0.0
 difh = (dif1-dif2)/dx + (dif3-dif4)/dy 

 azt = 0.25*(az(i,j,k)+az(i,j+1,k)+az(i-1,j,k)+az(i-1,j+1,k))
 dif1 = azt*(v(i-1,j,k)-v(i,j,k))/dz
 IF(dry(i-1,j,k)) dif1 = 0.0 
 azb = 0.25*(az(i,j,k)+az(i,j+1,k)+az(i+1,j,k)+az(i+1,j+1,k))
 dif2 = azb*(v(i,j,k)-v(i+1,j,k))/dz
 IF(i == ib(j,k) ) dif2 = tbv(j,k)
 difz = (dif1-dif2)/dz
 diffv = dti*(difh+difz)  

! diffusion of w
 ahe = 0.25*(ah(i,j,k)+ah(i-1,j,k)+ah(i,j,k+1)+ah(i-1,j,k+1))
 dif1 = ahe*(w(i,j,k+1)-w(i,j,k))/dx
 IF(dry(i,j,k+1)) dif1 = 0.0
 ahw = 0.25*(ah(i,j,k)+ah(i-1,j,k)+ah(i,j,k-1)+ah(i-1,j,k-1))
 dif2 = ahw*(w(i,j,k)-w(i,j,k-1))/dx
 IF(dry(i,j,k-1)) dif2 = 0.0
 ahn = 0.25*(ah(i,j,k)+ah(i-1,j,k)+ah(i,j+1,k)+ah(i-1,j+1,k))
 dif3 = ahn*(w(i,j+1,k)-w(i,j,k))/dy
 IF(dry(i,j+1,k)) dif3 = 0.0
 ahs = 0.25*(ah(i,j,k)+ah(i-1,j,k)+ah(i,j-1,k)+ah(i-1,j-1,k))
 dif4 = ahs*(w(i,j,k)-w(i,j-1,k))/dy
 IF(dry(i,j-1,k)) dif4 = 0.0
 difh = (dif1-dif2)/dx + (dif3-dif4)/dy

 azt = az(i-1,j,k)
 dif1 = azt*(w(i-1,j,k)-w(i,j,k))/dz
 IF(dry(i-1,j,k)) dif1 = 0.0
 azb = az(i,j,k)
 dif2 = azb*(w(i,j,k)-w(i+1,j,k))/dz
 IF(dry(i+1,j,k)) dif2 = 0.0 
 difz = (dif1-dif2)/dz
 diffw = dti*(difh+difz)  

  IF(wet(i,j,k))THEN
    um = 0.25*(u(i,j,k)+u(i,j,k-1)+u(i,j+1,k)+u(i,j+1,k-1))
    vm = 0.25*(v(i,j,k)+v(i,j,k+1)+v(i,j-1,k)+v(i,j-1,k+1))
    pressx = -drdx*(q(i,j,k+1)-q(i,j,k))-drdx*(p(i,j,k+1)-p(i,j,k))
    IF(wet(i,j,k+1)) du(i,j,k) = cos(alpha)*u(i,j,k) + sin(alpha)*vm + dti*pressx + advx(i,j,k) + diffu - u(i,j,k)
    pressy = -drdy*(q(i,j+1,k)-q(i,j,k))-drdy*(p(i,j+1,k)-p(i,j,k))
    IF(wet(i,j+1,k)) dv(i,j,k) = cos(alpha)*v(i,j,k) - sin(alpha)*um + dti*pressy + advy(i,j,k) + diffv - v(i,j,k)
    pressz = -drdz*(q(i-1,j,k)-q(i,j,k))
    IF(wet(i-1,j,k)) dw(i,j,k) = dti*pressz + advz(i,j,k) + diffw
  END IF

END DO
END DO
END DO

DO i = 1,nz
DO j = 1,ny
DO k = 1,nx
  ustar(i,j,k) = u(i,j,k)
  vstar(i,j,k) = v(i,j,k)
  wstar(i,j,k) = w(i,j,k)
END DO
END DO
END DO

!++++++++++++++++++++++++++++++
! start of short time stepping
!++++++++++++++++++++++++++++++

DO nrep = 1,nfrac

DO i = 1,nz
DO j = 1,ny
DO k = 1,nx
  ustar(i,j,k) = ustar(i,j,k)+du(i,j,k)/REAL(nfrac)
  vstar(i,j,k) = vstar(i,j,k)+dv(i,j,k)/REAL(nfrac)
  wstar(i,j,k) = wstar(i,j,k)+dw(i,j,k)/REAL(nfrac)
END DO
END DO
END DO

! boundary conditions
DO i = 1,nz
DO j = 0,ny+1
 ustar(i,j,0) = ustar(i,j,1) 
 ustar(i,j,nx+1) = ustar(i,j,nx) 
 vstar(i,j,0) = vstar(i,j,1)
 vstar(i,j,nx+1) = vstar(i,j,nx) 
 wstar(i,j,0) = wstar(i,j,1)
 wstar(i,j,nx+1) = wstar(i,j,nx) 
END DO
END DO

DO i = 1,nz
DO k = 0,nx+1
 ustar(i,0,k) = ustar(i,1,k)
 ustar(i,ny+1,k) = ustar(i,ny,k)
 vstar(i,0,k) = vstar(i,1,k) 
 vstar(i,ny+1,k) = vstar(i,ny,k) 
 wstar(i,0,k) = wstar(i,1,k)
 wstar(i,ny+1,k) = wstar(i,ny,k) 
END DO
END DO

! calculate right-hand side of Poisson equation
DO i = 1,nz
DO j = 1,ny
DO k = 1,nx
 qstar(i,j,k) = -1.*RHOREF/dte*(  &
 &  (ustar(i,j,k)-ustar(i,j,k-1))*dz + &
 &  (vstar(i,j,k)-vstar(i,j-1,k))*dx*dz/dy + &
 &  (wstar(i,j,k)-wstar(i+1,j,k))*dx )
END DO
END DO
END DO

! STEP 4: S.O.R. ITERATION

DO i = 0,nz+1
DO j = 0,ny+1
DO k = 0,nx+1
 dq(i,j,k) = 0.0
END DO
END DO
END DO

nstop = 8000

!*****************
DO nsor = 1,nstop
!*****************

perr = 0.0

! STEP 1: predict pressure correction
DO i = 1,nz
DO j = 1,ny
DO k = 1,nx

 IF(wet(i,j,k))THEN
 
 q1 = dq(i,j,k)
 term1 = qstar(i,j,k) + & 
  &      at(i,j,k)*dq(i-1,j,k) + ab(i,j,k)*dq(i+1,j,k) + & 
  &      aw(i,j,k)*dq(i,j,k-1) + ae(i,j,k)*dq(i,j,k+1) + &
  &      as(i,j,k)*dq(i,j-1,k) + an(i,j,k)*dq(i,j+1,k)
 q2 = (1.0-omega)*q1 + omega*term1/atot(i,j,k) 
 dq(i,j,k) = q2
 perr = MAX(ABS(q2-q1),perr)

 END IF

END DO
END DO
END DO

DO i = 1,nz
DO j = 1,ny
 dq(i,j,nx+1) = dq(i,j,nx)
END DO
END DO

DO i = 1,nz
DO k = 1,nx
 dq(i,ny+1,k) = dq(i,ny,k)
END DO
END DO

! STEP 2: predict new velocities 
DO i = 1,nz
DO j = 1,ny
DO k = 1,nx
  IF(wet(i,j,k))THEN
    pressx = -drdx*(dq(i,j,k+1)-dq(i,j,k))
    IF(wet(i,j,k+1))un(i,j,k) = ustar(i,j,k) + dte*pressx
    pressy = -drdy*(dq(i,j+1,k)-dq(i,j,k))
    IF(wet(i,j+1,k)) vn(i,j,k) = vstar(i,j,k) + dte*pressy
    pressz = -drdz*(dq(i-1,j,k)-dq(i,j,k))
    IF(wet(i-1,j,k)) wn(i,j,k) = wstar(i,j,k) + dte*pressz
  END IF
END DO
END DO
END DO

! STEP 3a: predict depth-integrated flow
DO j = 1,ny
DO k = 1,nx
  usum(j,k) = 0.
  vsum(j,k) = 0.
  DO i = 1,nz
    usum(j,k) = usum(j,k) + dz*un(i,j,k)
    vsum(j,k) = vsum(j,k) + dz*vn(i,j,k)
  END DO
END DO
END DO

! lateral boundary conditions
DO j = 1,ny 
 usum(j,0) = usum(j,1)
 usum(j,nx+1) = usum(j,nx)
END DO

DO k = 1,nx 
 vsum(0,k) = vsum(1,k)
 vsum(ny+1,k) = vsum(ny,k)
END DO

! STEP 3b: predict surface pressure field
DO j = 1,ny
DO k = 1,nx
 dq(0,j,k) = -dte*RHOREF*G*( (usum(j,k)-usum(j,k-1) )/dx + ( vsum(j,k)-vsum(j-1,k) )/dy )
END DO
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

DO i = 1,nz
DO j = 1,ny
DO k = 1,nx
  ustar(i,j,k) = un(i,j,k)
  vstar(i,j,k) = vn(i,j,k)
  wstar(i,j,k) = wn(i,j,k)
END DO
END DO
END DO


END DO
!++++++++++++++++++++++++++++
! end of short time stepping
!++++++++++++++++++++++++++++

! updating for next time step
DO i = 1,nz
DO j = 1,ny
DO k = 1,nx
  q(i,j,k) = q(i,j,k)+dq(i,j,k)
  u(i,j,k) = un(i,j,k)
  v(i,j,k) = vn(i,j,k)
  w(i,j,k) = wn(i,j,k)
END DO
END DO
END DO

DO j = 1,ny
DO k = 1,nx
  q(0,j,k) = q(0,j,k)+dq(0,j,k)
END DO
END DO

! lateral boundary conditions
DO j = 1,ny
 q(0,j,0) = q(0,j,1)
 q(0,j,nx+1) = q(0,j,nx)
 DO i = 1,nz
   u(i,j,0) = u(i,j,1) 
   u(i,j,nx+1) = u(i,j,nx) 
   v(i,j,0) = v(i,j,1)
   v(i,j,nx+1) = v(i,j,nx) 
   w(i,j,0) = w(i,j,1)
   w(i,j,nx+1) = w(i,j,nx)
   q(i,j,0) = q(i,j,1)
   q(i,j,nx+1) = q(i,j,nx)
 END DO
END DO

DO k = 1,nx
 q(0,1,k) = 0.0 ! keep sea level at zero along deep boundary
 q(0,0,k) = q(0,1,k)
 q(0,ny+1,k) = q(0,ny,k)
 DO i = 1,nz
   u(i,0,k) = u(i,1,k)
   u(i,ny+1,k) = u(i,ny,k) 
   v(i,0,k) = v(i,1,k)
   v(i,ny+1,k) = v(i,ny,k) 
   w(i,0,k) = w(i,1,k)
   w(i,ny+1,k) = w(i,ny,k)
   q(i,0,k) = q(i,1,k)
   q(i,ny+1,k) = q(i,ny,k)
 END DO
END DO

RETURN
END SUBROUTINE dyn

SUBROUTINE advect
! local parameters
REAL :: RxP(0:nz+1,0:ny+1,0:nx+1), RxN(0:nz+1,0:ny+1,0:nx+1)
REAL :: RyP(0:nz+1,0:ny+1,0:nx+1), RyN(0:nz+1,0:ny+1,0:nx+1)
REAL :: RzP(0:nz+1,0:ny+1,0:nx+1), RzN(0:nz+1,0:ny+1,0:nx+1)
REAL :: dB, term1, term2, term3, term4, term5, term6
REAL :: BwP, BwN, BeP, BeN
REAL :: BnP, BnN, BsP, BsN
REAL :: BbP, BbN, BtP, BtN 

DO i = 0,nz+1
DO j = 0,ny+1
DO k = 0,nx+1
  RxP(i,j,k) = 0.0
  RxN(i,j,k) = 0.0
  RyP(i,j,k) = 0.0
  RyN(i,j,k) = 0.0
  RzP(i,j,k) = 0.0
  RzN(i,j,k) = 0.0
END DO
END DO
END DO

DO i = 1,nz
DO j = 1,ny
DO k = 1,nx
  dB =  B(i,j,k+1)-B(i,j,k)
  IF(ABS(dB) > 0.0) RxP(i,j,k) = (B(i,j,k)-B(i,j,k-1))/dB
  dB =  B(i,j+1,k)-B(i,j,k)
  IF(ABS(dB) > 0.0) RyP(i,j,k) = (B(i,j,k)-B(i,j-1,k))/dB
  dB =  B(i-1,j,k)-B(i,j,k)
  IF(ABS(dB) > 0.0) RzP(i,j,k) = (B(i,j,k)-B(i+1,j,k))/dB
END DO
END DO
END DO

DO i = 1,nz
DO j = 1,ny
DO k = 0,nx-1
  dB =  B(i,j,k+1)-B(i,j,k)
  IF(ABS(dB) > 0.0) RxN(i,j,k) = (B(i,j,k+2)-B(i,j,k+1))/dB
END DO
END DO
END DO

DO i = 1,nz
DO j = 0,ny-1
DO k = 1,nx
  dB =  B(i,j+1,k)-B(i,j,k)
  IF(ABS(dB) > 0.0) RyN(i,j,k) = (B(i,j+2,k)-B(i,j+1,k))/dB
END DO
END DO
END DO

DO i = 2,nz+1
DO j = 1,ny
DO k = 1,nx
  dB =  B(i-1,j,k)-B(i,j,k)
  IF(ABS(dB) > 0.0) RzN(i,j,k) = (B(i-2,j,k)-B(i-1,j,k))/dB
END DO
END DO
END DO   

DO i = 1,nz
DO j = 1,ny
DO k = 1,nx

term1 = (1.0-CuP(i,j,k-1))*(B(i,j,k)-B(i,j,k-1))
BwP = B(i,j,k-1)+0.5*PSI(RxP(i,j,k-1))*term1

term1 = (1.0+CuN(i,j,k-1))*(B(i,j,k)-B(i,j,k-1))
BwN = B(i,j,k)-0.5*PSI(RxN(i,j,k-1))*term1

term1 = (1.0-CuP(i,j,k))*(B(i,j,k+1)-B(i,j,k))
BeP = B(i,j,k)+0.5*PSI(RxP(i,j,k))*term1

term1 = (1.0+CuN(i,j,k))*(B(i,j,k+1)-B(i,j,k))  
BeN = B(i,j,k+1)-0.5*PSI(RxN(i,j,k))*term1

term1 = (1.0-CvP(i,j-1,k))*(B(i,j,k)-B(i,j-1,k))
BsP = B(i,j-1,k)+0.5*PSI(RyP(i,j-1,k))*term1

term1 = (1.0+CvN(i,j-1,k))*(B(i,j,k)-B(i,j-1,k))
BsN = B(i,j,k)-0.5*PSI(RyN(i,j-1,k))*term1

term1 = (1.0-CvP(i,j,k))*(B(i,j+1,k)-B(i,j,k))
BnP = B(i,j,k)+0.5*PSI(RyP(i,j,k))*term1

term1 = (1.0+CvN(i,j,k))*(B(i,j+1,k)-B(i,j,k))  
BnN = B(i,j+1,k)-0.5*PSI(RyN(i,j,k))*term1

term1 = (1.0-CwP(i+1,j,k))*(B(i,j,k)-B(i+1,j,k))
BbP = B(i+1,j,k)+0.5*PSI(RzP(i+1,j,k))*term1

term1 = (1.0+CwN(i+1,j,k))*(B(i,j,k)-B(i+1,j,k))  
BbN = B(i,j,k)-0.5*PSI(RzN(i+1,j,k))*term1

term1 = (1.0-CwP(i,j,k))*(B(i-1,j,k)-B(i,j,k)) 
BtP = B(i,j,k)+0.5*PSI(RzP(i,j,k))*term1

term1 = (1.0+CwN(i,j,k))*(B(i-1,j,k)-B(i,j,k)) 
BtN = B(i-1,j,k)-0.5*PSI(RzN(i,j,k))*term1


term1 = CuP(i,j,k-1)*BwP+CuN(i,j,k-1)*BwN
term2 = CuP(i,j,k)*BeP+CuN(i,j,k)*BeN

term3 = CvP(i,j-1,k)*BsP+CvN(i,j-1,k)*BsN
term4 = CvP(i,j,k)*BnP+CvN(i,j,k)*BnN

term5 = CwP(i+1,j,k)*BbP+CwN(i+1,j,k)*BbN
term6 = CwP(i,j,k)*BtP+CwN(i,j,k)*BtN

BN(i,j,k) = term1-term2+term3-term4+term5-term6

END DO
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
REAL, PARAMETER :: vismin = 1.e-4
REAL, PARAMETER :: vismax = 0.05
REAL :: windmix

REAL :: utop, ubot, dudz
REAL :: vtop, vbot, dvdz
REAL :: eddyx, eddyy, stab, wurz
REAL :: c2

windmix = vismax

c2 = 0.15
c2 = c2*c2*dz*dz

DO k = 1,nx
DO j = 1,ny
DO i = 1,nz

utop = 0.5*(u(i-1,j,k)+u(i-1,j,k-1))
ubot = 0.5*(u(i+1,j,k)+u(i+1,j,k-1))
dudz = (utop-ubot)/(2.0*dz)
eddyx = dudz*dudz

vtop = 0.5*(v(i-1,j,k)+v(i-1,j-1,k))
vbot = 0.5*(v(i+1,j,k)+v(i+1,j-1,k))
dvdz = (vtop-vbot)/(2.0*dz)
eddyy = dvdz*dvdz

stab = -G/RHOREF*(rho(i-1,j,k)-rho(i+1,j,k))/(2.0*dz)
wurz = AMAX1(0.0, eddyx+eddyy-stab)
wurz = AMAX1(c2*SQRT(wurz),vismin)
wurz = AMIN1(wurz,vismax)
IF(stab < 0)wurz = vismax
IF(i<2) wurz = windmix
kz(i,j,k) = wurz

END DO
END DO
END DO

DO j = 1,ny
DO k = 1,nx
  kz(0,j,k) = kz(1,j,k)
  kz(nz+1,j,k) = kz(nz,j,k)
END DO
END DO

DO i = 0,nz+1
DO k = 1,nx
  kz(i,0,k) = kz(i,1,k)
  kz(i,ny+1,k) = kz(i,ny,k)
END DO
END DO

DO i = 0,nz+1
DO j = 1,ny
  kz(i,j,0) = kz(i,j,1)
  kz(i,j,nx+1) = kz(i,j,nx)
END DO
END DO

DO i = 0,nz+1
DO j = 0,ny+1
DO k = 0,nx+1
  az(i,j,k) = kz(i,j,k)
END DO
END DO
END DO

END SUBROUTINE eddy


END MODULE sub