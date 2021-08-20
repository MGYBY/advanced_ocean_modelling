MODULE sub
USE param

CONTAINS

!******************
SUBROUTINE init
!local parameters
REAL :: ccc, cmc

! set all arrays to zero
DO i = 0,nz+1
DO k = 0,nx+1
 dp(i,k) = 0.0
 u(i,k) = 0.0
 un(i,k) = 0.0
 w(i,k) = 0.0
 wn(i,k) = 0.0
END DO
END DO

DO k = 0,nx+1
 q(k) = 0.0
END DO

! grid parameters
dx = 5.0
dz = 2.0
dt = 0.05

! coefficients for SOR
omega = 1.4
peps = 1.e-2

DO i = 1,nz
DO k = 1,nx
  ct(i,k) = dx/dz
  cb(i,k) = dx/dz
  ce(i,k) = dz/dx
  cw(i,k) = dz/dx
  IF(k==1) cw(i,k) = 0.0
  IF(k==nx) ce(i,k) = 0.0
  IF(i==nz) cb(i,k) = 0.0
  ctot(i,k) = cb(i,k)+ct(i,k)+ce(i,k)+cw(i,k)
END DO
END DO

END SUBROUTINE init

!*******************
SUBROUTINE dyn

! local parameters and arrays
REAL :: ustar(0:nz+1,0:nx+1), wstar(0:nz+1,0:nx+1)
REAL :: pressx, drdxh, pressz, drdzh
REAL :: dpstore(nx+1), perr, p1, p2, term1
INTEGER :: nsor, nstop

drdxh = 0.5/(RHO*dx)
drdzh = 0.5/(RHO*dz)

! sea-level forcing
dp(0,2) = ad*RHO*G

! STEP 1: Store surface pressure field
DO k = 0,nx+1
 dpstore(k) = dp(0,k)
END DO

! STEP 2: Calculate ustar and vstar 
DO i = 1,nz
DO k = 1,nx
  pressx = -drdxh*(dp(i,k+1)-dp(i,k))
  ustar(i,k) = u(i,k) + dt*pressx
  pressz = -drdzh*(dp(i-1,k)-dp(i,k))
  wstar(i,k) = w(i,k) + dt*pressz
END DO
END DO

! lateral boundary conditions
DO i = 1,nz
 ustar(i,nx) = 0.0
 ustar(i,nx+1) = ustar(i,nx)
 ustar(i,0) = 0.0 
END DO

! bottom boundary condition
DO k = 1,nx
 wstar(nz+1,k) = 0.0
END DO

! STEP 3: calculate right-hand side of poisson equation
DO i = 1,nz
DO k = 1,nx
 pstar(i,k) = -2.0*RHO/dt*(  &
 &  (ustar(i,k)-u(i,k)-ustar(i,k-1)+u(i,k-1))*dz + &
 &  (wstar(i,k)-w(i,k)-wstar(i+1,k)+w(i+1,k))*dx )
END DO
END DO

! STEP 4: S.O.R. ITERATION

nstop = 5000

!*****************
DO nsor = 1,5000
!*****************

perr = 0.0

! STEP 4.1: predict new pressure
DO i = 1,nz
DO k = 1,nx
 p1 = dp(i,k)
 term1 = pstar(i,k) + & 
  &      ct(i,k)*dp(i-1,k) + cb(i,k)*dp(i+1,k) + & 
  &      cw(i,k)*dp(i,k-1) + ce(i,k)*dp(i,k+1)  
 p2 = (1.0-omega)*p1 + omega*term1/ctot(i,k) 
 dp(i,k) = p2
 perr = MAX(ABS(p2-p1),perr)
END DO
END DO

DO i = 1,nz
  dp(i,0) = dp(i,1)
  dp(i,nx+1) = dp(i,nx)
END DO

! STEP 4.2: predict new velocities 
DO i = 1,nz
DO k = 1,nx
  pressx = -drdxh*(dp(i,k+1)-dp(i,k))
  un(i,k) = ustar(i,k) + dt*pressx
  pressz = -drdzh*(dp(i-1,k)-dp(i,k))
  wn(i,k) = wstar(i,k) + dt*pressz
END DO
END DO

DO i = 1,nz
 un(i,nx) = 0.0
 un(i,nx+1) = un(i,nx)
 un(i,0) = 0.0 
END DO

! STEP 4.3: predict depth-integrated flow
DO k = 1,nx
  q(k) = 0.
  DO i = 1,nz
    q(k) = q(k) + dz*un(i,k)
  END DO
END DO

! lateral boundary conditions
q(0) = 0.0
q(nx) = 0.0
q(nx+1) = q(nx)

! STEP 4.4: predict surface pressure field
DO k = 1,nx
 dp(0,k) = dpstore(k)-dt*RHO*G*(q(k)-q(k-1))/dx
END DO

IF(perr <= peps)THEN
  nstop = nsor
  GOTO 33
END IF

!********************
END DO
!********************

 32 GOTO 34
 33 WRITE(*,*) "No. of Interactions =>", nstop
 34 CONTINUE

! updating for next time step
DO i = 1,nz
DO k = 1,nx
  u(i,k) = un(i,k)
  w(i,k) = wn(i,k)
END DO
END DO

! lateral boundary conditions
DO i = 1,nz
 u(i,0) = 0.0
 u(i,nx) = 0.0
 u(i,nx+1) = 0.0
END DO

! vertical boundary conditions
DO k = 1,nx
 w(nz+1,k) = 0.0
END DO

RETURN
END SUBROUTINE dyn

END MODULE sub