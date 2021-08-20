MODULE sub
USE param

CONTAINS

!******************
SUBROUTINE init
!local parameters
REAL :: ccc, cmc, depth(nx)
INTEGER :: nb

! set all arrays to zero
DO i = 0,nz+1
DO k = 0,nx+1
 dp(i,k) = 0.0
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
 q(k) = 0.0
END DO

! grid parameters
dx = 5.0
dz = 2.0
dt = 0.05

DO i = 0,nz+1
 wet(i,0) = .false.
 wet(i,nx+1) = .false.
END DO

! variable bottom topography
DO k = 1,51
 depth(k) =  100.0-95.*REAL(k)/REAL(51)
END DO
DO k = 52,nx
 depth(k) = 5.0+95.0*REAL(k-51)/REAL(nx-51)
END DO

OPEN(50,file ='h.dat',form='formatted')
  WRITE(50,'(101F12.6)')(depth(k),k=1,nx)
CLOSE(50)

! wet/dry pointer
DO k = 1,nx
  nb = INT(depth(k)/dz)
  nbot = MIN(nbot,nz+1)
  DO i = nb,nz+1
    wet(i,k) = .false.
  END DO
END DO

! coefficients for SOR
omega = 1.4
peps = 1.e-2

DO i = 1,nz
DO k = 1,nx
  ct(i,k) = dx/dz
  cb(i,k) = dx/dz
  ce(i,k) = dz/dx
  cw(i,k) = dz/dx
  IF(.not.wet(i,k-1)) cw(i,k) = 0.0
  IF(.not.wet(i,k+1)) ce(i,k) = 0.0
  IF(.not.wet(i+1,k)) cb(i,k) = 0.0
  IF(.not.wet(i-1,k)) ct(i,k) = 0.0
  ctot(i,k) = cb(i,k)+ct(i,k)+ce(i,k)+cw(i,k)
END DO
END DO

END SUBROUTINE init

!*******************
SUBROUTINE dyn

! local parameters
REAL :: pressx, drdxh, pressz, drdzh
REAL :: dpstore(0:nx+1), perr, p1, p2, term1
INTEGER :: nsor, nstop

drdxh = 0.5/(RHO*dx)
drdzh = 0.5/(RHO*dz)

! sea-level forcing
dp(0,1) = ad*RHO*G

! STEP 1: Store surface pressure field
DO k = 0,nx+1
 dpstore(k) = dp(0,k)
END DO

! STEP 2: Calculate ustar and vstar 
DO i = 1,nz
DO k = 1,nx
  IF(wet(i,k))THEN
  pressx = -drdxh*(dp(i,k+1)-dp(i,k))
  IF(wet(i,k+1)) ustar(i,k) = u(i,k) + dt*pressx
  pressz = -drdzh*(dp(i-1,k)-dp(i,k))
  IF(wet(i-1,k)) wstar(i,k) = w(i,k) + dt*pressz
  END IF
END DO
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

 IF(wet(i,k))THEN
 
 p1 = dp(i,k)
 term1 = pstar(i,k) + & 
  &      ct(i,k)*dp(i-1,k) + cb(i,k)*dp(i+1,k) + & 
  &      cw(i,k)*dp(i,k-1) + ce(i,k)*dp(i,k+1)  
 p2 = (1.0-omega)*p1 + omega*term1/ctot(i,k) 
 dp(i,k) = p2
 perr = MAX(ABS(p2-p1),perr)

 END IF

END DO
END DO

DO i = 1,nz
  dp(i,0) = dp(i,1)
  dp(i,nx+1) = dp(i,nx)
END DO

! STEP 4.2: predict new velocities 
DO i = 1,nz
DO k = 1,nx
  IF(wet(i,k))THEN

  pressx = -drdxh*(dp(i,k+1)-dp(i,k))
  IF(wet(i,k+1)) un(i,k) = ustar(i,k) + dt*pressx
  pressz = -drdzh*(dp(i-1,k)-dp(i,k))
  IF(wet(i-1,k)) wn(i,k) = wstar(i,k) + dt*pressz

  END IF
END DO
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

RETURN
END SUBROUTINE dyn

END MODULE sub