MODULE param

INTEGER(4), PARAMETER :: nx = 101 ! horizontal dimension	
INTEGER(4), PARAMETER :: nz = 51 ! vertical dimension
REAL, PARAMETER :: G = 9.81 ! acceleration due to gravity
REAL, PARAMETER :: RHO = 1028.0 ! reference density
REAL, PARAMETER :: PI = 3.1416 ! pi

REAL :: dt, dx, dz, ad
INTEGER :: i,k

REAL :: dp(0:nz+1,0:nx+1) ! dynamic pressure
REAL :: u(0:nz+1,0:nx+1), un(0:nz+1,0:nx+1) ! horizontal velocity 
REAL :: w(0:nz+1,0:nx+1), wn(0:nz+1,0:nx+1) ! vertical velocity 
REAL :: q(0:nx+1) ! depth-integrated flow

REAL :: period, amplitude ! forcing parameters

! matrices and parameters for S.O.R iteration
REAL :: omega 
REAL :: pstar(nz,nx)
REAL :: ct(nz,nx), cb(nz,nx) 
REAL :: ce(nz,nx), cw(nz,nx)
REAL :: ctot(nz,nx) 
REAL :: peps ! pressure accuracy

END MODULE param