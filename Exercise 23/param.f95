MODULE param

INTEGER(4), PARAMETER :: nx = 100 ! number of grid points x-direction
INTEGER(4), PARAMETER :: ny = 50 ! number of grid points y-direction
INTEGER(4), PARAMETER :: nz = 11 ! number of grid points z-direction
REAL, PARAMETER :: G = 9.81 ! acceleration due to gravity
REAL, PARAMETER :: RHOREF = 1028.0 ! reference density
REAL, PARAMETER :: PI = 3.14159265359 ! pi

REAL :: dt, dx, dy, dz ! grid parameters
REAL :: drdz, drdy, drdx ! derived parameters
REAL :: ad ! adjustment parameter
REAL :: kh(0:nz+1,0:ny+1,0:nx+1), kz(0:nz+1,0:ny+1,0:nx+1) ! density diffusivities
REAL :: ah(0:nz+1,0:ny+1,0:nx+1), az(0:nz+1,0:ny+1,0:nx+1) ! eddy viscosities
REAL :: r ! bottom friction coefficient
REAL :: depth(0:ny+1,0:nx+1)

INTEGER :: n,i,j,k ! time level and grid indices

INTEGER :: ib(0:ny+1,0:nx+1)  ! bottom cell pointers 
REAL :: tbu(ny,nx), tbv(ny,nx) ! bottom-stress parameters
REAL :: taux ! wind-stress in x-direction

! density arrays
REAL :: rho(0:nz+1,0:ny+1,0:nx+1),rhon(0:nz+1,0:ny+1,0:nx+1), rhozero(0:nz+1)
! pressure arrays
REAL :: p(0:nz+1,0:ny+1,0:nx+1) ! hydrostatic pressure
REAL :: dq(0:nz+1,0:ny+1,0:nx+1) ! pressure correction
REAL :: q(0:nz+1,0:ny+1,0:nx+1) ! nonhydrostatic pressure
! velocity arrays
REAL :: u(0:nz+1,0:ny+1,0:nx+1), un(0:nz+1,0:ny+1,0:nx+1) 
REAL :: v(0:nz+1,0:ny+1,0:nx+1), vn(0:nz+1,0:ny+1,0:nx+1) 
REAL :: w(0:nz+1,0:ny+1,0:nx+1), wn(0:nz+1,0:ny+1,0:nx+1)
! concentration array  
REAL :: c(0:nz+1,0:ny+1,0:nx+1), cn(0:nz+1,0:ny+1,0:nx+1)  

REAL :: usum(0:ny+1,0:nx+1), vsum(0:ny+1,0:nx+1) ! depth-integrated flow
LOGICAL :: wet(0:nz+1,0:ny+1,0:nx+1) ! wet pointer
LOGICAL :: dry(0:nz+1,0:ny+1,0:nx+1) ! dry pointer

! parameters for S.O.R iteration
REAL :: ustar(0:nz+1,0:ny+1,0:nx+1), vstar(0:nz+1,0:ny+1,0:nx+1)
REAL :: wstar(0:nz+1,0:ny+1,0:nx+1), qstar(nz,ny,nx)
REAL :: omega 
REAL :: at(nz,ny,nx), ab(nz,ny,nx), ae(nz,ny,nx), aw(nz,ny,nx)
REAL :: atot(nz,ny,nx), an(nz,ny,nx), as(nz,ny,nx) 
REAL :: peps 

! parameters for Coriolis force
REAL :: f, alpha

! arrays for advection subroutine
REAL :: CuP(0:nz+1,0:ny+1,0:nx+1), CuN(0:nz+1,0:ny+1,0:nx+1)
REAL :: CvP(0:nz+1,0:ny+1,0:nx+1), CvN(0:nz+1,0:ny+1,0:nx+1)
REAL :: CwP(0:nz+1,0:ny+1,0:nx+1), CwN(0:nz+1,0:ny+1,0:nx+1)
REAL :: B(0:nz+1,0:ny+1,0:nx+1), BN(0:nz+1,0:ny+1,0:nx+1)

END MODULE param