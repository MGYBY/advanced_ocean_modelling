MODULE param

INTEGER(4), PARAMETER :: nz = 500 ! number of vertical grid points

REAL :: u(0:nz+1), v(0:nz+1) ! horizontal velocity at time level n
REAL :: un(0:nz+1), vn(0:nz+1) ! horizontal velocity at time level n+1
REAL :: az(0:nz+1) ! vertical eddy viscosity
REAL :: alpha, beta ! parameters for semi-implicit treatment of the Coriolis force

REAL :: dt,dz ! time step and grid spacing
INTEGER :: i ! vertical cell index
INTEGER :: mode ! various choices for turbulence closure

REAL :: taux, tauy ! wind stress vector
REAL :: rho ! constant density
REAL :: f ! Coriolis parameter
REAL :: az0 ! uniform background viscosity

END MODULE param