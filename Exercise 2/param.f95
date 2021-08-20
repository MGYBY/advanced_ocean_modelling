MODULE param

INTEGER(4), PARAMETER :: nz = 500	

REAL :: u(0:nz+1), v(0:nz+1) ! actual velocity
REAL :: un(0:nz+1), vn(0:nz+1) ! velocity at time level n+1
REAL :: az(0:nz+1) ! eddy viscosity
REAL :: alpha, beta
REAL :: ugeo, vgeo


REAL :: dt,dz ! time step and grid spacing
INTEGER :: i ! cell index
INTEGER :: mode ! various choices for turbulence closure

REAL :: taux, tauy ! wind stress vector
REAL :: rho ! constant density
REAL :: f ! Coriolis parameter
REAL :: az0 ! uniform background viscosity

REAL :: rx, ry ! needed for semi-implicit Coriolis force

END MODULE param