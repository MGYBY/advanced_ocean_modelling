MODULE sub
USE param

CONTAINS

!=======================
SUBROUTINE init

! set parameters
dz = 1.0 ! vertical grid spacing is 1 metre
dt = 5.0 ! time step is 5 seconds
rho = 1028.0 ! typical seawater density
taux = 0.0 ! no east-west wind stress
tauy = 0.5 ! southerly wind stress is 0.5 Pa
f = 1.e-4 ! Coriolis parameter

alpha = dt*f
beta = 0.25*alpha*alpha

! initial conditions
ugeo = 0.0
vgeo = 0.1

DO i = 0,nz+1
  u(i) = 0.0
  v(i) = 0.0
END DO

END SUBROUTINE init

!=======================
SUBROUTINE eddy
! local parameters
REAL :: az0, azmix, azmin
REAL :: zlength, dudz, dvdz, shear

az0 = 5.e-3 ! uniform eddy viscosity
winmix = 1.e-1 ! value near sea surface
azmin = 4.e-3 ! local minimum value of eddy viscosity
zlength = 2.0 ! Prandtl mixing length

! settings for Mode = 1
DO i = 0,nz+1
  az(i) = az0
END DO

END SUBROUTINE eddy

!=======================
SUBROUTINE dyn
! local parameters
REAL :: diffu, diffv, a, b, c
REAL :: atop, abot ! needed for eddy viscosity averaging

a = sin(f*dt)
b = cos(f*dt)

! surface boundary conditions

u(0) = 0.0
v(0) = 0.0
u(nz+1) = -ugeo
v(nz+1) = -vgeo

DO i = 1,nz
  atop = 0.5*(az(i-1)+az(i))
  abot = 0.5*(az(i)+az(i+1))
  diffu = dt*(atop*(u(i-1)-u(i))/dz-abot*(u(i)-u(i+1))/dz)/dz
  diffv = dt*(atop*(v(i-1)-v(i))/dz-abot*(v(i)-v(i+1))/dz)/dz
  un(i) = b*u(i)+a*v(i)+diffu
  vn(i) = b*v(i)-a*u(i)+diffv
END DO

un(nz+1) = u(nz+1)
vn(nz+1) = v(nz+1)

! updating
DO i = 1,nz+1
   u(i) = un(i)
   v(i) = vn(i)
END DO

END SUBROUTINE dyn

END MODULE sub