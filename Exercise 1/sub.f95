MODULE sub
USE param

CONTAINS

!=======================
SUBROUTINE init

! set parameters
dz = 1.0 ! vertical grid spacing is 1 metre
dt = 5.0 ! time step is 5 seconds
rho = 1028.0 ! typical seawater density
taux = 0.0 ! east-west wind stress disables
tauy = 0.5 ! southerly wind stress is set to 0.5 Pa
f = 1.e-4 ! Coriolis parameter

alpha = dt*f
beta = 0.25*alpha*alpha

! initial conditions
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

az0 = 5.e-2 ! uniform eddy viscosity
winmix = 1.e-1 ! value near sea surface
azmin = 4.e-3 ! local minimum value of eddy viscosity
zlength = 2.0 ! Prandtl mixing length

! settings for Mode = 1
DO i = 0,nz+1
  az(i) = az0
END DO

! settings for Mode = 2
IF(mode == 2)THEN
 az(19)=azmin 
 az(20)=azmin
 az(21)=azmin
END IF

! settings for Mode = 3
IF(mode == 3)THEN
az(0) = winmix
DO i = 1,nz
  dudz = (u(i-1)-u(i+1))/(2.*dz)
  dvdz = (v(i-1)-v(i+1))/(2.*dz)
  shear = sqrt(dudz*dudz+dvdz*dvdz)
  az(i) = AMAX1(zlength*zlength*shear,azmin)
END DO
END IF

END SUBROUTINE eddy

!=======================
SUBROUTINE dyn
! local parameters
REAL :: diffu, diffv, a, b, c
REAL :: atop, abot ! needed for eddy viscosity averaging

a = alpha
b = 1.0-beta
c = 1.0+beta

! surface boundary conditions

atop = 0.5*(az(0)+az(1))
u(0) = u(1)+dz*taux/rho/atop
v(0) = v(1)+dz*tauy/rho/atop

DO i = 1,nz
  atop = 0.5*(az(i-1)+az(i))
  abot = 0.5*(az(i)+az(i+1))
  diffu = dt*(atop*(u(i-1)-u(i))/dz-abot*(u(i)-u(i+1))/dz)/dz
  diffv = dt*(atop*(v(i-1)-v(i))/dz-abot*(v(i)-v(i+1))/dz)/dz
  un(i) = (b*u(i)+a*v(i)+0.5*a*diffv+diffu)/c
  vn(i) = (b*v(i)-a*u(i)-0.5*a*diffu+diffv)/c
END DO

! bottom boundary conditions (no stress)
un(nz+1) = un(nz)
vn(nz+1) = vn(nz)

! updating
DO i = 1,nz+1
   u(i) = un(i)
   v(i) = vn(i)
END DO

END SUBROUTINE dyn

END MODULE sub