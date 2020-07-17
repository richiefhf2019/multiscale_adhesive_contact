module constitutive

use math
implicit none

real(8) :: e, v, bulk, viscosity, c1, c2

contains

subroutine Calculate_PK1(matid,nmat,vmat,f,df, pk1)

implicit none

integer, intent( in) :: matid, nmat
real(8), intent( in) :: vmat(nmat), f(3,3), df(3,3)
real(8), intent(out) :: pk1(3,3)

select case(matid)
  case(1)
    bulk = vmat(1); viscosity = vmat(2)
    call Newtonian_Fluid(bulk,viscosity,f,df, pk1)
  case(2)
    e = vmat(1); v = vmat(2)
    call Linear_Elastic3d(e,v,f, pk1)
  case default
    write(*,*) 'Wrong material id in Cacluate_PK1, Quit !'
    stop
end select


return
end subroutine Calculate_PK1

subroutine Linear_Elastic2d(e, v, f, pk1)

implicit none

real(8), intent( in) :: e, v,  f(3,3)
real(8), intent(out) :: pk1(3,3)
! local variables
real(8) :: sigma(3), sts_cauchy(3,3), eps(3)
real(8) :: aa(3,3), factor, finv(3,3), ff(3,3), dudx(3,3)
integer :: i, j, k

! Tangent Modulus
aa = 0.d0; aa(1,1) = 1.d0; aa(2,2) = 1.d0
aa(3,3) = (1.d0-v)/2.d0; aa(2,1) = v; aa(1,2) = v
! For Plane Stress
factor = e/(1-v**2); aa = factor*aa 

! Get the Engineering Strain
ff = f
do i = 1, 3
  ff(i,i) = ff(i,i) - 1.d0
end do
finv = inverse(f)
dudx = matrixmul(ff,finv)

! Arrange the Engineering Strain in Voight Notation
eps(1) = dudx(1,1); eps(2) = dudx(2,2); eps(3) = dudx(1,2) + dudx(2,1)
! Stress in Voight Notation
sigma = matmul(aa,eps)

! Cauchy Stress
sts_cauchy = 0.d0
sts_cauchy(1,1) = sigma(1); sts_cauchy(2,2) = sigma(2)
sts_cauchy(1,2) = sigma(3); sts_cauchy(2,1) = sigma(3)

! PK1 Stress
sts_cauchy = sts_cauchy*determinant(f)
ff = transpose(f)
finv = inverse(ff)
pk1 = matrixmul(sts_cauchy,finv)

return
end subroutine Linear_Elastic2d



subroutine Linear_Elastic3d(e, v, f, pk1)

!------------------------------------------------------------------------------
!
! e-- Young's Modulus, v-- Possion's Ratio
! f-- deformation gradient
! pk1-- first Piola-Kirchhoff stress
!
!------------------------------------------------------------------------------

implicit none

real(8), intent( in) :: e, v, f(3,3)
real(8), intent(out) :: pk1(3,3)
! local variables
real(8) :: sigma(6), sts_cauchy(3,3), eps(6), lambda, mu
real(8) :: factor, finv(3,3), ff(3,3), dudx(3,3)
integer :: i, j, k

! Lame Constants
lambda = v*E/((1.d0+v)*(1-2.d0*v)); mu = E/2.d0/(1.d0+v)

! Get the Engineering Strain
ff = f
do i = 1, 3
  ff(i,i) = ff(i,i) - 1.d0
end do
finv = inverse(f)
dudx = matrixmul(ff,finv)

! Arrange the Engineering Strain in Voight Notation
eps(1) = dudx(1,1); eps(2) = dudx(2,2); eps(3) = dudx(3,3)
eps(4) = dudx(1,2) + dudx(2,1); eps(5) = dudx(1,3) + dudx(3,1)
eps(6) = dudx(2,3) + dudx(3,2)

! 

! Stress in Voight Notation
sigma(1) = lambda*(eps(1)+eps(2)+eps(3)) + 2.d0*mu*eps(1)
sigma(2) = lambda*(eps(1)+eps(2)+eps(3)) + 2.d0*mu*eps(2)
sigma(3) = lambda*(eps(1)+eps(2)+eps(3)) + 2.d0*mu*eps(3)
sigma(4) = mu*eps(4)
sigma(5) = mu*eps(5)
sigma(6) = mu*eps(6)

! Cauchy Stress
sts_cauchy(1,1) = sigma(1); sts_cauchy(1,2) = sigma(4); sts_cauchy(1,3) = sigma(5)
sts_cauchy(2,1) = sigma(4); sts_cauchy(2,2) = sigma(2); sts_cauchy(2,3) = sigma(6)
sts_cauchy(3,1) = sigma(5); sts_cauchy(3,2) = sigma(6); sts_cauchy(3,3) = sigma(3)

! PK1 Stress
sts_cauchy = sts_cauchy*determinant(f)
ff = transpose(f)
finv = inverse(ff)
pk1 = matrixmul(sts_cauchy,finv)

return
end subroutine Linear_Elastic3d


subroutine Newtonian_Fluid(bulk, viscosity, f, df, pk1)

implicit none

real(8), intent( in) :: bulk, viscosity, f(3,3), df(3,3)
real(8), intent(out) :: pk1(3,3)
! local variables
integer :: i
real(8) :: finv(3,3), lij(3,3), sts_cauchy(3,3), det

det = determinant(f)
finv = inverse(f)
lij = matrixmul(df,finv) ! L = \dot{F} F^{-1}

! Calculate the Cauchy Stress
sts_cauchy = 0.d0
do i = 1, 3
  sts_cauchy(i,i) = -(1.0-det)*bulk
end do
sts_cauchy = sts_cauchy + viscosity*(lij + transpose(lij))

! PK1 = J \sigma F^{-T}
pk1 = det*matrixmul(sts_cauchy,transpose(finv))

return
end subroutine Newtonian_Fluid

subroutine Mooney_Rivlin(c1, c2, f, pk1)
!------------------------------------------------
! Hyperelastic Mooney-Rivlin Model
! W = c1(I1 - 3(I3)^(1/3)) + c2(I2 - 3(I3)^(1/3))
! input--c1, c2, F
! output--pk1-->1st PK stress
!------------------------------------------------

implicit none

real(8), intent( in) :: c1, c2, f(3,3)
real(8), intent(out) :: pk1(3,3)

integer :: i
real(8) :: C(3,3), det, trace


pk1 = 0.d0

! Right-Cauchy Green Tensor
C = matmul(transpose(f),f)
det = determinant(C); trace = C(1,1) + C(2,2) + C(3,3)


pk1 = 2.0 * add(-c2*C, -(c1*det**(1.0/3.0) + 2.0 * c2 * det**(2.0/3.0)) * inverse(C))

do i = 1,3
  pk1(i,i) = pk1(i,i) + 2.0 * (c1 + c2 * trace)
end do

return
end subroutine Mooney_Rivlin

!subroutine Neo_Hookean(lambda,mu,defgrad, sigma,aa)
!!------------------------------------------------------------------------------
!! Hyperelastic Neo_Hookean Model
!! W = \frac{\Lambda}{2}(\ln J)^2 + \frac{mu}{2}(I_1-3) - \mu \ln J
!! input: lambda, mu, defgrad
!! output: Cauchy stress
!!------------------------------------------------------------------------------
!
!implicit none
!
!real(8), intent( in) :: lambda, mu, defgrad(3,3)
!real(8), intent(out) :: sigma(6), aa(6,6)
!
!! page 25 in Roger's thesis
!
!return
!end subroutine Neo_Hookean


end module constitutive