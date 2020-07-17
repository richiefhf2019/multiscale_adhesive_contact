module fem

implicit none

contains

subroutine int2d(l, lint,sg)
!
!--------------------------------------------------------------------
! Purpose: Form Gauss points and weights for two dimensions
! 
! Inputs:
!    l          - Number of points/direction
! 
! Outputs:
!    lint       - Total number of points
!    sg(3,lint) - Array of points and weights
!--------------------------------------------------------------------
!
use config, only : nen
implicit  none

integer  :: l, lint

integer  ::  i,j,k,lr(9),lz(9),lw(9)
real(8)  ::  g,h, third, sg(3,*),ss(5),ww(5)

data      lr/-1,1,1,-1,0,1,0,-1,0/,lz/-1,-1,1,1,-1,0,1,0,0/
data      lw/4*25,4*40,64/
data      third / 0.3333333333333333d0 /

! Set number of total points

lint = l*l

! 5 pt. integration

if(l.eq.0) then

  lint = 5
  g    = sqrt(0.6d0)
  do i = 1,4
    sg(1,i) = g*lr(i)
    sg(2,i) = g*lz(i)
    sg(3,i) = 5.d0/9.d0
  end do ! i

  sg(1,5) = 0.0d0
  sg(2,5) = 0.0d0
  sg(3,5) = 16.d0/9.d0

! 1x1 integration

elseif(l.eq.1) then
  sg(1,1) = 0.d0
  sg(2,1) = 0.d0
  if(nen.eq.3) sg(2,1) = -third
  sg(3,1) = 4.d0

! 2x2 integration

elseif(l.eq.2) then
  g = sqrt(third)
  do i = 1,4
    sg(1,i) = g*lr(i)
    sg(2,i) = g*lz(i)
    sg(3,i) = 1.d0
  end do ! i

! 3x3 integration

elseif(l.eq.3) then
  g = sqrt(0.6d0)
  h = 1.d0/81.d0
  do i = 1,9
    sg(1,i) = g*lr(i)
    sg(2,i) = g*lz(i)
    sg(3,i) = h*lw(i)
  end do ! i

! 4x4 integration

elseif(l.eq.4) then
  g     = sqrt(4.8d0)
  h     = third/g
  ss(1) = sqrt((3.d0+g)/7.d0)
  ss(4) = - ss(1)
  ss(2) = sqrt((3.d0-g)/7.d0)
  ss(3) = -ss(2)
  ww(1) = 0.5d0 - h
  ww(2) = 0.5d0 + h
  ww(3) = 0.5d0 + h
  ww(4) = 0.5d0 - h
  i = 0
  do j = 1,4
    do k = 1,4
      i = i + 1
      sg(1,i) = ss(k)
      sg(2,i) = ss(j)
      sg(3,i) = ww(j)*ww(k)
    end do ! k
  end do ! i

! 5x5 integration

elseif(l.eq.5) then

  g     =  sqrt(1120.d0)
  ss(1) =  sqrt((70.d0 + g)/126.d0)
  ss(2) =  sqrt((70.d0 - g)/126.d0)
  ss(3) =  0.0d0
  ss(4) = -ss(2)
  ss(5) = -ss(1)

  ww(1) =  (21.d0*g + 117.6d0)/(g*(70.d0 + g))
  ww(2) =  (21.d0*g - 117.6d0)/(g*(70.d0 - g))
  ww(3) =  2.d0*(1.d0 - ww(1) - ww(2))
  ww(4) =  ww(2)
  ww(5) =  ww(1)

  i = 0
  do j = 1,5
    do k = 1,5
      i = i + 1
      sg(1,i) = ss(k)
      sg(2,i) = ss(j)
      sg(3,i) = ww(j)*ww(k)
    end do ! k
  end do ! j

! Error

else

  write(*,2000) l
  stop

endif

! Format

2000  format(' *ERROR* INT2D: Illegal quadrature order =',i16)

return
end subroutine int2d

subroutine shp2d(xi,eta,xl,IsNatural, shp,xsj)

implicit none

! global variables
real(8), intent( in) :: xi, eta, xl(:,:)
logical, intent( in) :: IsNatural
real(8), intent(out) :: shp(3,4), xsj
! local variables
integer :: i, j, k
real(8) :: xim, xip, etam, etap
real(8) :: jac(2,2), jacinv(2,2), shpl(2,4)

! define working variables
xim=0.25d0*(1.d0-xi); etam=0.25d0*(1.d0-eta) 
xip=0.25d0*(1.d0+xi); etap=0.25d0*(1.d0+eta)

! shape functions (N^[1-4]) evaluated at (xi, eta)
shp(3,1) = 4.d0*xim*etam; shp(3,2) = 4.d0*xip*etam
shp(3,3) = 4.d0*xip*etap; shp(3,4) = 4.d0*xim*etap

! natural derivatives of shape fuctions evaluated at (xi, eta)
shp(1,1)=-etam; shp(1,2)=etam; shp(1,3)=etap; shp(1,4)=-etap
shp(2,1)=-xim;  shp(2,2)=-xip; shp(2,3)=xip;  shp(2,4)=xim

! loop to the find the jacobian matrix
jac = 0.d0
do i = 1, 2
  do j = 1, 2
    do k = 1, 4 ! k is the node index
      jac(i,j) = jac(i,j) + shp(j,k)*xl(i,k)
    end do
  end do
end do
! compute the determinant of the jacobi matrix
xsj = jac(1,1)*jac(2,2) - jac(1,2)*jac(2,1)

! if one is looking for the global derivatives
if(.not. IsNatural) then
  ! store the local derivatives
  do i = 1, 2
    do k = 1, 4
      shpl(i,k) = shp(i,k); shp(i,k) = 0.d0
    end do
  end do
  ! inverse of the jacobi matrix
  jacinv(1,1) = jac(2,2)/xsj; jacinv(1,2) = -jac(1,2)/xsj
  jacinv(2,1) = -jac(2,1)/xsj; jacinv(2,2) = jac(1,1)/xsj
  ! get the global derivatives
  do k = 1, 4
    do i = 1, 2
      do j = 1, 2
        shp(i,k) = shp(i,k) + shpl(j,k)*jacinv(j,i)
      end do
    end do
  end do
endif

return
end subroutine shp2d

subroutine shp2ds(xi,eta,xl,IsNatural, shp,xsj,jacinv)

implicit none

! global variables
real(8), intent( in) :: xi, eta, xl(:,:)
logical, intent( in) :: IsNatural
real(8), intent(out) :: shp(3,4), jacinv(2,2), xsj
! local variables
integer :: i, j, k
real(8) :: xim, xip, etam, etap
real(8) :: jac(2,2), shpl(2,4)

! define working variables
xim=0.25d0*(1.d0-xi); etam=0.25d0*(1.d0-eta) 
xip=0.25d0*(1.d0+xi); etap=0.25d0*(1.d0+eta)

! shape functions (N^[1-4]) evaluated at (xi, eta)
shp(3,1) = 4.d0*xim*etam; shp(3,2) = 4.d0*xip*etam
shp(3,3) = 4.d0*xip*etap; shp(3,4) = 4.d0*xim*etap

! natural derivatives of shape fuctions evaluated at (xi, eta)
shp(1,1)=-etam; shp(1,2)=etam; shp(1,3)=etap; shp(1,4)=-etap
shp(2,1)=-xim;  shp(2,2)=-xip; shp(2,3)=xip;  shp(2,4)=xim

! loop to the find the jacobian matrix
jac = 0.d0
do i = 1, 2
  do j = 1, 2
    do k = 1, 4 ! k is the node index
      jac(i,j) = jac(i,j) + shp(j,k)*xl(i,k)
    end do
  end do
end do
! compute the determinant of the jacobi matrix
xsj = jac(1,1)*jac(2,2) - jac(1,2)*jac(2,1)


! inverse of the jacobi matrix
jacinv(1,1) = jac(2,2)/xsj; jacinv(1,2) = -jac(1,2)/xsj
jacinv(2,1) = -jac(2,1)/xsj; jacinv(2,2) = jac(1,1)/xsj

! if one is looking for the global derivatives
if(.not. IsNatural) then
  ! store the local derivatives
  do i = 1, 2
    do k = 1, 4
      shpl(i,k) = shp(i,k); shp(i,k) = 0.d0
    end do
  end do
  ! get the global derivatives
  do k = 1, 4
    do i = 1, 2
      do j = 1, 2
        shp(i,k) = shp(i,k) + shpl(j,k)*jacinv(j,i)
      end do
    end do
  end do
endif

return
end subroutine shp2ds

subroutine int3d(ll, lint,s)
![--.----+----.----+----.-----------------------------------------]
! Purpose: Gauss quadrature for 3-d element
!
! Inputs:
!    ll     - Number of points/direction
!
! Outputs:
!    lint   - Total number of quadrature points
!    s(4,*) - Gauss points (1-3) and weights (4)
![--.----+----.----+----.-----------------------------------------]
implicit  none

integer   i,j,k,ll,lint, ig(4),jg(4)
real*8    g, s(4,8), sg(5),wg(5)

save

data      ig/-1,1,1,-1/,jg/-1,-1,1,1/

! 1 pt. quadrature

if(ll.eq.1) then

  lint = 1
  do i = 1,3
    s(i,1) = 0.0d0
  end do ! i
  s(4,1) = 8.0d0

! 2 x 2 x 2 pt. quadrature

elseif(ll.eq.2) then

  lint = 8
  g    = 1.d0/sqrt(3.0d0)
  do i = 1,4
    s(1,i)   = ig(i)*g
    s(1,i+4) = s(1,i)
    s(2,i)   = jg(i)*g
    s(2,i+4) = s(2,i)
    s(3,i)   =  -g !  Modified here, for surface formulation
    s(3,i+4) =  g  !
    s(4,i)   = 1.d0
    s(4,i+4) = 1.d0
  end do ! i

! Error
else

  write(*,2000) ll

endif

2000  format('  *ERROR* INT3D: Illegal quadrature order =',i16)

return
end subroutine int3d

subroutine shp3d(xi,eta,zeta,xl,IsNatural, shp,xsj)

use math
implicit none

! global variables
real(8), intent( in) :: xi, eta,zeta, xl(:,:)
logical, intent( in) :: IsNatural
real(8), intent(out) :: shp(4,8), xsj
! local variables
integer :: i, j, k
real(8) :: xim, xip, etam, etap, zetam, zetap
real(8) :: jac(3,3), jacinv(3,3), shpl(3,8)
real(8), parameter :: pt125 = 0.125d0

! define working variables
xim=(1.d0-xi); etam=(1.d0-eta); zetam=(1.d0-zeta)
xip=(1.d0+xi); etap=(1.d0+eta); zetap=(1.d0+zeta)

! shape functions (N^[1-8]) evaluated at (xi, eta, zeta)
shp(4,1) = pt125*xim*etam*zetam; shp(4,2) = pt125*xip*etam*zetam
shp(4,3) = pt125*xip*etap*zetam; shp(4,4) = pt125*xim*etap*zetam
shp(4,5) = pt125*xim*etam*zetap; shp(4,6) = pt125*xip*etam*zetap
shp(4,7) = pt125*xip*etap*zetap; shp(4,8) = pt125*xim*etap*zetap

! natural derivatives of shape fuctions evaluated at (xi, eta, zeta)
shp(1,1)= -pt125*etam*zetam; shp(1,2)=  pt125*etam*zetam
shp(1,3)=  pt125*etap*zetam; shp(1,4)= -pt125*etap*zetam
shp(1,5)= -pt125*etam*zetap; shp(1,6)=  pt125*etam*zetap
shp(1,7)=  pt125*etap*zetap; shp(1,8)= -pt125*etap*zetap

shp(2,1)= -pt125*xim*zetam; shp(2,2)= -pt125*xip*zetam
shp(2,3)=  pt125*xip*zetam; shp(2,4)=  pt125*xim*zetam
shp(2,5)= -pt125*xim*zetap; shp(2,6)= -pt125*xip*zetap
shp(2,7)=  pt125*xip*zetap; shp(2,8)=  pt125*xim*zetap
 
shp(3,1)= -pt125*xim*etam; shp(3,2)= -pt125*xip*etam
shp(3,3)= -pt125*xip*etap; shp(3,4)= -pt125*xim*etap
shp(3,5)=  pt125*xim*etam; shp(3,6)=  pt125*xip*etam
shp(3,7)=  pt125*xip*etap; shp(3,8)=  pt125*xim*etap

! loop to the find the jacobian matrix
jac = 0.d0
do i = 1, 3
  do j = 1, 3
    do k = 1, 8 ! k is the node index
      jac(i,j) = jac(i,j) + shp(j,k)*xl(i,k)
    end do
  end do
end do
! compute the determinant of the jacobi matrix
xsj = determinant(jac)

! if one is looking for the global derivatives
if(.not. IsNatural) then
  ! store the local derivatives
  do i = 1, 3
    do k = 1, 8
      shpl(i,k) = shp(i,k); shp(i,k) = 0.d0
    end do
  end do
  ! inverse of the jacobi matrix
  jacinv = inverse(jac)
  ! get the global derivatives
  do k = 1, 8
    do i = 1, 3
      do j = 1, 3
        shp(i,k) = shp(i,k) + shpl(j,k)*jacinv(j,i)
      end do
    end do
  end do
endif

return
end subroutine shp3d

subroutine shp3ds(xi,eta,zeta,xl,IsNatural, shp,xsj,jacinv)

!---------------------------------------------------------------------
! Here we treat the jacobian matrix [J] as the deformation graident
! To find the surface area in the reference configuration, one may 
! use da/dA = J\sqrt{N \cdot C^{-1}N}
! where N is the initial out normal at the surface Gauss point
! cjj stores the diagonal of the matrix(J^T*J)^{-1}
!---------------------------------------------------------------------

use math
implicit none

! global variables
real(8), intent( in) :: xi, eta,zeta, xl(:,:)
logical, intent( in) :: IsNatural
real(8), intent(out) :: shp(4,8), xsj
real(8), intent(out) :: jacinv(3,3)
! local variables
integer :: i, j, k
real(8) :: xim, xip, etam, etap, zetam, zetap
real(8) :: jac(3,3), shpl(3,8)
real(8), parameter :: pt125 = 0.125d0

! define working variables
xim=(1.d0-xi); etam=(1.d0-eta); zetam=(1.d0-zeta)
xip=(1.d0+xi); etap=(1.d0+eta); zetap=(1.d0+zeta)

! shape functions (N^[1-8]) evaluated at (xi, eta, zeta)
shp(4,1) = pt125*xim*etam*zetam; shp(4,2) = pt125*xip*etam*zetam
shp(4,3) = pt125*xip*etap*zetam; shp(4,4) = pt125*xim*etap*zetam
shp(4,5) = pt125*xim*etam*zetap; shp(4,6) = pt125*xip*etam*zetap
shp(4,7) = pt125*xip*etap*zetap; shp(4,8) = pt125*xim*etap*zetap

! natural derivatives of shape fuctions evaluated at (xi, eta, zeta)
shp(1,1)= -pt125*etam*zetam; shp(1,2)=  pt125*etam*zetam
shp(1,3)=  pt125*etap*zetam; shp(1,4)= -pt125*etap*zetam
shp(1,5)= -pt125*etam*zetap; shp(1,6)=  pt125*etam*zetap
shp(1,7)=  pt125*etap*zetap; shp(1,8)= -pt125*etap*zetap

shp(2,1)= -pt125*xim*zetam; shp(2,2)= -pt125*xip*zetam
shp(2,3)=  pt125*xip*zetam; shp(2,4)=  pt125*xim*zetam
shp(2,5)= -pt125*xim*zetap; shp(2,6)= -pt125*xip*zetap
shp(2,7)=  pt125*xip*zetap; shp(2,8)=  pt125*xim*zetap
 
shp(3,1)= -pt125*xim*etam; shp(3,2)= -pt125*xip*etam
shp(3,3)= -pt125*xip*etap; shp(3,4)= -pt125*xim*etap
shp(3,5)=  pt125*xim*etam; shp(3,6)=  pt125*xip*etam
shp(3,7)=  pt125*xip*etap; shp(3,8)=  pt125*xim*etap

! loop to the find the jacobian matrix
jac = 0.d0
do i = 1, 3
  do j = 1, 3
    do k = 1, 8 ! k is the node index
      jac(i,j) = jac(i,j) + shp(j,k)*xl(i,k)
    end do
  end do
end do
! compute the determinant of the jacobi matrix
xsj = determinant(jac)

! inverse of the jacobi matrix
jacinv = inverse(jac)

! if one is looking for the global derivatives
if(.not. IsNatural) then
  ! store the local derivatives
  do i = 1, 3
    do k = 1, 8
      shpl(i,k) = shp(i,k); shp(i,k) = 0.d0
    end do
  end do
!  ! inverse of the jacobi matrix
!  jacinv = inverse(jac)
  ! get the global derivatives
  do k = 1, 8
    do i = 1, 3
      do j = 1, 3
        shp(i,k) = shp(i,k) + shpl(j,k)*jacinv(j,i)
      end do
    end do
  end do
endif

return
end subroutine shp3ds


subroutine assemble_vector(ndf,nen,numnp,icn,feint, fint)

implicit none

! global variables, input
integer :: ndf, nen, numnp
integer :: icn(nen)
real(8) :: feint(ndf*nen)
! global variables, output
real(8) :: fint(ndf*numnp)

! local variables
integer :: i, j

do i = 1, nen
  do j = 1, ndf
    fint((icn(i)-1)*ndf+j) = fint((icn(i)-1)*ndf+j) + feint((i-1)*ndf+j)
  end do
end do

return
end subroutine assemble_vector

subroutine assemble_matrix(ndf,nen,numnp,icn,ke, kint)

implicit none

! global variables, input
integer :: ndf, nen, numnp
integer :: icn(nen)
real(8) :: ke(ndf*nen,ndf*nen)
! global variables, ouput
real(8) :: kint(ndf*numnp,ndf*numnp)

! local variables
integer :: i, j, ii, jj

do i = 1, nen
  do j = 1, nen
    do ii = 1, ndf
      do jj = 1, ndf
        kint((icn(i)-1)*ndf+ii,(icn(j)-1)*ndf+jj) = &
		kint((icn(i)-1)*ndf+ii,(icn(j)-1)*ndf+jj) + ke((i-1)*ndf+ii,(j-1)*ndf+jj)
      end do
    end do
  end do
end do

return
end subroutine assemble_matrix

subroutine assemble_vectorij(ndf,nen,numnp,icn,jcn,fecontij, fcont)

implicit none

! global variables, input
integer :: ndf, nen, numnp
integer :: icn(nen), jcn(nen)
real(8) :: fecontij(2*ndf*nen)
! global variables, output
real(8) :: fcont(ndf*numnp)

! local variables
integer :: i, j

! contact force in element iel
do i = 1, nen
  do j = 1, ndf
    fcont((icn(i)-1)*ndf+j) = fcont((icn(i)-1)*ndf+j) + fecontij((i-1)*ndf+j)
    fcont((jcn(i)-1)*ndf+j) = fcont((jcn(i)-1)*ndf+j) + fecontij((i+nen-1)*ndf+j)
  end do
end do

return
end subroutine assemble_vectorij

subroutine assemble_matrixij(ndf,nen,numnp,icn,jcn,ke, kint)

implicit none

! global variables, input
integer, intent(in) :: ndf, nen, numnp
integer, intent(in) :: icn(nen), jcn(nen)
real(8), intent(in) :: ke(2*ndf*nen,2*ndf*nen)
! global variables, ouput
real(8), intent(inout) :: kint(ndf*numnp,ndf*numnp)

! local variables
integer :: i, j, ii, jj

do i = 1, nen
  do j = 1, nen
    do ii = 1, ndf
      do jj = 1, ndf
        ! (iel, iel)
        kint((icn(i)-1)*ndf+ii,(icn(j)-1)*ndf+jj) = &
		kint((icn(i)-1)*ndf+ii,(icn(j)-1)*ndf+jj) + ke((i-1)*ndf+ii,(j-1)*ndf+jj)
        ! (jel, jel)
        kint((jcn(i)-1)*ndf+ii,(jcn(j)-1)*ndf+jj) = &
		kint((jcn(i)-1)*ndf+ii,(jcn(j)-1)*ndf+jj) + ke((i+nen-1)*ndf+ii,(j+nen-1)*ndf+jj)
        ! (iel, jel)
        kint((icn(i)-1)*ndf+ii,(jcn(j)-1)*ndf+jj) = &
		kint((icn(i)-1)*ndf+ii,(jcn(j)-1)*ndf+jj) + ke((i-1)*ndf+ii,(j+nen-1)*ndf+jj)
        ! (jel, iel)
        kint((jcn(i)-1)*ndf+ii,(icn(j)-1)*ndf+jj) = &
		kint((jcn(i)-1)*ndf+ii,(icn(j)-1)*ndf+jj) + ke((j+nen-1)*ndf+ii,(j+nen-1)*ndf+jj)
      end do
    end do
  end do
end do

return
end subroutine assemble_matrixij

end module fem