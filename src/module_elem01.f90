module elem01

use fem
implicit none

contains

subroutine stif01(ndim,ndf,nen,nst,icn,d,ul,xl, feint,keint)
!
!----------------------------------------------------------------------------------------
!
!  stif01: Compute the element stiffness matrix and the internal element force vector
!          three-dimensional brick element, Neo-Hookean material model
!   
!  INPUT:
!         ndim               Number of space dimensions (=3)
!         ndf                Number of degrees of freedom per node (=3)
!         nen                Number of node per element (=8)
!         nst                Number of degrees of freedom per element(=24)
!         icn(nen)           Element connnectivity
!         d(:)               Material properties
!         ul(ndf,nen)        Nodal displacements
!         xl(ndim,nen)       Nodal coordinates in the reference configuration
!
!
!  OUTPUT:
!         keint(nst,nst)     Element stiffness
!         feint(nst)         Element internal force vector
!
!  OTHER:
!         sg(1-3,ig)         Isoparameteric coordinate of quadrature point ig
!         sg(4,ig)           weight of quadrature point ig
!
!
!----------------------------------------------------------------------------------------
!

use config, only : ngqp
implicit none

! global variables, input
integer :: ndim, ndf, nen, nst
integer :: icn(nen)
real(8) :: d(:), ul(ndf,nen), xl(ndim,nen)

! global varaibles, output
real(8) :: keint(nst, nst), feint(nst)

! local variables
integer :: i, j, ii, jj
integer :: ig, nip
real(8) :: sg(4,125), shps(4,nen)
real(8) :: xu(ndim,nen)
real(8) :: xsj, dvol
real(8) :: F(ndim,ndim), sig(6), aa(6,6)
real(8) :: bbd(3,6), kij
! real(8) :: bmat(6,24), keint_test(nst, nst)


keint = 0.d0; feint = 0.d0 !; keint_test = 0.d0
!------------------------------------------------------------------------------
! Set up parameters
!------------------------------------------------------------------------------
! integration point in each direction

call int3d(ngqp,nip,sg)
! nodal coordinates in the current configuration
xu = xl + ul
!------------------------------------------------------------------------------
! Loop over the quadrature points in the element
!------------------------------------------------------------------------------
do ig = 1, nip
  ! shape functions and their cartesian derivatives wrt the current configuration 
  call shp3d(sg(1,ig),sg(2,ig),sg(3,ig),xu,.false., shps,xsj)
  ! compute the volume element
  dvol = xsj*sg(4,ig)
  ! get the deformation gradient at ig
  call kine01(ndim,ndf,nen,shps,ul, F)
  ! get the Cauchy stress sig and spatial tagent modulus c at ig
  call modl01(ndim,d,F, sig,aa)
  ! multiply tagent moduli and Cauchy stress by the volume element
  sig = sig*dvol
  aa = aa*dvol
  !------------------------------------------------------------------------------------------------
  !! bmatrix of the element
  !bmat = 0.d0
  !do i = 1,nen
  !  bmat(1,(i-1)*3+1) = shps(1,i); bmat(2,(i-1)*3+2) = shps(2,i); bmat(3,(i-1)*3+3) = shps(3,i)
  !  bmat(4,(i-1)*3+1) = shps(2,i); bmat(4,(i-1)*3+2) = shps(1,i); bmat(5,(i-1)*3+2) = shps(3,i)
  !  bmat(5,(i-1)*3+3) = shps(2,i); bmat(6,(i-1)*3+1) = shps(3,i); bmat(6,(i-1)*3+3) = shps(1,i)
  !end do
  !keint_test = keint_test + matmul(transpose(bmat),matmul(aa,bmat))
  !------------------------------------------------------------------------------------------------
  ! the element internal force vector
  do i = 1, nen
    feint((i-1)*ndf+1) = feint((i-1)*ndf+1) + shps(1,i)*sig(1) + shps(2,i)*sig(4) + shps(3,i)*sig(6)
    feint((i-1)*ndf+2) = feint((i-1)*ndf+2) + shps(2,i)*sig(2) + shps(1,i)*sig(4) + shps(3,i)*sig(5)
    feint((i-1)*ndf+3) = feint((i-1)*ndf+3) + shps(3,i)*sig(3) + shps(2,i)*sig(5) + shps(1,i)*sig(6)
  end do
  ! accumulate the stiffness
  do i = 1, nen
    ! compute bmat-t*aa*dvol (3 x 6 matrix in 3D case)
    do jj = 1, 6
      bbd(1,jj) = shps(1,i)*aa(1,jj) + shps(2,i)*aa(4,jj) + shps(3,i)*aa(6,jj)
      bbd(2,jj) = shps(2,i)*aa(2,jj) + shps(1,i)*aa(4,jj) + shps(3,i)*aa(5,jj)
      bbd(3,jj) = shps(3,i)*aa(3,jj) + shps(2,i)*aa(5,jj) + shps(1,i)*aa(6,jj)
    end do
    ! compute material tangent stiffness
    do j = 1, i ! only need to get the lower part of the matrix, symmetric
      ! compute the factor for the geometric stiffness, corresponding to node(i,j), i.e., \int_{\Omega_t}(Ni;k)*\sigma_{kl}*(Nj;l)dv
      kij =  shps(1,i)*sig(1)*shps(1,j) + shps(2,i)*sig(2)*shps(2,j) + shps(3,i)*sig(3)*shps(3,j) &
           + shps(1,i)*sig(4)*shps(2,j) + shps(2,i)*sig(5)*shps(3,j) + shps(3,i)*sig(6)*shps(1,j) &
           + shps(2,i)*sig(4)*shps(1,j) + shps(3,i)*sig(5)*shps(2,j) + shps(1,i)*sig(6)*shps(3,j)
      ! accumulate the geometric stiffness matrix, a diagonal identity matrix with factor kij
      do ii = 1, ndf
        keint((i-1)*ndf+ii,(j-1)*ndf+ii) = keint((i-1)*ndf+ii,(j-1)*ndf+ii) + kij
      end do
      ! accumulte the material stiffness matrix, submatrix Bt*D*B
      do ii = 1, ndf
        keint((i-1)*ndf+ii,(j-1)*ndf+1) = keint((i-1)*ndf+ii,(j-1)*ndf+1) &
		+ bbd(ii,1)*shps(1,j) + bbd(ii,4)*shps(2,j) + bbd(ii,6)*shps(3,j)
        keint((i-1)*ndf+ii,(j-1)*ndf+2) = keint((i-1)*ndf+ii,(j-1)*ndf+2) &
		+ bbd(ii,2)*shps(2,j) + bbd(ii,4)*shps(1,j) + bbd(ii,5)*shps(3,j)
        keint((i-1)*ndf+ii,(j-1)*ndf+3) = keint((i-1)*ndf+ii,(j-1)*ndf+3) &
		+ bbd(ii,3)*shps(3,j) + bbd(ii,5)*shps(2,j) + bbd(ii,6)*shps(1,j)
      end do
    end do ! end of the loop over node j
  end do ! end of the loop over node i
  
!------------------------------------------------------------------------------
! End of the loop over the quadraure points in the element
!------------------------------------------------------------------------------
end do

! compute the upper part of the stiffness matrix by symmetry
do j = 1, nst
  do i = 1, j-1
    keint(i,j) = keint(j,i)
  end do
end do

return
end subroutine stif01

subroutine kine01(ndim,ndf,nen,shps,ul, F)
!
!----------------------------------------------------------------------------------------
!
!  kine01: Compute the deformation gradient F
!  
!  INPUT:
!         ndim               Number of space dimension (=3)
!         ndf                Number of degrees of freedom per node (=3)
!         nen                Number of node per element (=8)
!         ul(3,i)            Displacements of node i 
!         shps(1,i)          x-derivative of shape function of node i
!         shps(2,i)          y-derivative of shape function of node i
!         shps(3,i)          z-derivative of shape function of node i
!         shps(4,i)          shape function of node i
!
!  OUTPUT:
!         F(ndim,ndim)       Deformation gradient
!
!----------------------------------------------------------------------------------------
!
use math
implicit none

! global variables, input
integer :: ndim, ndf, nen
real(8) :: shps(4,nen), ul(ndf,nen)
! global variables, output
real(8) :: F(ndim,ndim)

! local variables
integer :: i, j, k
real(8) :: gradu(ndim,ndim)

! get the displacement gradient at ig
gradu = 0.d0
do i = 1, nen
  do j = 1, ndim
    do k = 1, ndim
      gradu(j,k) = gradu(j,k) + shps(k,i)*ul(j,i)
    end do
  end do
end do
! get the inverse deformation gradient at ig, F^{-1} = I - \nabla u
do i = 1, ndim
  do j = 1, ndim
    if(i == j) then
      F(i,j) = 1.d0 - gradu(i,j)
    else
      F(i,j) = - gradu(i,j)
    endif
  end do
end do
! invert F^{-1} to get F, overritten by F
call invert(F,ndim,ndim)

return
end subroutine kine01

subroutine modl01(ndim,d,F, sig,aa)
!
!----------------------------------------------------------------------------------------
!
!  modl01: Compute the tangent modulus and the Cauchy stress, 
!          for three-dimensional, Neo-Hookean material model
!          W = \frac{\Lambda}{2}(\ln J)^2 + \frac{mu}{2}(I_1-3) - \mu \ln J
!          S = \frac{\partial W}{\partial C} = \Lambda \ln J C^{-1} + \mu(I - C^{-1})
!          \sigma = \frac{1}{J}(\lambda \ln(J) I - \mu (b-I))
!          c = \frac{1}{J}\left[2(\mu - \Lambda\ln(J))II + \Lambda(I \otimes I)\right]
!
!  INPUT: 
!         ndim               Number of space dimension (=3)
!         d(:)               Material properties
!         F(ndim,ndim)       Deformation gradient
!
!  OUPUT:
!         aa(6,6)            Spatial tagent modulus of Neo-Hookean Model
!         sig(6)             Cauchy stress in voight notation
!
!  OTHER VARIABLES:
!         xsj                Determinant of the Jacobian
!         b(ndim,ndim)       Left Cauchy Green tensor
!         oneg               Lame constant \mu
!         xlambda            Lame constant \Lambda
!         onegp              Modified Lame constant, (\mu-\lambda*log(J))/J
!         xlambdap           Modified Lame constant, \lambda/J
!
!----------------------------------------------------------------------------------------
!
use math
implicit none
! global variables, input
integer :: ndim
real(8) :: d(:), F(ndim,ndim)
! global variables, output
real(8) :: sig(6), aa(6,6)

! local variables
real(8) :: oneg, xlambda
real(8) :: onegp, twogp, xlambdap
real(8) :: xsj, b(ndim,ndim)
real(8) :: Ev, nu

! initialization
sig = 0.d0; aa = 0.d0

! Youngs modulus and possion's ratio
Ev = d(1); nu = d(2)
! lame constants
oneg = Ev/(2.d0+2.d0*nu); xlambda = oneg*2.d0*nu/(1.d0-2.d0*nu)

! determinant of the deformation gradient
xsj = determinant(F)

if(xsj .lt. 1.0d-25) then
  write(*,*) 'Negative Jacobian in subroutine modl01, quit!'
  write(*,*) 'det[F] = ',xsj
  stop
endif

! get the modified lame constants, corresponding to the linear elastic case
! xlambdap = \lambda/J; onegp = (\mu - \lambda*log(J) )/J
xlambdap = xlambda/xsj
onegp = (oneg-xlambda*log(xsj))/xsj
twogp = onegp*2.0d0

! get the spatial tagent c_{ijkl}, in voight notation
aa(1,1) = xlambdap + twogp
aa(2,1) = xlambdap 
aa(3,1) = xlambdap 

aa(1,2) = xlambdap 
aa(2,2) = xlambdap + twogp
aa(3,2) = xlambdap 

aa(1,3) = xlambdap 
aa(2,3) = xlambdap 
aa(3,3) = xlambdap + twogp

aa(4,4) = onegp
aa(5,5) = onegp
aa(6,6) = onegp

! compute the left Cauchy Green tensor
b = matmul(F,transpose(F))

! get the Cauchy stress, in voight notation
sig(1) = oneg/xsj*b(1,1) - onegp
sig(2) = oneg/xsj*b(2,2) - onegp
sig(3) = oneg/xsj*b(3,3) - onegp
sig(4) = oneg/xsj*b(1,2)
sig(5) = oneg/xsj*b(2,3)
sig(6) = oneg/xsj*b(1,3)


return
end subroutine modl01


subroutine stre01(ndim,ndf,nen,nst,icn,d,ul,xl, lpm,stp)
!
!----------------------------------------------------------------------------------------
!
!  stre01: get the elemental contribution to the lumped projection matrix lpm and the
!          projected stresses stp(to be divided by lpm)
!   
!  INPUT:
!         ndim               Number of space dimensions (=3)
!         ndf                Number of degrees of freedom per node (=3)
!         nen                Number of node per element (=8)
!         nst                Number of degrees of freedom per element(=24)
!         icn(nen)           Element connnectivity
!         d(:)               Material properties
!         ul(ndf,nen)        Nodal displacements
!         xl(ndim,nen)       Nodal coordinates in the reference configuration
!
!
!  OUTPUT:
!         keint(nst,nst)     Element stiffness
!         feint(nst)         Element internal force vector
!
!  OTHER:
!         sg(1-3,ig)         Isoparameteric coordinate of quadrature point ig
!         sg(4,ig)           weight of quadrature point ig
!
!
!----------------------------------------------------------------------------------------
!

use config, only : ngqp
implicit none

! global variables, input
integer :: ndim, ndf, nen, nst
integer :: icn(nen)
real(8) :: d(:), ul(ndf,nen), xl(ndim,nen)

! global varaibles, output
real(8) :: lpm(:), stp(:,:)

! local variables
integer :: i, j, ii, jj
integer :: ig, nip
real(8) :: sg(4,125), shps(4,nen)
real(8) :: xu(ndim,nen)
real(8) :: xsj, dvol
real(8) :: F(ndim,ndim), sig(6), aa(6,6)
real(8) :: shear, normal

!------------------------------------------------------------------------------
! Set up parameters
!------------------------------------------------------------------------------
! integration point in each direction
call int3d(ngqp,nip,sg)
! nodal coordinates in the current configuration
xu = xl + ul
!------------------------------------------------------------------------------
! Loop over the quadrature points in the element
!------------------------------------------------------------------------------
do ig = 1, nip
  ! shape functions and their cartesian derivatives wrt the current configuration 
  call shp3d(sg(1,ig),sg(2,ig),sg(3,ig),xu,.false., shps,xsj)
  ! compute the volume element
  dvol = xsj*sg(4,ig)
  ! get the deformation gradient at ig
  call kine01(ndim,ndf,nen,shps,ul, F)
  ! get the Cauchy stress sig and spatial tagent modulus c at ig
  call modl01(ndim,d,F, sig,aa)
  ! get the von-mises (\tau) and the normal stress (tr(\sigma))
  shear =  sqrt(((sig(1)-sig(2))**2 + (sig(2)-sig(3))**2 + (sig(1)-sig(3))**2 &
              + 6.d0*(sig(4)**2 + sig(5)**2 + sig(6)**2))/2.d0)
  normal = sig(1) + sig(2) + sig(3)
  ! lumped projection matrix lpm(numnp), and projected stresse( p to division by lmp)
  do i = 1, nen
    j = icn(i)
    lpm(j) = lpm(j) + shps(4,i)*dvol
    stp(1,j) = stp(1,j) + shps(4,i)*shear*dvol
    stp(2,j) = stp(2,j) + shps(4,i)*normal*dvol
  end do
end do

return
end subroutine stre01


subroutine stifcont01(ndim,ndf,nen,nip,icn,jcn, fecontij,kecontij)

use math
use config, only : vcntp, xr, disp, ngqp
implicit none

! global variables, input
integer :: icn(8), jcn(8)

! global variables, output
real(8) :: kecontij(:,:), fecontij(:)

integer :: ndim, ndf, nen, nip
integer :: ig, jg, i, j, ii, jj
real(8) :: xsj, dvoli, dvolj
real(8) :: sg(4,125), shpsi(4,8), shpsj(4,8)
real(8) :: xli(3,8), uli(3,8), xui(3,8)
real(8) :: xlj(3,8), ulj(3,8), xuj(3,8)
real(8) :: xig(3), xjg(3), rij(3)
real(8) :: r, r0, betai, betaj, epsilon, tmp, dfdr,df2dr2
real(8) :: forceij(3), kii(3,3)

!------------------------------------------------------------------------------
! Set up parameters
!------------------------------------------------------------------------------
! contact parameters
betai = vcntp(1); betaj = vcntp(2)
r0 = vcntp(3); epsilon = vcntp(4)
! integration point in each direction
call int3d(ngqp,nip,sg)

!----------------------------------------------------------------------------
! information related to element iel
!----------------------------------------------------------------------------
! coordinate and displacement
xli = xr(:, icn)
do i = 1, nen
  do j = 1, ndf
    uli(j,i) = disp((icn(i)-1)*ndf + j)
  end do
end do
xui = xli + uli

!----------------------------------------------------------------------------
! information related to element jel
!----------------------------------------------------------------------------
! coordinate and displacement
xlj = xr(:, jcn)
do i = 1, nen
  do j = 1, ndf
    ulj(j,i) = disp((jcn(i)-1)*ndf + j)
  end do
end do
xuj = xlj + ulj

!------------------------------------------------------------------------------
! Loop over the quadraure points in the element iel
!------------------------------------------------------------------------------
do ig = 1, nip
  ! shape functions and their cartesian derivatives wrt the current configuration 
  call shp3d(sg(1,ig),sg(2,ig),sg(3,ig),xui,.false., shpsi,xsj)
  ! volume element associated with ig
  dvoli = xsj*sg(4,ig)
  ! current location of quadrature point ig
  xig = 0.d0
  do i = 1, ndim
    do j = 1, nen
      xig(i) = xig(i) + shpsi(4,j)*xui(i,j) 
    end do
  end do
  ! loop over the quadrature points in jel
  do jg = 1, nip
    ! shape functions and their cartesian derivatives wrt the current configuration 
    call shp3d(sg(1,jg),sg(2,jg),sg(3,jg),xuj,.false., shpsj,xsj)
    ! volume element associated with jg
    dvolj = xsj*sg(4,jg)
    ! current location of quadrature point jg
    xjg = 0.d0
    do i = 1, ndim
      do j = 1, nen
        xjg(i) = xjg(i) + shpsj(4,j)*xuj(i,j) 
      end do
    end do
    ! vector of the two quadrature points
    rij = xig - xjg
    ! distance between the two quadrature points
    r = sqrt(dot_product(rij,rij))
    ! interaction only happens within the cutoff distance
    if(r >= 3.d0*r0) cycle
    ! intermediate variable
    tmp = (r0/r)**6
    ! first derivative of the potential function
    dfdr = 12.d0*epsilon*tmp*(-tmp + 1.d0)/r
    ! force vector \beta_i\beta_j\frac{\partial \phi}{\partial x_i}dv_j dv_i
    forceij(:) = betai*betaj*dfdr*rij(:)*dvoli*dvolj/r
    ! distribute forceij, -forceij to the nodes of element iel, jel, respectively
    do i = 1, nen
      do j = 1, ndf
        fecontij((i-1)*ndf+j) = fecontij((i-1)*ndf+j) + shpsi(4,i)*forceij(j)
        fecontij((i+nen-1)*ndf+j) = fecontij((i+nen-1)*ndf+j) - shpsj(4,i)*forceij(j)
      end do
    end do
    ! second derivative of the potential function
    df2dr2 = 12.d0*epsilon*tmp*(13.d0*tmp-7.d0)/r/r
    ! stiffness matrix \beta_i\beta_j\frac{\partial^2 \phi}{\partial x_i \partial x_i}dv_j dv_i 
    kii = (df2dr2 - dfdr/r)*tensor_product(rij/r,rij/r)
    do i = 1, ndim
      kii(i,i) = kii(i,i) + dfdr/r
    end do  
    kii = betai*betaj*dvoli*dvolj*kii
    ! contribution to the element contact stiffness matrix along the diagonal, need summation over jel
    ! part1, (iel,iel), (jel, jel), (iel, jel), (jel, iel)
    do i = 1, nen
      do j = 1, nen
        do ii = 1, ndf
          do jj = 1, ndf
            kecontij((i-1)*ndf+ii, (j-1)*ndf+jj) = &
			kecontij((i-1)*ndf+ii, (j-1)*ndf+jj) + shpsi(4,i)*shpsi(4,j)*kii(ii,jj)
            kecontij((i+nen-1)*ndf+ii, (j+nen-1)*ndf+jj) = &
			kecontij((i+nen-1)*ndf+ii, (j+nen-1)*ndf+jj) + shpsj(4,i)*shpsj(4,j)*kii(ii,jj)
            kecontij((i-1)*ndf+ii, (j+nen-1)*ndf+jj) = &
			kecontij((i-1)*ndf+ii, (j+nen-1)*ndf+jj) - shpsi(4,i)*shpsj(4,j)*kii(ii,jj)
            kecontij((i+nen-1)*ndf+ii, (j-1)*ndf+jj) = &
			kecontij((i+nen-1)*ndf+ii, (j-1)*ndf+jj) - shpsi(4,j)*shpsj(4,i)*kii(ii,jj)
          end do
        end do
      end do
    end do
    !
  end do ! end of jg
!------------------------------------------------------------------------------
! end of the loop over the quadraure points in the element iel
!------------------------------------------------------------------------------
end do ! end of ig


return
end subroutine stifcont01

subroutine getsurfinfo !(numels,stob,srfid,cn,xr, dA,sN,dNdX)
!
! purpose: 1. Obtain the unit surface norm N and the associated area dA numerically, 
!          based on the Nanson's Formula
!          2. Compute and store the derivative of the surface out normal wrt the
!          reference configuration---dN/dX
!
use math
use config , only : numels, stob, srfid, cn, xr, dA, sN, dNdX, ngqps
implicit none

!! global variables, input
!integer :: numels
!integer :: stob(numels), srfid(numels)
!integer :: cn(:,:)
!real(8) :: xr(:,:)
!! global variables, output
!real(8) :: dA(:,:), sN(:,:,:), dNdX(:,:,:,:)

! local variables
integer :: i, j, k, ii, jj, kk
integer :: iels, iel, igs, nips, icn(8)
real(8) :: xli(3,8), points(2,25), weights(25), n(3), xc, xx(3)
real(8) :: sg(3), shpsi(4,8), xsj, jacinv(3,3), cc(3,3)


! get the surface quadrature points and integration weight
call sint3d(ngqps, points,weights); nips = ngqps*ngqps
!------------------------------------------------------------------------------
! Loop over the surface element numels in the system
!------------------------------------------------------------------------------
do iels = 1, numels
  ! corresponding bulk element index
  iel = stob(iels)
  ! bulk element connectivity
  icn = cn(:,iel)
  ! coordinate
  xli = xr(:, icn)
  !------------------------------------------------------------------------------
  ! Loop over the quadrature points in the surface element iels
  !------------------------------------------------------------------------------
  do igs = 1, nips
    ! coordinate of the surface quadrature point igs (of dimension [ndim-1]) in the ndim space
    select case(srfid(iels)) 
      case(1)
        sg(1) = points(1,igs); sg(2) = points(2,igs); sg(3) = -1.d0 ! face 1-4-3-2
      case(2)
        sg(1) = points(1,igs); sg(2) = points(2,igs); sg(3) =  1.d0 ! face 5-6-7-8
      case(3)
        sg(1) = points(1,igs); sg(3) = points(2,igs); sg(2) = -1.d0 ! face 1-2-6-5
      case(4)
        sg(1) = points(1,igs); sg(3) = points(2,igs); sg(2) =  1.d0 ! face 3-4-7-8
      case(5)
        sg(2) = points(1,igs); sg(3) = points(2,igs); sg(1) = -1.d0 ! face 1-5-8-4
      case(6)
        sg(2) = points(1,igs); sg(3) = points(2,igs); sg(1) =  1.d0 ! face 2-3-7-6
      case default
        write(*,*) 'Wrong surface element ID in getsurfinfo, Quit!'
        stop
    end select
    ! shape functions and their cartesian derivatives wrt the reference configuration 
    call shp3ds(sg(1),sg(2),sg(3),xli,.true., shpsi,xsj,jacinv)
    ! inverse of C = inv(jac^T*jac)
    cc = matmul(jacinv,transpose(jacinv))
    ! area element associated with igs
    select case(srfid(iels)) ! n = F^{-T}*N/{*}
      case(1)
        n = -jacinv(3,:); sN(:,igs,iels) = normalize(n)
        dA(igs,iels) = xsj*weights(igs)*sqrt(abs(cc(3,3)))
      case(2)
        n =  jacinv(3,:); sN(:,igs,iels) = normalize(n)
        dA(igs,iels) = xsj*weights(igs)*sqrt(abs(cc(3,3)))
      case(3)
        n = -jacinv(2,:); sN(:,igs,iels) = normalize(n)
        dA(igs,iels) = xsj*weights(igs)*sqrt(abs(cc(2,2)))
      case(4)
        n =  jacinv(2,:); sN(:,igs,iels) = normalize(n)
        dA(igs,iels) = xsj*weights(igs)*sqrt(abs(cc(2,2)))
      case(5)
        n = -jacinv(1,:); sN(:,igs,iels) = normalize(n)
        dA(igs,iels) = xsj*weights(igs)*sqrt(abs(cc(1,1)))
      case(6)
        n =  jacinv(1,:); sN(:,igs,iels) = normalize(n)
        dA(igs,iels) = xsj*weights(igs)*sqrt(abs(cc(1,1)))
      case default
        write(*,*) 'Wrong surface element ID in getsurfinfo, Quit!'
        stop
    end select
  end do
!------------------------------------------------------------------------------
! End of the loop over the surface element numels in the system
!------------------------------------------------------------------------------
end do

return
end subroutine getsurfinfo

subroutine sint3d(ngp, points,weights)

implicit none

! global variables, input
integer :: ngp

! global variables, output
real(8) :: points(2,ngp*ngp), weights(ngp*ngp)


! local variables
integer :: i, j, ik
real(8) :: x(ngp), w(ngp)

select case(ngp)
  case(1)
    x(1) = 0.d0; w(1) = 2.d0
  case(2)
    x(1) = -1.d0/sqrt(3.d0); w(1) = 1.d0
    x(2) =  1.d0/sqrt(3.d0); w(2) = 1.d0
  case(3)
    x(1) = -sqrt(3.d0/5.d0); w(1) = 5.d0/9.d0
    x(2) = 0.d0; w(2) = 8.d0/9.d0
    x(3) = sqrt(3.d0/5.d0); w(3) = 5.d0/9.d0
  case(4)
    x(1) = -sqrt((3.d0+2.d0*sqrt(6.d0/5.d0))/7.d0); w(1) = (18.d0-sqrt(30.d0))/36.d0
    x(2) = -sqrt((3.d0-2.d0*sqrt(6.d0/5.d0))/7.d0); w(2) = (18.d0+sqrt(30.d0))/36.d0
    x(3) =  sqrt((3.d0-2.d0*sqrt(6.d0/5.d0))/7.d0); w(3) = (18.d0+sqrt(30.d0))/36.d0
    x(4) =  sqrt((3.d0+2.d0*sqrt(6.d0/5.d0))/7.d0); w(4) = (18.d0-sqrt(30.d0))/36.d0
  case(5)
    x(1) = -sqrt(5+2.d0*sqrt(10.d0/7.d0))/3.d0; w(1) = (322.d0-13.d0*sqrt(70.d0))/900.d0
    x(2) = -sqrt(5-2.d0*sqrt(10.d0/7.d0))/3.d0; w(2) = (322.d0+13.d0*sqrt(70.d0))/900.d0
    x(3) = 0.d0                               ; w(3) = 128.d0/225.d0
    x(4) =  sqrt(5-2.d0*sqrt(10.d0/7.d0))/3.d0; w(4) = (322.d0+13.d0*sqrt(70.d0))/900.d0
    x(5) =  sqrt(5+2.d0*sqrt(10.d0/7.d0))/3.d0; w(5) = (322.d0-13.d0*sqrt(70.d0))/900.d0
  case default
    write(*,*) 'Not Supported ngp = ', ngp
end select

ik = 0
do i = 1, ngp
  do j = 1, ngp
    ik = ik + 1
    points(1,ik) = x(i); points(2,ik) = x(j); weights(ik) = w(i)*w(j)
  end do
end do


return
end subroutine sint3d


subroutine stifconts01(ndim,ndf,nen,nip,icn,jcn,iels,jels, fecontij,kecontij)

use math
use config, only : vcntp, xr, disp, srfid, dA, sN, dNdX, ngqps, ngqp
implicit none

! global variables, input
integer :: iels, jels
integer :: icn(nen), jcn(nen)

! global variables, output
real(8) :: kecontij(:,:), fecontij(:)

integer :: ndim, ndf, nen, nip
integer :: igs, jgs, i, j, k, ii, jj, kk
real(8) :: xsj, dvoli, dvolj
real(8) :: sg(3), shpsi(4,8), shpsj(4,8)
real(8) :: xli(3,8), uli(3,8), xui(3,8)
real(8) :: xlj(3,8), ulj(3,8), xuj(3,8)
real(8) :: xigs(3), xjgs(3), rij(3)
real(8) :: r, r0, betai, betaj, epsilon, tmp, psi, dpsidr, factor

integer :: nips
real(8) :: points(2,25), weights(25)
real(8) :: Fi(3,3), Fj(3,3), invFit(3,3), invFjt(3,3)
real(8) :: sNi(3), sNj(3), dNidX(3,3), dNjdX(3,3), dAi, dAj
real(8) :: forcei(3), forcej(3)
real(8) :: kii(3,3), kij(3,3), kjj(3,3), kji(3,3)
!------------------------------------------------------------------------------
! Set up parameters
!------------------------------------------------------------------------------
! contact parameters
betai = vcntp(1); betaj = vcntp(2)
r0 = vcntp(3); epsilon = vcntp(4)
! integration point in each direction
call sint3d(ngqps,points,weights); nips = ngqps*ngqps
!----------------------------------------------------------------------------
! information related to surface element iels, belonging to iel
!----------------------------------------------------------------------------
! coordinate and displacement
xli = xr(:, icn)
do i = 1, nen
  do j = 1, ndf
    uli(j,i) = disp((icn(i)-1)*ndf + j)
  end do
end do
xui = xli + uli

!----------------------------------------------------------------------------
! information related to surface element jels, belonging to jel
!----------------------------------------------------------------------------
! coordinate and displacement
xlj = xr(:, jcn)
do i = 1, nen
  do j = 1, ndf
    ulj(j,i) = disp((jcn(i)-1)*ndf + j)
  end do
end do
xuj = xlj + ulj

!------------------------------------------------------------------------------
! Loop over the quadraure points in the surface element iels
!------------------------------------------------------------------------------
do igs = 1, nips
  ! coordinate of the surface quadrature point igs (of dimension [ndim-1]) in the ndim space
  select case(srfid(iels)) 
    case(1)
      sg(1) = points(1,igs); sg(2) = points(2,igs); sg(3) = -1.d0 ! face 1-4-3-2
    case(2)
      sg(1) = points(1,igs); sg(2) = points(2,igs); sg(3) =  1.d0 ! face 5-6-7-8
    case(3)
      sg(1) = points(1,igs); sg(3) = points(2,igs); sg(2) = -1.d0 ! face 1-2-6-5
    case(4)
      sg(1) = points(1,igs); sg(3) = points(2,igs); sg(2) =  1.d0 ! face 3-4-7-8
    case(5)
      sg(2) = points(1,igs); sg(3) = points(2,igs); sg(1) = -1.d0 ! face 1-5-8-4
    case(6)
      sg(2) = points(1,igs); sg(3) = points(2,igs); sg(1) =  1.d0 ! face 2-3-7-6
    case default
      write(*,*) 'Wrong surface element ID in stifconts01, Quit!'
      stop
  end select
  ! shape functions and their cartesian derivatives wrt the reference configuration 
  call shp3d(sg(1),sg(2),sg(3),xli,.false., shpsi,xsj)
  ! current location of surface quadrature point igs
  xigs = 0.d0
  do i = 1, ndim
    do j = 1, nen
      xigs(i) = xigs(i) + shpsi(4,j)*xui(i,j) 
    end do
  end do
  ! get the deformation gradient at igs
  Fi = 0.d0
  do i = 1, nen
    do j = 1, ndim
      do k = 1, ndim
        Fi(j,k) = Fi(j,k) + shpsi(k,i)*xui(j,i)
      end do
    end do
  end do
  ! inverse transpose of the deformation gradient
  invFit = transpose(inverse(Fi))
  ! surface element area, unit normal, and its first derivative
  dAi = dA(igs,iels); sNi(:) = sN(:,igs,iels); dNidX(:,:) = dNdX(:,:,igs,iels)
  ! loop over the surface quadrature points in jels
  do jgs = 1, nips
    ! coordinate of the surface quadrature point jgs (of dimension [ndim-1]) in the ndim space
    select case(srfid(jels)) 
      case(1)
        sg(1) = points(1,jgs); sg(2) = points(2,jgs); sg(3) = -1.d0 ! face 1-4-3-2
      case(2)
        sg(1) = points(1,jgs); sg(2) = points(2,jgs); sg(3) =  1.d0 ! face 5-6-7-8
      case(3)
        sg(1) = points(1,jgs); sg(3) = points(2,jgs); sg(2) = -1.d0 ! face 1-2-6-5
      case(4)
        sg(1) = points(1,jgs); sg(3) = points(2,jgs); sg(2) =  1.d0 ! face 3-4-7-8
      case(5)
        sg(2) = points(1,jgs); sg(3) = points(2,jgs); sg(1) = -1.d0 ! face 1-5-8-4
      case(6)
        sg(2) = points(1,jgs); sg(3) = points(2,jgs); sg(1) =  1.d0 ! face 2-3-7-6
      case default
        write(*,*) 'Wrong surface element ID in stifconts01, Quit!'
        stop
    end select
    ! shape functions and their cartesian derivatives wrt the reference configuration 
    call shp3d(sg(1),sg(2),sg(3),xlj,.false., shpsj,xsj)
    ! current location of quadrature point jgs
    xjgs = 0.d0
    do i = 1, ndim
      do j = 1, nen
        xjgs(i) = xjgs(i) + shpsj(4,j)*xuj(i,j) 
      end do
    end do
    ! vector of the two quadrature points
    rij = xigs - xjgs
    ! distance between the two quadrature points
    r = sqrt(dot_product(rij,rij))
    ! interaction only happens within the cutoff distance
    if(r >= 3.d0*r0) cycle
    ! get the deformation gradient at igs
    Fj = 0.d0
    do i = 1, nen
      do j = 1, ndim
        do k = 1, ndim
          Fj(j,k) = Fj(j,k) + shpsj(k,i)*xuj(j,i)
        end do
      end do
    end do
    ! inverse of the deormaiton graident
    invFjt = transpose(inverse(Fj))
    ! surface element area, unit normal, and its first derivative
    dAj = dA(jgs,jels); sNj(:) = sN(:,jgs,jels); dNjdX(:,:) = dNdX(:,:,jgs,jels)
    ! intermediate variable
    tmp = (r0/r)**6
    ! obtain \psi(s) = , to be changed
    psi = 2.d0/3.d0*epsilon*tmp*(1.d0/6.d0*tmp-1.d0)
    ! force vector \beta_i\beta_j\frac{\partial \phi}{\partial x_i}dv_j dv_i
    factor = betai*betaj*psi*dAi*dAj !*detFi*detFj
    forcei(:) = factor*dot_product( rij,matmul(invFit,sNi))*matmul(invFjt,sNj)
    forcej(:) = factor*dot_product(-rij,matmul(invFjt,sNj))*matmul(invFit,sNi)
    ! distribute forcei, forcej to the nodes of element iel, jel, respectively
    do i = 1, nen
      do j = 1, ndf
        fecontij((i-1)*ndf+j) = fecontij((i-1)*ndf+j) + shpsi(4,i)*forcei(j)
        fecontij((i+nen-1)*ndf+j) = fecontij((i+nen-1)*ndf+j) + shpsj(4,i)*forcej(j)
      end do
    end do
    ! first derivative of the potential function psi(r)
    dpsidr = 4.d0*epsilon*tmp*(-tmp/3.d0 + 1.d0)/r
    ! stiffness kii
    kii = psi*tensor_product(matmul(invFjt,sNj), matmul(invFit,sNi)+ &
	matmul( matmul(invFit,transpose(dNidX)),matmul(transpose(invFit),rij)))
    kii = kii + dpsidr/r*dot_product(matmul(invFit,sNi),rij)*tensor_product(matmul(invFjt,sNj),rij(:))
    kii = betai*betaj*dAi*dAj*kii
    ! stiffness kij
    kij = -psi*tensor_product(matmul(invFjt,sNj), matmul(invFit,sNi)) 
    kij = kij + psi*dot_product(matmul(invFit,sNi),rij)*matmul(matmul(invFjt,dNjdX),transpose(invFjt))
    kij = kij + dpsidr/r*dot_product(matmul(invFit,sNi),rij)*tensor_product(matmul(invFjt,sNj),-rij(:))
    kij = betai*betaj*dAi*dAj*kij
    ! stiffness kji
    kji = -psi*tensor_product(matmul(invFit,sNi), matmul(invFjt,sNj)) 
    kji = kji + psi*dot_product(matmul(invFjt,sNj),-rij)*matmul(matmul(invFit,dNidX),transpose(invFit))
    kji = kji + dpsidr/r*dot_product(matmul(invFjt,sNj),-rij)*tensor_product(matmul(invFit,sNi),rij(:))
    kji = betai*betaj*dAi*dAj*kji
    ! stiffness kjj
    kjj = psi*tensor_product(matmul(invFit,sNi), matmul(invFjt,sNj)+&
	matmul( matmul(invFjt,transpose(dNjdX)),matmul(transpose(invFjt),-rij)))
    kjj = kjj + dpsidr/r*dot_product(matmul(invFjt,sNj),-rij)*tensor_product(matmul(invFit,sNi),-rij(:))
    kjj = betai*betaj*dAi*dAj*kjj
    ! contribution to the element contact stiffness matrix 
    do i = 1, nen
      do j = 1, nen
        do ii = 1, ndf
          do jj = 1, ndf
            kecontij((i-1)*ndf+ii, (j-1)*ndf+jj) = &
			kecontij((i-1)*ndf+ii, (j-1)*ndf+jj) + shpsi(4,i)*shpsi(4,j)*kii(ii,jj)
            kecontij((i+nen-1)*ndf+ii, (j+nen-1)*ndf+jj) = &
			kecontij((i+nen-1)*ndf+ii, (j+nen-1)*ndf+jj) + shpsj(4,i)*shpsj(4,j)*kjj(ii,jj)
            kecontij((i-1)*ndf+ii, (j+nen-1)*ndf+jj) = &
			kecontij((i-1)*ndf+ii, (j+nen-1)*ndf+jj) + shpsi(4,i)*shpsj(4,j)*kij(ii,jj)
            kecontij((i+nen-1)*ndf+ii, (j-1)*ndf+jj) = &
			kecontij((i+nen-1)*ndf+ii, (j-1)*ndf+jj) + shpsj(4,i)*shpsi(4,j)*kji(ii,jj)
          end do
        end do
      end do
    end do
    !
  end do ! end of jgs
!------------------------------------------------------------------------------
! end of the loop over the quadraure points in the element iels
!------------------------------------------------------------------------------
end do ! end of igs

return
end subroutine stifconts01


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  ROUTINES FOR DYNAMIC  SIMULATIONS
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine mass01(ndim,nen,d,xl, mm)

use config, only : ngqp

implicit none

! global variables, input
integer :: ndim, nen
real(8) :: d(:),  xl(ndim,nen)

! global, output
real(8) :: mm(nen)

! local variables
integer :: i, j, ii, jj
integer :: ig, nip
real(8) :: sg(4,125), shps(4,nen)
real(8) :: xsj, dvm, rho

!------------------------------------------------------------------------------
! Set up parameters
!------------------------------------------------------------------------------
! density
rho = d(3)
! integration point in each direction
call int3d(ngqp,nip,sg)
! nodal coordinates in the reference configuration
!------------------------------------------------------------------------------
! Loop over the quadrature points in the element
!------------------------------------------------------------------------------
do ig = 1, nip
  ! shape functions and their cartesian derivatives wrt the reference configuration 
  call shp3d(sg(1,ig),sg(2,ig),sg(3,ig),xl,.true., shps,xsj)
  ! compute mass element
  dvm = xsj*sg(4,ig)*rho
  ! accumulate the stiffness
  do i = 1, nen
    mm(i) = mm(i) + shps(4,i)*dvm
  end do ! end of the loop over node i
!------------------------------------------------------------------------------
! End of the loop over the quadraure points in the element
!------------------------------------------------------------------------------
end do

return
end subroutine mass01


subroutine fconts01(ndim,ndf,nen,nip,icn,jcn,iels,jels, fecontij)

use math
use config, only : vcntp, xr, disp, srfid, dA, sN, dNdX, ngqps, ngqp, energy_cont
implicit none

! global variables, input
integer :: iels, jels
integer :: icn(nen), jcn(nen)

! global variables, output
real(8) :: fecontij(:)

integer :: ndim, ndf, nen, nip
integer :: igs, jgs, i, j, k, ii, jj, kk
real(8) :: xsj, dvoli, dvolj
real(8) :: sg(3), shpsi(4,8), shpsj(4,8)
real(8) :: xli(3,8), uli(3,8), xui(3,8)
real(8) :: xlj(3,8), ulj(3,8), xuj(3,8)
real(8) :: xigs(3), xjgs(3), rij(3)
real(8) :: r, r0, betai, betaj, epsilon, tmp, psi, dpsidr, factor

integer :: nips
real(8) :: points(2,25), weights(25)
real(8) :: Fi(3,3), Fj(3,3), invFit(3,3), invFjt(3,3)
real(8) :: sNi(3), sNj(3), dNidX(3,3), dNjdX(3,3), dAi, dAj
real(8) :: forcei(3), forcej(3)
real(8) :: ur
!------------------------------------------------------------------------------
! Set up parameters
!------------------------------------------------------------------------------
! contact parameters
betai = vcntp(1); betaj = vcntp(2)
r0 = vcntp(3); epsilon = vcntp(4)
! integration point in each direction
call sint3d(ngqps,points,weights); nips = ngqps*ngqps
!----------------------------------------------------------------------------
! information related to surface element iels, belonging to iel
!----------------------------------------------------------------------------
! coordinate and displacement
xli = xr(:, icn)
do i = 1, nen
  do j = 1, ndf
    uli(j,i) = disp((icn(i)-1)*ndf + j)
  end do
end do
xui = xli + uli

!----------------------------------------------------------------------------
! information related to surface element jels, belonging to jel
!----------------------------------------------------------------------------
! coordinate and displacement
xlj = xr(:, jcn)
do i = 1, nen
  do j = 1, ndf
    ulj(j,i) = disp((jcn(i)-1)*ndf + j)
  end do
end do
xuj = xlj + ulj

!------------------------------------------------------------------------------
! Loop over the quadraure points in the surface element iels
!------------------------------------------------------------------------------
do igs = 1, nips
  ! coordinate of the surface quadrature point igs (of dimension [ndim-1]) in the ndim space
  select case(srfid(iels)) 
    case(1)
      sg(1) = points(1,igs); sg(2) = points(2,igs); sg(3) = -1.d0 ! face 1-4-3-2
    case(2)
      sg(1) = points(1,igs); sg(2) = points(2,igs); sg(3) =  1.d0 ! face 5-6-7-8
    case(3)
      sg(1) = points(1,igs); sg(3) = points(2,igs); sg(2) = -1.d0 ! face 1-2-6-5
    case(4)
      sg(1) = points(1,igs); sg(3) = points(2,igs); sg(2) =  1.d0 ! face 3-4-7-8
    case(5)
      sg(2) = points(1,igs); sg(3) = points(2,igs); sg(1) = -1.d0 ! face 1-5-8-4
    case(6)
      sg(2) = points(1,igs); sg(3) = points(2,igs); sg(1) =  1.d0 ! face 2-3-7-6
    case default
      write(*,*) 'Wrong surface element ID in stifconts01, Quit!'
      stop
  end select
  ! shape functions and their cartesian derivatives wrt the reference configuration 
  call shp3d(sg(1),sg(2),sg(3),xli,.false., shpsi,xsj)
  ! current location of surface quadrature point igs
  xigs = 0.d0
  do i = 1, ndim
    do j = 1, nen
      xigs(i) = xigs(i) + shpsi(4,j)*xui(i,j) 
    end do
  end do
  ! get the deformation gradient at igs
  Fi = 0.d0
  do i = 1, nen
    do j = 1, ndim
      do k = 1, ndim
        Fi(j,k) = Fi(j,k) + shpsi(k,i)*xui(j,i)
      end do
    end do
  end do
  ! inverse transpose of the deformation gradient
  invFit = transpose(inverse(Fi))
  ! surface element area, unit normal, and its first derivative
  dAi = dA(igs,iels); sNi(:) = sN(:,igs,iels) !; dNidX(:,:) = dNdX(:,:,igs,iels)
  ! loop over the surface quadrature points in jels
  do jgs = 1, nips
    ! coordinate of the surface quadrature point jgs (of dimension [ndim-1]) in the ndim space
    select case(srfid(jels)) 
      case(1)
        sg(1) = points(1,jgs); sg(2) = points(2,jgs); sg(3) = -1.d0 ! face 1-4-3-2
      case(2)
        sg(1) = points(1,jgs); sg(2) = points(2,jgs); sg(3) =  1.d0 ! face 5-6-7-8
      case(3)
        sg(1) = points(1,jgs); sg(3) = points(2,jgs); sg(2) = -1.d0 ! face 1-2-6-5
      case(4)
        sg(1) = points(1,jgs); sg(3) = points(2,jgs); sg(2) =  1.d0 ! face 3-4-7-8
      case(5)
        sg(2) = points(1,jgs); sg(3) = points(2,jgs); sg(1) = -1.d0 ! face 1-5-8-4
      case(6)
        sg(2) = points(1,jgs); sg(3) = points(2,jgs); sg(1) =  1.d0 ! face 2-3-7-6
      case default
        write(*,*) 'Wrong surface element ID in stifconts01, Quit!'
        stop
    end select
    ! shape functions and their cartesian derivatives wrt the reference configuration 
    call shp3d(sg(1),sg(2),sg(3),xlj,.false., shpsj,xsj)
    ! current location of quadrature point jgs
    xjgs = 0.d0
    do i = 1, ndim
      do j = 1, nen
        xjgs(i) = xjgs(i) + shpsj(4,j)*xuj(i,j) 
      end do
    end do
    ! vector of the two quadrature points
    rij = xigs - xjgs
    ! distance between the two quadrature points
    r = sqrt(dot_product(rij,rij))
    ! interaction only happens within the cutoff distance
    if(r >= 3.d0*r0) cycle
    ! get the deformation gradient at igs
    Fj = 0.d0
    do i = 1, nen
      do j = 1, ndim
        do k = 1, ndim
          Fj(j,k) = Fj(j,k) + shpsj(k,i)*xuj(j,i)
        end do
      end do
    end do
    ! inverse of the deormaiton graident
    invFjt = transpose(inverse(Fj))
    ! surface element area, unit normal, and its first derivative
    dAj = dA(jgs,jels); sNj(:) = sN(:,jgs,jels) 
    ! intermediate variable
    tmp = (r0/r)**6
    ! obtain \psi(s) = , to be changed
    psi = 2.d0/3.d0*epsilon*tmp*(1.d0/6.d0*tmp-1.d0)
    ! check for penetration
    if(r <= r0) then
      if(dot_product(rij,matmul(invFjt,sNj)) .lt. 1.d-25) then
        psi = -psi
      endif
    endif
    ! obtain ur, for contact interaction energy
    ur = epsilon*r0*r0/6.d0*(-tmp/15.d0+1.d0)*(r0/r)**4
    ! force vector \beta_i\beta_j\frac{\partial \phi}{\partial x_i}dv_j dv_i
    factor = betai*betaj*psi*dAi*dAj !*detFi*detFj
    forcei(:) = factor*dot_product( rij,matmul(invFit,sNi))*matmul(invFjt,sNj)
    forcej(:) = factor*dot_product(-rij,matmul(invFjt,sNj))*matmul(invFit,sNi)
    ! distribute forcei, forcej to the nodes of element iel, jel, respectively
    do i = 1, nen
      do j = 1, ndf
        fecontij((i-1)*ndf+j) = fecontij((i-1)*ndf+j) + shpsi(4,i)*forcei(j)
        fecontij((i+nen-1)*ndf+j) = fecontij((i+nen-1)*ndf+j) + shpsj(4,i)*forcej(j)
      end do
    end do
    ! accumulate the contact interaction energy
    energy_cont = energy_cont + betai*betaj*ur*dAi*dAj*dot_product(matmul(invFit,sNi),matmul(invFjt,sNj))
  end do ! end of jgs
!------------------------------------------------------------------------------
! end of the loop over the quadraure points in the element iels
!------------------------------------------------------------------------------
end do ! end of igs

return
end subroutine fconts01

subroutine fcont01(ndim,ndf,nen,nip,icn,jcn, fecontij)

use math
use config, only : vcntp, xr, disp, ngqp
implicit none

! global variables, input
integer :: icn(8), jcn(8)

! global variables, output
real(8) :: fecontij(:)

integer :: ndim, ndf, nen, nip
integer :: ig, jg, i, j, ii, jj
real(8) :: xsj, dvoli, dvolj
real(8) :: sg(4,125), shpsi(4,8), shpsj(4,8)
real(8) :: xli(3,8), uli(3,8), xui(3,8)
real(8) :: xlj(3,8), ulj(3,8), xuj(3,8)
real(8) :: xig(3), xjg(3), rij(3)
real(8) :: r, r0, betai, betaj, epsilon, tmp, dfdr,df2dr2
real(8) :: forceij(3), kii(3,3)

!------------------------------------------------------------------------------
! Set up parameters
!------------------------------------------------------------------------------
! contact parameters
betai = vcntp(1); betaj = vcntp(2)
r0 = vcntp(3); epsilon = vcntp(4)
! integration point in each direction
call int3d(ngqp,nip,sg)

!----------------------------------------------------------------------------
! information related to element iel
!----------------------------------------------------------------------------
! coordinate and displacement
xli = xr(:, icn)
do i = 1, nen
  do j = 1, ndf
    uli(j,i) = disp((icn(i)-1)*ndf + j)
  end do
end do
xui = xli + uli

!----------------------------------------------------------------------------
! information related to element jel
!----------------------------------------------------------------------------
! coordinate and displacement
xlj = xr(:, jcn)
do i = 1, nen
  do j = 1, ndf
    ulj(j,i) = disp((jcn(i)-1)*ndf + j)
  end do
end do
xuj = xlj + ulj

!------------------------------------------------------------------------------
! Loop over the quadraure points in the element iel
!------------------------------------------------------------------------------
do ig = 1, nip
  ! shape functions and their cartesian derivatives wrt the current configuration 
  call shp3d(sg(1,ig),sg(2,ig),sg(3,ig),xui,.false., shpsi,xsj)
  ! volume element associated with ig
  dvoli = xsj*sg(4,ig)
  ! current location of quadrature point ig
  xig = 0.d0
  do i = 1, ndim
    do j = 1, nen
      xig(i) = xig(i) + shpsi(4,j)*xui(i,j) 
    end do
  end do
  ! loop over the quadrature points in jel
  do jg = 1, nip
    ! shape functions and their cartesian derivatives wrt the current configuration 
    call shp3d(sg(1,jg),sg(2,jg),sg(3,jg),xuj,.false., shpsj,xsj)
    ! volume element associated with jg
    dvolj = xsj*sg(4,jg)
    ! current location of quadrature point jg
    xjg = 0.d0
    do i = 1, ndim
      do j = 1, nen
        xjg(i) = xjg(i) + shpsj(4,j)*xuj(i,j) 
      end do
    end do
    ! vector of the two quadrature points
    rij = xig - xjg
    ! distance between the two quadrature points
    r = sqrt(dot_product(rij,rij))
    ! interaction only happens within the cutoff distance
    if(r >= 3.d0*r0) cycle
    ! intermediate variable
    tmp = (r0/r)**6
    ! first derivative of the potential function
    dfdr = 12.d0*epsilon*tmp*(-tmp + 1.d0)/r
    ! force vector \beta_i\beta_j\frac{\partial \phi}{\partial x_i}dv_j dv_i
    forceij(:) = betai*betaj*dfdr*rij(:)*dvoli*dvolj/r
    ! distribute forceij, -forceij to the nodes of element iel, jel, respectively
    do i = 1, nen
      do j = 1, ndf
        fecontij((i-1)*ndf+j) = fecontij((i-1)*ndf+j) + shpsi(4,i)*forceij(j)
        fecontij((i+nen-1)*ndf+j) = fecontij((i+nen-1)*ndf+j) - shpsj(4,i)*forceij(j)
      end do
    end do
    !
  end do ! end of jg
!------------------------------------------------------------------------------
! end of the loop over the quadraure points in the element iel
!------------------------------------------------------------------------------
end do ! end of ig


return
end subroutine fcont01


subroutine fint01(ndim,ndf,nen,nst,icn,d,ul,xl, feint)
!
!----------------------------------------------------------------------------------------
!
!  stif04: Compute the element stiffness matrix and the internal element force vector
!          two-dimensional quad element, Neo-Hookean material model, plain strain case
!   
!  INPUT:
!         ndim               Number of space dimensions (=2)
!         ndf                Number of degrees of freedom per node (=2)
!         nen                Number of node per element (=4)
!         nst                Number of degrees of freedom per element(=8)
!         icn(nen)           Element connnectivity
!         d(:)               Material properties
!         ul(ndf,nen)        Nodal displacements
!         xl(ndim,nen)       Nodal coordinates in the reference configuration
!
!
!  OUTPUT:
!         keint(nst,nst)     Element stiffness
!         feint(nst)         Element internal force vector
!
!  OTHER:
!         sg(1-2,ig)         Isoparameteric coordinate of quadrature point ig
!         sg(3,ig)           weight of quadrature point ig
!
!
!----------------------------------------------------------------------------------------
!
use config, only : ngqp
implicit none

! global variables, input
integer :: ndim, ndf, nen, nst
integer :: icn(nen)
real(8) :: d(:), ul(ndf,nen), xl(ndim,nen)

! global varaibles, output
real(8) :: feint(nst)

! local variables
integer :: i, j, ii, jj
integer :: ig, nip
real(8) :: sg(4,125), shps(4,nen)
real(8) :: xu(ndim,nen)
real(8) :: xsj, dvol
real(8) :: F(ndim,ndim), sig(6), aa(6,6)
real(8) :: bbd(3,6), kij
! real(8) :: bmat(6,24), keint_test(nst, nst)

feint = 0.d0 
!------------------------------------------------------------------------------
! Set up parameters
!------------------------------------------------------------------------------
! integration point in each direction

call int3d(ngqp,nip,sg)
! nodal coordinates in the current configuration
xu = xl + ul
!------------------------------------------------------------------------------
! Loop over the quadrature points in the element
!------------------------------------------------------------------------------
do ig = 1, nip
  ! shape functions and their cartesian derivatives wrt the current configuration 
  call shp3d(sg(1,ig),sg(2,ig),sg(3,ig),xu,.false., shps,xsj)
  ! compute the volume element
  dvol = xsj*sg(4,ig)
  ! get the deformation gradient at ig
  call kine01(ndim,ndf,nen,shps,ul, F)
  ! get the Cauchy stress sig and spatial tagent modulus c at ig
  call modl01(ndim,d,F, sig,aa)
  ! multiply tagent moduli and Cauchy stress by the volume element
  sig = sig*dvol
  aa = aa*dvol
  ! the element internal force vector
  do i = 1, nen
    feint((i-1)*ndf+1) = feint((i-1)*ndf+1) + shps(1,i)*sig(1) + shps(2,i)*sig(4) + shps(3,i)*sig(6)
    feint((i-1)*ndf+2) = feint((i-1)*ndf+2) + shps(2,i)*sig(2) + shps(1,i)*sig(4) + shps(3,i)*sig(5)
    feint((i-1)*ndf+3) = feint((i-1)*ndf+3) + shps(3,i)*sig(3) + shps(2,i)*sig(5) + shps(1,i)*sig(6)
  end do
!------------------------------------------------------------------------------
! End of the loop over the quadraure points in the element
!------------------------------------------------------------------------------
end do

return
end subroutine fint01

end module elem01