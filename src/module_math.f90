module math

use constant
implicit none

contains

function identity_matrix(nn)
   implicit none
   INTEGER, INTENT(IN) :: nn
   REAL(8) , DIMENSION(nn,nn) :: identity_matrix
									! This function returns an n X n array.
   INTEGER :: j            ! local loop index
   
   identity_matrix = 0.0   ! Set each element to zero

   DO j = 1,nn                      ! Change value of each element
      identity_matrix(j,j) = 1.0    ! along diagonal to zero
   END DO
end function identity_matrix

function trace(a) result(b)

implicit none
real(8), intent(in) :: a(:,:)
real(8) :: b

integer :: i

b = 0.d0
do i = 1, ubound(a,1)
  b = b + a(i,i)
end do

end function trace

function norm(a) result(b)

implicit none
real(8), intent(in) :: a(:)
real(8) :: b
b = sqrt(dot_product(a,a))

end function norm

function tensor_product(a,b) result(c)
!--------------------------------------------------------------------
! Purpose: Return the tensor product of two double precision
!          3 dimensional vectors
! by Houfu Fan, 4.18.2013
!--------------------------------------------------------------------
implicit none
real(8), intent(in) :: a(:), b(:)
real(8) :: c(ubound(a,1),ubound(b,1))
integer :: i, j, ma, mb

ma = ubound(a,1); mb = ubound(b,1)
do i = 1, ma
  do j = 1, mb
    c(i,j) = a(i)*b(j)
  end do
end do

end function tensor_product

function cross_product(a,b) result(c)
!--------------------------------------------------------------------
! Purpose: Return the cross product of two double precision
!          3 dimensional vectors
! by Houfu Fan, 4.18.2013
!--------------------------------------------------------------------
implicit none
real(8), intent(in) :: a(3), b(3)
real(8) :: c(3)
integer :: i, j

c(1) = a(2)*b(3) - a(3)*b(2)
c(2) = a(3)*b(1) - a(1)*b(3)
c(3) = a(1)*b(2) - a(2)*b(1)

end function cross_product

function normalize(a) result(b)

implicit none
real(8), intent(in) :: a(:)
real(8) :: b(ubound(a,1))
real(8) :: s

s = sqrt(dot_product(a,a))
b = a/s
 
end function normalize

function inverse(a) result(b)
!
! this subroutine inverts a small square matrix a to b.
!
implicit none
real(8), intent(in) :: a(:,:)
real(8) :: b(ubound(a, 1), ubound(a, 2))
real(8) :: det, j11, j12, j13, j21, j22, j23, j31, j32, j33, con
integer :: ndim, i, k
ndim=ubound(a, 1)
if (ndim == 2)then
  det = a(1, 1) * a(2, 2) - a(1, 2) * a(2, 1)
  b(1, 1) = a(2, 2)
  b(2, 2) = a(1, 1)
  b(1, 2) = -a(1, 2)
  b(2, 1) = -a(2, 1)
  b = b/det
else if (ndim == 3)then
  det = a(1, 1)*(a(2, 2) * a(3, 3) - a(3, 2) * a(2, 3))
  det = det - a(1, 2)*(a(2, 1) * a(3, 3) - a(3, 1) * a(2, 3))
  det = det + a(1, 3)*(a(2, 1) * a(3, 2) - a(3, 1) * a(2, 2))
  j11 = a(2, 2) * a(3, 3) - a(3, 2) * a(2, 3)
  j21 = -a(2, 1) * a(3, 3) + a(3, 1) * a(2, 3)
  j31 = a(2, 1) * a(3, 2) - a(3, 1) * a(2, 2)
  j12 = -a(1, 2) * a(3, 3) + a(3, 2) * a(1, 3)
  j22 = a(1, 1) * a(3, 3) - a(3, 1) * a(1, 3)
  j32 = -a(1, 1) * a(3, 2) + a(3, 1) * a(1, 2)
  j13 = a(1, 2) * a(2, 3) - a(2, 2) * a(1, 3)
  j23 = -a(1, 1) * a(2, 3) + a(2, 1) * a(1, 3)
  j33 = a(1, 1) * a(2, 2) - a(2, 1) * a(1, 2)
  b(1, 1) = j11
  b(1, 2) = j12
  b(1, 3) = j13
  b(2, 1) = j21
  b(2, 2) = j22
  b(2, 3) = j23
  b(3, 1) = j31
  b(3, 2) = j32
  b(3, 3) = j33
  b = b/det
else
  b = a
  do k = 1, ndim
    con = b(k, k)
    b(k, k) = 1.0
    b(k,:) = b(k,:)/con
    do i = 1, ndim
      if (i /= k)then
          con = b(i, k)
          b(i, k) = 0.0
          b(i,:) = b(i,:) - b(k,:) * con
      end if
    end do
  end do
end if

end function inverse

function determinant(a) result(b)
!--------------------------------------------------------------------
! This function returns the determinant
! of a matrix a, 3x3 or 2x2
!--------------------------------------------------------------------

implicit none

real(8), intent(in) :: a(:,:)
real(8) :: b
integer :: ndim

ndim = ubound(a, 1)
if (ndim == 2) then
    b = a(1, 1) * a(2, 2) - a(1, 2) * a(2, 1)
elseif (ndim == 3) then
    b = a(1, 1)*(a(2, 2) * a(3, 3) - a(3, 2) * a(2, 3))
    b = b - a(1, 2)*(a(2, 1) * a(3, 3) - a(3, 1) * a(2, 3))
    b = b + a(1, 3)*(a(2, 1) * a(3, 2) - a(3, 1) * a(2, 2))
else
    write(*, *) 'Warning, matrix dimension not acceptable in determinant, exit!'
    stop
end if

end function determinant


function add(a, b) result(c)

implicit none

real(8), intent(in) :: a(:,:), b(:,:)
real(8) :: c(ubound(a, 1), ubound(a, 2))
integer :: ia, ja, ib, jb, i, j

ia = ubound(a, 1); ja = ubound(a, 2)
ib = ubound(b, 1); jb = ubound(b, 2)

if (ia /= ib .or. ja /= jb) then
    write(*, *) 'Matrix dimension not conform in add, exit!'
    stop
end if

do i =  1, ia
    do j = 1, ja
        c(i, j) = a(i, j) + b(i, j)
    end do
end do

end function add

function matrixmul(a, b) result(c)

implicit none

real(8), intent(in) :: a(:,:), b(:,:)
real(8) :: c(ubound(a, 1), ubound(a, 2))
integer :: ia, ja, ib, jb, i, j, k

ia = ubound(a, 1); ja = ubound(a, 2)
ib = ubound(b, 1); jb = ubound(b, 2)

if (ja /= ib) then
    write(*, *) 'Matrix dimension not conform in matrixmul, exit!'
    stop
end if

c = 0.d0
do i = 1, ia
  do j = 1, jb
    do k = 1, ja
      c(i,j) = c(i,j) + a(i,k)*b(k,j)
    end do
  end do
end do

end function matrixmul

subroutine invert(a,nmax,ndm)
!--------------------------------------------------------------------
! Purpose: Invert small square matrix
! 
! Inputs:
!    a(ndm,*) - Matrix to be inverted
!    nmax     - Size of upper submatrix to invert
!    ndm      - Dimension of array
! 
! Outputs:
!    a(ndm,*) - Submatrix replaces original terms, others not
!               changed
!--------------------------------------------------------------------
implicit  none

integer ::   i,j,n,ndm,nmax
real(8) ::   d, a(:,:)

do n = 1,nmax
  if(a(n,n).ne.0.0d0) then
    d = 1.d0/a(n,n)
    do j = 1,nmax
      a(n,j) = -a(n,j)*d
    end do

    do i = 1,nmax
      if(n.ne.i) then
        do j = 1,nmax
          if(n.ne.j) a(i,j) = a(i,j) + a(i,n)*a(n,j)
        end do
      endif
      a(i,n) = a(i,n)*d
    end do
    a(n,n) = d
  else
    write(*,*) ' *WARNING* Zero pivot in INVERT'
    stop
  endif
end do

return
end subroutine invert

subroutine Get_Surface_Operator(sn, PP)

implicit none

integer :: ii
real(8) :: sn(3), PP(3,3)

PP(:,:) = -tensor_product(sn,sn)
do ii = 1, 3
  PP(ii,ii) = 1.d0 + PP(ii,ii)
end do 

return
end subroutine Get_Surface_Operator

subroutine conjgrad(a,b,x)

  implicit none

  real(8), intent(in) :: a(:,:), b(:)
  real(8), intent(inout) :: x(:)

  integer :: i, imax
  real(8) :: epsilon, alpha, beta, delta0, delta_old, delta_new
  real(8) :: r(size(b)), d(size(b)), q(size(b))

  imax = 10000
  epsilon = 1.0d-5
  r = b - matmul(a,x)
  d = r
  delta_new = dot_product(r,r)
  delta0 = delta_new

  do i = 1, imax
     if(delta_new <= epsilon * epsilon * delta0) exit
     q = matmul(a,d)
     alpha = delta_new/dot_product(d,q)
     x = x + alpha * d
     if(mod(i,50) == 0) then
        r = b - matmul(a,x)
     else
        r = r - alpha*q
     end if
     delta_old = delta_new
     delta_new = dot_product(r,r)
     beta = delta_new/delta_old
     d = r + beta * d
  end do

  write(*,*) '# of interation needed in the conjgate gradient:', i-1
  if(i == imax + 1) then 
     write(*,*) 'not covergent after', imax, 'iterations in cg'
     write(*,*) 'ratio =', delta_new/delta0 ,'quit!'
     stop
  end if

  return
end subroutine conjgrad

!----------------------------------------------------------------------------------------
!
!                cholesky decomposition : choldc
!                and its back substitution : cholsl
!              
!----------------------------------------------------------------------------------------
subroutine choldc(a,p)

implicit none

real(8), intent(inout) :: a(:,:)
real(8), intent(out) :: p(:)
integer :: i,n
real(8) :: summ
n=size(p)
do i=1,n
	summ=a(i,i)-dot_product(a(i,1:i-1),a(i,1:i-1))
	if (summ <= 0.0) then
	  write(*,*) 'choldc failed'
	  stop
	endif
	p(i)=sqrt(summ)
	a(i+1:n,i)=(a(i,i+1:n)-matmul(a(i+1:n,1:i-1),a(i,1:i-1)))/p(i)
end do

end subroutine choldc

subroutine cholsl(a,p,b,x)

implicit none
real(8), intent(in) :: a(:,:)
real(8), intent(in) :: p(:),b(:)
real(8), intent(inout) :: x(:)
integer :: i,n
n=size(p)
do i=1,n
	x(i)=(b(i)-dot_product(a(i,1:i-1),x(1:i-1)))/p(i)
end do
do i=n,1,-1
	x(i)=(x(i)-dot_product(a(i+1:n,i),x(i+1:n)))/p(i)
end do

end subroutine cholsl

!----------------------------------------------------------------------------------------
!
!                LU decomposition : ludcmp
!                and its back substitution : lubksb
!              
!----------------------------------------------------------------------------------------
subroutine ludcmp(a,indx,d)

implicit none
real(8), intent(inout) :: a(:,:)
integer, intent(out) :: indx(:) ! store the pivoting
real(8), intent(out) :: d

real(8) :: vv(size(a,1)), tmp(size(a,2))
real(8), parameter :: tiny=1.0d-20
integer :: j,n,imax

n=size(indx) 
d=1.d0
vv=maxval(abs(a),dim=2)
if (any(vv == 0.d0)) then
  write(*,*) 'singular matrix in ludcmp'
  stop
endif
vv=1.d0/vv
do j=1,n
	imax=(j-1)+imaxloc(vv(j:n)*abs(a(j:n,j)))
	if (j /= imax) then
		tmp(:) = a(imax,:)
		a(imax,:) = a(j,:)
		a(j,:) = tmp(:)
		d=-d
		vv(imax)=vv(j)
	end if
	indx(j)=imax
	if (a(j,j) == 0.d0) a(j,j)=tiny
	a(j+1:n,j)=a(j+1:n,j)/a(j,j)
	a(j+1:n,j+1:n)=a(j+1:n,j+1:n)-tensor_product(a(j+1:n,j),a(j,j+1:n))
end do
end subroutine ludcmp

function imaxloc(arr)
real(8), dimension(:), intent(in) :: arr
integer :: imaxloc
integer, dimension(1) :: imax
imax=maxloc(arr(:))
imaxloc=imax(1)
end function imaxloc

subroutine lubksb(a,indx,b)

implicit none
real(8), intent(in) :: a(:,:)
integer, intent(in) :: indx(:)
real(8), intent(inout) :: b(:)
integer :: i,n,ii,ll
real(8) :: summ
n=size(indx)
ii=0
do i=1,n
	ll=indx(i)
	summ=b(ll)
	b(ll)=b(i)
	if (ii /= 0) then
		summ=summ-dot_product(a(i,ii:i-1),b(ii:i-1))
	else if (summ /= 0.0) then
		ii=i
	end if
	b(i)=summ
end do
do i=n,1,-1
	b(i) = (b(i)-dot_product(a(i,i+1:n),b(i+1:n)))/a(i,i)
end do
end subroutine lubksb


end module math