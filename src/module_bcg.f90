module bcg

! basical type definition of the sparse matrix
type sprs2
	integer :: n,len
	real(8), dimension(:), pointer :: val
	integer, dimension(:), pointer :: irow
	integer, dimension(:), pointer :: jcol
end type sprs2

contains

subroutine linbcg(sa,b,x,itol,tol,itmax,iter,err)
!
!----------------------------------------------------------------------------------------
!
!  biconjugate gradient solution of sparse systems
!
!----------------------------------------------------------------------------------------
!
implicit none
type(sprs2), intent(in) :: sa
real(8), dimension(:), intent(in) :: b
real(8), dimension(:), intent(inout) :: x
integer, intent(in) :: itol,itmax
real(8), intent(in) :: tol
integer, intent(out) :: iter
real(8), intent(out) :: err
real(8), parameter :: eps=1.0e-14
integer :: n
real(8) :: ak,akden,bk,bkden,bknum,bnrm,dxnrm,xnrm,zm1nrm,znrm
real(8), dimension(size(b)) :: p,pp,r,rr,z,zz

n=size(x)
iter=0
call atimes(sa,x,r,0)
r=b-r
rr=r
call atimes(sa,r,rr,0)
select case(itol)
	case(1)
		bnrm=snrm(b,itol)
		call asolve(sa,r,z,0)
	case(2)
		call asolve(sa,b,z,0)
		bnrm=snrm(z,itol)
		call asolve(sa,r,z,0)
	case(3:4)
		call asolve(sa,b,z,0)
		bnrm=snrm(z,itol)
		call asolve(sa,r,z,0)
		znrm=snrm(z,itol)
	case default
		write(*,*) 'illegal itol in linbcg'
        stop
end select
do
	if (iter > itmax) exit
	iter=iter+1
	call asolve(sa,rr,zz,1)
	bknum=dot_product(z,rr)
	if (iter == 1) then
		p=z
		pp=zz
	else
		bk=bknum/bkden
		p=bk*p+z
		pp=bk*pp+zz
	end if
	bkden=bknum
	call atimes(sa,p,z,0)
	akden=dot_product(z,pp)
	ak=bknum/akden
	call atimes(sa,pp,zz,1)
	x=x+ak*p
	r=r-ak*z
	rr=rr-ak*zz
	call asolve(sa,r,z,0)
	select case(itol)
		case(1)
			err=snrm(r,itol)/bnrm
		case(2)
			err=snrm(z,itol)/bnrm
		case(3:4)
			zm1nrm=znrm
			znrm=snrm(z,itol)
			if (abs(zm1nrm-znrm) > eps*znrm) then
				dxnrm=abs(ak)*snrm(p,itol)
				err=znrm/abs(zm1nrm-znrm)*dxnrm
			else
				err=znrm/bnrm
				cycle
			end if
			xnrm=snrm(x,itol)
			if (err <= 0.5*xnrm) then
				err=err/xnrm
			else
				err=znrm/bnrm
				cycle
			end if
	end select
	if (err <= tol) then
      write (*,*) 'number of biconjugate iteractons =',iter,' err=',err
      exit
    endif
end do
return
end subroutine linbcg

subroutine atimes(sa,x,r,itrnsp)
!
! used by linbcg for sparse multiplication
!
type(sprs2), intent(in) :: sa
real(8), dimension(:), intent(in) :: x
real(8), dimension(:), intent(out) :: r
integer, intent(in) :: itrnsp
integer :: n
n=size(r)
if (itrnsp == 0) then
	call sprsax(sa,x,r)
else
	call sprstx(sa,x,r)
end if
return
end subroutine atimes

subroutine asolve(sa,b,x,itrnsp)
!
! used by linbcg for preconditioner
!
type(sprs2), intent(in) :: sa
real(8), dimension(:), intent(in) :: b
real(8), dimension(:), intent(out) :: x
integer, intent(in) :: itrnsp
integer :: ndum
ndum=size(x)
call sprsdiag(sa,x)
if (any(x == 0.0)) then
  write(*,*) 'asolve: singular diagonal matrix'
  stop
endif
x=b/x
return
end subroutine asolve

function snrm(sx,itol)
!
! used by linbcg for vector norm
!
implicit none
real(8), dimension(:), intent(in) :: sx
integer, intent(in) :: itol
real(8) :: snrm
if (itol <= 3) then
	snrm=sqrt(dot_product(sx,sx))
else
	snrm=maxval(abs(sx))
end if
end function snrm

subroutine sprsax(sa,x,b)
!
! product of sparse matrix and vector, used by atimes
!
implicit none
type(sprs2), intent(in) :: sa
real(8), dimension (:), intent(in) :: x
real(8), dimension (:), intent(out) :: b
integer :: ndum
ndum=size(x)
b=0.d0
call scatter_add(b,sa%val*x(sa%jcol),sa%irow)
return
end subroutine sprsax

subroutine sprstx(sa,x,b)
!
! product of transpose sparse matrix and vector, used by atimes
!
implicit none
type(sprs2), intent(in) :: sa
real(8), dimension (:), intent(in) :: x
real(8), dimension (:), intent(out) :: b
integer :: ndum
ndum=size(x)
b=0.d0
call scatter_add(b,sa%val*x(sa%irow),sa%jcol)
return
end subroutine sprstx

subroutine scatter_add(dest,source,dest_index)
real(8), dimension(:), intent(out) :: dest
real(8), dimension(:), intent(in) :: source
integer, dimension(:), intent(in) :: dest_index
integer :: m,n,j,i
n=size(source)
m=size(dest)
do j=1,n
	i=dest_index(j)
	if (i > 0 .and. i <= m) dest(i)=dest(i)+source(j)
end do
return
end subroutine scatter_add

subroutine sprsdiag(sa,b)
implicit none
type(sprs2), intent(in) :: sa
real(8), dimension(:), intent(out) :: b
real(8), dimension(size(b)) :: val
integer :: k,l,ndum,nerr
integer, dimension(size(b)) :: i
logical, dimension(:), allocatable :: mask
ndum=size(b)
l=sa%len
allocate(mask(l))
mask = (sa%irow(1:l) == sa%jcol(1:l))
call array_copy(pack(sa%val(1:l),mask),val,k,nerr)
i(1:k)=pack(sa%irow(1:l),mask)
deallocate(mask)
b=0.0
b(i(1:k))=val(1:k)
return
end subroutine sprsdiag

subroutine array_copy(src,dest,n_copied,n_not_copied)
real(8), dimension(:), intent(in) :: src
real(8), dimension(:), intent(out) :: dest
integer, intent(out) :: n_copied, n_not_copied
n_copied=min(size(src),size(dest))
n_not_copied=size(src)-n_copied
dest(1:n_copied)=src(1:n_copied)
return
end subroutine array_copy
	
subroutine sprsin(a,thresh,sa)
! convert matrix to sparse format
implicit none
real(8), dimension(:,:), intent(in) :: a
real(8), intent(in) :: thresh
type(sprs2), intent(out) :: sa
integer :: n,len
logical, dimension(size(a,1),size(a,2)) :: mask
n=size(a,1)
mask=abs(a)>thresh
len=count(mask)
allocate(sa%val(len),sa%irow(len),sa%jcol(len))
sa%n=n
sa%len=len
sa%val=pack(a,mask)
sa%irow=pack(spread(arth(1,1,n),2,n),mask)
sa%jcol=pack(spread(arth(1,1,n),1,n),mask)
return
end subroutine sprsin

function arth(first,increment,n)
integer, parameter :: npar_arth=16,npar2_arth=8
integer, intent(in) :: first,increment,n
integer, dimension(n) :: arth
integer :: k,k2,temp
if (n > 0) arth(1)=first
if (n <= npar_arth) then
	do k=2,n
		arth(k)=arth(k-1)+increment
	end do
else
	do k=2,npar2_arth
		arth(k)=arth(k-1)+increment
	end do
	temp=increment*npar2_arth
	k=npar2_arth
	do
		if (k >= n) exit
		k2=k+k
		arth(k+1:min(k2,n))=temp+arth(1:min(k,n-k))
		temp=temp+temp
		k=k2
	end do
end if
end function arth

end module bcg