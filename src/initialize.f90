subroutine initialize

use config
implicit none

integer :: i, j, k
integer :: isok
real(8) :: youngs, nu
real(8) :: ratex, ratey, ratez
character(200) :: aline

open(unit = ior, file = 'cgcm2d-dynamic.inp', status = 'unknown', action = 'read')
open(unit = iow, file = 'cgcm2d-dynamic.out', status = 'unknown', action = 'write')
open(unit = ivs, file = 'dynamic-cont.dat', status = 'unknown', action = 'write')

100 FORMAT('Error in Memory Allocation for Array',1x,A)

! integration method
call readaline(ior,iow,aline)
read(aline,*) intmet

! number of Gauss Quadrature Points in each direction
! ngqp, ngqps
call readaline(ior,iow,aline)
read(aline,*) ngqp, ngqps

! master
call readaline(ior,iow,aline)
read(aline,*) ndim, nen, ndf, nip;  nst = ndf*nen

! mesh
call readaline(ior,iow,aline)
read(aline,*) numnp, numel; neq = ndf*numnp

! allocate arrays for mesh
allocate( xr(ndim,numnp), STAT = isok); if(isok /= 0) then; write(*,100) 'xr'; stop; endif
allocate( cn(nen,numel), STAT = isok); if(isok /= 0) then; write(*,100) 'cn'; stop; endif
allocate( etype(numel), STAT = isok); if(isok /= 0) then; write(*,100) 'etype'; stop; endif
allocate( ntype(numnp), STAT = isok); if(isok /= 0) then; write(*,100) 'ntype'; stop; endif

! t, dt, nload, iprint
call readaline(ior,iow,aline)
read(aline,*) t, dt, nload, iprint


! number of parameters for material properties
call readaline(ior,iow,aline)
read(aline,*) nmatp
! allocate material properties array
allocate( vmatp(nmatp), STAT = isok); if(isok /= 0) then; write(*,100) 'vmatp'; stop; endif
! read material properties
call readaline(ior,iow,aline)
read(aline,*) (vmatp(i),i = 1, nmatp)

! number of contact parameters
call readaline(ior,iow,aline)
read(aline,*) ncntp
! allocate contact parameter array
allocate( vcntp(ncntp), STAT = isok); if(isok /= 0) then; write(*,100) 'vcntp'; stop; endif
! read contact parameters
call readaline(ior,iow,aline)
read(aline,*) (vcntp(i), i = 1, ncntp)

! convergence tolerance
call readaline(ior,iow,aline)
read(aline,*) tolerance

! coordinate rates
call readaline(ior,iow,aline)
read(aline,*) ratex, ratey

! node coordinates
do i = 1, numnp
  call readaline(ior,iow,aline)
  read(aline,*) j, (xr(k,i), k = 1, ndim)
  xr(1,i) = xr(1,i)*ratex
  xr(2,i) = xr(2,i)*ratey
end do

! mesh connectivity
do i = 1, numel
  call readaline(ior,iow,aline)
  read(aline,*) j, (cn(k,i), k = 1, nen), etype(i)
  ntype(cn(:,i)) = etype(i)
end do

! surface element information
call readaline(ior, iow, aline)
read(aline,*) numels, nips
if(numels > 0) then
  ! allocate the matrix
  allocate(stob(numels),srfid(numels))
  do i = 1, numels
    call readaline(ior, iow, aline)
    read(aline,*) j, stob(i), srfid(i)
  end do
endif

! manipulate the initial gap
do i = 1, numnp
  !if(ntype(i) == 1) then
  !  xr(1,i) = xr(1,i) + 1.d0 
  !endif

  !if(ntype(i) == 2) then
  !  xr(2,i) = xr(2,i) - 30.d0 
  !endif
end do

! number of displacement boundary
call readaline(ior,iow,aline)
read(aline,*) ndispb

! number of unknown variables
nequ = neq-ndispb

! allocate displacement boundary arrays
allocate( bid(ndispb), STAT = isok); if(isok /= 0) then; write(*,100) 'bid'; stop; endif
allocate( bdir(ndispb), STAT = isok); if(isok /= 0) then; write(*,100) 'bdir'; stop; endif
allocate( bvl(ndispb), STAT = isok); if(isok /= 0) then; write(*,100) 'bvl'; stop; endif

! read displacement boundary 
do i = 1, ndispb
  call readaline(ior,iow,aline)
  read(aline,*) j, bid(i), bdir(i), bvl(i)
end do 


! number of force boundary
call readaline(ior,iow,aline)
read(aline,*) nforceb

if(nforceb > 0) then
  ! allocate force boundary arrays
  allocate( bfid(nforceb), STAT = isok); if(isok /= 0) then; write(*,100) 'bfid'; stop; endif
  allocate( bfdir(nforceb), STAT = isok); if(isok /= 0) then; write(*,100) 'bfdir'; stop; endif
  allocate( bfvl(nforceb), STAT = isok); if(isok /= 0) then; write(*,100) 'bfvl'; stop; endif
  
  ! read force boundary
  do i = 1, nforceb
    call readaline(ior,iow,aline)
    read(aline,*) j, bfid(i), bfdir(i), bfvl(i)
  end do
endif

! allocate global variables
call allocate_variables


fext = 0.d0
! external force
do i = 1, nforceb
  j = (bfid(i)-1)*ndf + bfdir(i) ! global index of the dof i
  fext(j) = bfvl(i) ! value of the driven displacement
end do

! preprocess the essential boundary condition
do i = 1, ndispb
  j = (bid(i)-1)*ndf + bdir(i) ! global index of the dof i
  idd(i) = j ! indices for the driven displacement
  dud(i) = bvl(i) ! value of the driven displacement
end do

isok = 1; k = 0
do i = 1, neq
  do j = 1, ndispb
    if(i /= idd(j)) then
      cycle
    else 
      isok = 0
      exit
    endif
  end do
  if(isok == 1) then
    k = k + 1
    idf(k) = i ! find the indices for the unkown dispacement
  endif
  isok = 1
end do

close(ior)
close(iow)

return
end subroutine initialize