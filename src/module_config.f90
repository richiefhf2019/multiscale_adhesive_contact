module config

implicit none

! file channels for input & output
integer, parameter :: ior = 10
integer, parameter :: iow = 11
integer, parameter :: ivs = 12 ! tecplot file channel
integer, parameter :: ich = 13
integer, parameter :: ict = 14

! output id for tecplot
integer :: ioutput, iprint

! integration method, 1--body-body, 2--surface-surface
integer, save :: intmet 

! master
integer, save :: numnp, numel
integer, save :: ndim, nen, ndf, nip, nst, neq, nequ ! nequ--number of equations, unknown

! time
real(8), save :: t, dt

! number of Gauss Quadrature Point in each direction
integer, save :: ngqp, ngqps

! mesh
real(8), allocatable, save :: xr(:,:)
integer, allocatable, save :: cn(:,:), etype(:), ntype(:)

! variables related to surface approach
integer, save :: numels ! number of surface elements
integer, save :: nips ! number of surface quadrature points
integer, allocatable, save :: stob(:), srfid(:)
real(8), allocatable, save :: dA(:,:)        ! surface element area, associated with a quadrature point
real(8), allocatable, save :: sN(:,:,:)      ! unit out normal in the reference configuration
real(8), allocatable, save :: dNdX(:,:,:,:)  ! Gradient of the reference unit out normal

! material parameters
integer, save :: nmatp
real(8), allocatable, save :: vmatp(:)

! contact parameters
integer, save :: ncntp
real(8), allocatable, save :: vcntp(:)

! energy conservation
real(8), save :: energy_kin, energy_int, energy_cont

! stiffness matrices and force vectors
real(8), allocatable, save :: ktot(:,:)
real(8), allocatable, save :: ftot(:)
real(8), allocatable, save :: fext(:)

! solution variables
real(8), allocatable, save :: disp(:), residual(:), dud(:)
integer, allocatable, save :: indx(:) ! for LU decomposition

! for dynamic analysis
real(8), allocatable, save :: vel(:), vel_half(:), acc(:), mass(:)

! local element variables
integer, allocatable :: icn(:), jcn(:)
real(8), allocatable :: xl(:,:), ul(:,:), xu(:,:)
real(8), allocatable :: feint(:), keint(:,:)
real(8), allocatable :: fecontij(:), kecontij(:,:)

! boundary condition --  essential boundary condition
integer, save :: ndispb
integer, allocatable, save :: bid(:)    ! global node index
integer, allocatable, save :: bdir(:)   ! displacement direction
real(8), allocatable, save :: bvl(:)    ! dispalcement value

! boundary condition -- force boundary condition
integer, save :: nforceb
integer, allocatable, save :: bfid(:)    ! global node index
integer, allocatable, save :: bfdir(:)   ! displacement direction
real(8), allocatable, save :: bfvl(:)    ! dispalcement value
integer, save :: nload ! number of load steps

! static condensation
integer, allocatable, save :: idd(:), idf(:)
real(8), allocatable, save :: kff(:,:), rf(:)

! lumped projection matrix, projected stresses
real(8), allocatable, save :: lpm(:), stp(:,:)

! convergence tolerance, for r*du
real(8), save :: tolerance

! working array for pbicg
real(8), allocatable, save :: dxf(:), dxfr(:)
real(8), save, allocatable :: wrk(:,:)

contains

subroutine allocate_variables

implicit none

integer :: isok

100 FORMAT('Error in Memory Allocation for Array',1x,A, 'in allocate_variables')

! global force vector, displacement, stiffness matrix
allocate( ftot(neq), STAT = isok); if(isok /= 0) then; write(*,100) 'ftot'; stop; endif
allocate( fext(neq), STAT = isok); if(isok /= 0) then; write(*,100) 'fext'; stop; endif
!allocate( ktot(neq, neq), STAT = isok); if(isok /= 0) then; write(*,100) 'ktot'; stop; endif
allocate( residual(neq), STAT = isok); if(isok /= 0) then; write(*,100) 'residual'; stop; endif
allocate( disp(neq), STAT = isok); if(isok /= 0) then; write(*,100) 'disp'; stop; endif

ftot = 0.d0
!ktot = 0.d0
disp = 0.d0;  residual = 0.d0

! for dynamic analysis
allocate(vel(neq), acc(neq), vel_half(neq))
allocate(mass(numnp))

vel = 0.d0; vel_half = 0.d0; acc = 0.d0; mass = 0.d0



! element force vector, displacement, stiffness matrix
allocate( icn(nen), STAT = isok); if(isok /= 0) then; write(*,100) 'icn'; stop; endif
allocate( jcn(nen), STAT = isok); if(isok /= 0) then; write(*,100) 'jcn'; stop; endif
allocate( xl(ndim, nen), STAT = isok); if(isok /= 0) then; write(*,100) 'xl'; stop; endif
allocate( ul(ndf, nen), STAT = isok); if(isok /= 0) then; write(*,100) 'ul'; stop; endif
allocate( xu(ndim, nen), STAT = isok); if(isok /= 0) then; write(*,100) 'xu'; stop; endif
allocate( feint(nst), STAT = isok); if(isok /= 0) then; write(*,100) 'feint'; stop; endif
allocate( keint(nst, nst), STAT = isok); if(isok /= 0) then; write(*,100) 'keint'; stop; endif
allocate( fecontij(2*nst), STAT = isok); if(isok /= 0) then; write(*,100) 'fecontij'; stop; endif
allocate( kecontij(2*nst, 2*nst), STAT = isok); if(isok /= 0) then; write(*,100) 'kcontij'; stop; endif

icn = 0; jcn = 0; xl = 0.d0; ul = 0.d0; xu = 0.d0
feint = 0.d0; keint = 0.d0
fecontij = 0.d0; kecontij = 0.d0

! indices for essential boundary condition, static condenstation
allocate( idd(ndispb), STAT = isok); if(isok /= 0) then; write(*,100) 'idd'; stop; endif
allocate( dud(ndispb), STAT = isok); if(isok /= 0) then; write(*,100) 'dud'; stop; endif
allocate( idf(nequ), STAT = isok); if(isok /= 0) then; write(*,100) 'idf'; stop; endif
!allocate( kff(nequ, nequ), STAT = isok); if(isok /= 0) then; write(*,100) 'kff'; stop; endif
!allocate( rf(nequ), STAT = isok); if(isok /= 0) then; write(*,100) 'rf'; stop; endif
!allocate( indx(nequ), STAT = isok); if(isok /= 0) then; write(*,100) 'indx'; stop; endif

! lumped mass matrix and projected stress
allocate(lpm(numnp),stp(2,numnp))
lpm = 0.d0; stp = 0.d0

!!
!allocate( wrk(nequ,12), STAT = isok); if(isok /= 0) then; write(*,100) 'wrk'; stop; endif
!allocate( dxf(nequ), STAT = isok); if(isok /= 0) then; write(*,100) 'dxf'; stop; endif
!allocate( dxfr(nequ), STAT = isok); if(isok /= 0) then; write(*,100) 'dxfr'; stop; endif
!
!idd = 0; idf = 0
!kff = 0.d0
!rf = 0.d0; indx = 0
!dxf = 0.d0; dxfr = 0.d0; wrk = 0.d0


! matrices related to the surface elements
allocate(dA(nips,numels))        ! surface element area, associated with a quadrature point
allocate(sN(ndim,nips,numels))      ! unit out normal in the reference configuration
allocate(dNdX(ndim,ndim,nips,numels))  ! Gradient of the reference unit out normal

dA = 0.d0; sN = 0.d0; dNdX = 0.d0


return
end subroutine allocate_variables

end module config