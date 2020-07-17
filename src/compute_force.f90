subroutine compute_force
! purpose: compute the force vector

use elem04
use elem02
use fem
use config
implicit none

integer :: iel, jel, iels, jels, i, j
real(8) :: vd(3)
real(8) :: factor = 10.d0

real(8) :: start, finish

vd = vmatp; vd(1) = vmatp(1)*1.d5
!----------------------------------------------------------------------------
! loop over the volume elements to compute fint, using d^{(i+1)}
!----------------------------------------------------------------------------
ftot = 0.d0
do iel = 1, numel
  ! connectivity
  icn = cn(:,iel)
  ! coordinate and displacement
  xl = xr(:, icn) 
  do i = 1, nen
    do j = 1, ndf
      ul(j,i) = disp((icn(i)-1)*ndf + j)
    end do
  end do
  ! initialize the internal force vector and matrix for iel
  feint = 0.d0
  ! get the element force vector and stiffness matrix
  !if(etype(iel) == 1) then
  !  call fint04(ndim,ndf,nen,nst,icn,vd,ul,xl, feint)
  !else
  !  call fint04(ndim,ndf,nen,nst,icn,vmatp,ul,xl, feint)
  !endif
  call fint04(ndim,ndf,nen,nst,icn,vmatp,ul,xl, feint)
  ! assemble feint into the global force vector fint
  call assemble_vector(ndf,nen,numnp,icn,feint, ftot)
!----------------------------------------------------------------------------
! end of the loop over the volume elements to compute fint and kint
!----------------------------------------------------------------------------
end do ! end of the loop over iel
! compute the internal energy
energy_int = 0.d0
do i = 1, numnp
  do j = 1, ndf
    energy_int = energy_int + 0.5d0*ftot((i-1)*ndf+j)*disp((i-1)*ndf+j)
  end do
end do
!----------------------------------------------------------------------------
! loop over the volume element pairs to compute contact force and stiffness
!----------------------------------------------------------------------------
if (intmet == 2) then
  !call cpu_time(start)
  energy_cont = 0.d0
  do iels = 1, numels-1
    ! bulk element index, corresponding to surface element iels
    iel = stob(iels); icn = cn(:,iel)
    do jels = iels+1, numels ! don't over count the interaction force
      ! bulk element index, corresponding to surface element jels
      jel = stob(jels); jcn = cn(:,jel)
      ! contact forces come from inter-body interaction
      if(etype(jel) == etype(iel)) cycle
      ! initialize the contact force vector and matrix between element iel and element jel
      fecontij = 0.d0;
      ! get the element-element contact force vector and stiffness matrix
      call fconts04(ndim,ndf,nen,nip,icn,jcn,iels,jels, fecontij)
      ! assemble fecontij into the global force vector fcont
      call assemble_vectorij(ndf,nen,numnp,icn,jcn,fecontij, ftot)
    end do ! end of jels
  !----------------------------------------------------------------------------
  ! end of the loop to compute contact force and stiffness
  !----------------------------------------------------------------------------
  end do
  !call cpu_time(finish)
   ! write(*,*) 'surface integral', finish-start
elseif (intmet == 1) then
  !call cpu_time(start)
  do iel = 1, numel-1
    ! bulk element index, corresponding to surface element iel
    icn = cn(:,iel)
    do jel = iel+1, numel ! don't over count the interaction force
      ! bulk element index, corresponding to surface element jel
      jcn = cn(:,jel)
      ! contact forces come from inter-body interaction
      if(etype(jel) == etype(iel)) cycle
      ! initialize the contact force vector and matrix between element iel and element jel
      fecontij = 0.d0
      ! get the element-element contact force vector and stiffness matrix
      call fcont04(ndim,ndf,nen,nip,icn,jcn, fecontij)
      ! assemble fecontij into the global force vector fcont
      call assemble_vectorij(ndf,nen,numnp,icn,jcn,fecontij, ftot)
    end do ! end of jel
  !----------------------------------------------------------------------------
  ! end of the loop to compute contact force and stiffness
  !----------------------------------------------------------------------------
  end do
  !call cpu_time(finish)
  !write(*,*) 'volume integral', finish-start
endif

!stop

return
end subroutine compute_force