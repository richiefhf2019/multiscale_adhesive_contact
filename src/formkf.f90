subroutine formkf
! purpose: form the total stiffness matrix and force vector

use elem04
use elem02
use fem
use config
implicit none

integer :: iel, jel, iels, jels, i, j
!----------------------------------------------------------------------------
! loop over the volume elements to compute fint and kint, using d^{(i+1)}
!----------------------------------------------------------------------------
ftot = 0.d0; ktot = 0.d0
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
  feint = 0.d0; keint = 0.d0
  ! get the element force vector and stiffness matrix
  call stif04(ndim,ndf,nen,nst,icn,vmatp,ul,xl, feint,keint)
  ! assemble feint into the global force vector fint
  call assemble_vector(ndf,nen,numnp,icn,feint, ftot)
  ! assemble keint into the global stiffness matrix kint
  call assemble_matrix(ndf,nen,numnp,icn,keint, ktot)
!----------------------------------------------------------------------------
! end of the loop over the volume elements to compute fint and kint
!----------------------------------------------------------------------------
end do ! end of the loop over iel
!----------------------------------------------------------------------------
! loop over the volume element pairs to compute contact force and stiffness
!----------------------------------------------------------------------------
if (intmet == 2) then
  do iels = 1, numels-1
    ! bulk element index, corresponding to surface element iels
    iel = stob(iels); icn = cn(:,iel)
    do jels = iels+1, numels ! don't over count the interaction force
      ! bulk element index, corresponding to surface element jels
      jel = stob(jels); jcn = cn(:,jel)
      ! contact forces come from inter-body interaction
      if(etype(jel) == etype(iel)) cycle
      ! initialize the contact force vector and matrix between element iel and element jel
      fecontij = 0.d0; kecontij = 0.d0
      ! get the element-element contact force vector and stiffness matrix
      call stifconts04(ndim,ndf,nen,nip,icn,jcn,iels,jels, fecontij,kecontij)
      ! assemble fecontij into the global force vector fcont
      call assemble_vectorij(ndf,nen,numnp,icn,jcn,fecontij, ftot)
      ! assemble kecontij into the global stiffness matrix kcont
      call assemble_matrixij(ndf,nen,numnp,icn,jcn,kecontij, ktot)
    end do ! end of jels
  !----------------------------------------------------------------------------
  ! end of the loop to compute contact force and stiffness
  !----------------------------------------------------------------------------
  end do
elseif (intmet == 1) then
  do iel = 1, numel-1
    ! bulk element index, corresponding to surface element iel
    icn = cn(:,iel)
    do jel = iel+1, numel ! don't over count the interaction force
      ! bulk element index, corresponding to surface element jel
      jcn = cn(:,jel)
      ! contact forces come from inter-body interaction
      if(etype(jel) == etype(iel)) cycle
      ! initialize the contact force vector and matrix between element iel and element jel
      fecontij = 0.d0; kecontij = 0.d0
      ! get the element-element contact force vector and stiffness matrix
      call stifcont02(ndim,ndf,nen,nip,icn,jcn, fecontij,kecontij)
      ! assemble fecontij into the global force vector fcont
      call assemble_vectorij(ndf,nen,numnp,icn,jcn,fecontij, ftot)
      ! assemble kecontij into the global stiffness matrix kcont
      call assemble_matrixij(ndf,nen,numnp,icn,jcn,kecontij, ktot)
    end do ! end of jel
  !----------------------------------------------------------------------------
  ! end of the loop to compute contact force and stiffness
  !----------------------------------------------------------------------------
  end do
endif

return
end subroutine formkf