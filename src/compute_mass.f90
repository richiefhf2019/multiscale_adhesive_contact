subroutine compute_mass

use config
use elem04
implicit none
integer :: iel
real(8), allocatable :: mm(:)

allocate(mm(nen))

mass = 0.d0
do iel = 1, numel
  icn = cn(:,iel)
  xl = xr(:,icn)
  mm = 0.d0
  call mass04(ndim,nen,vmatp,xl, mm)
  mass(icn) = mass(icn) + mm
end do

return
end subroutine compute_mass