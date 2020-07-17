subroutine output

use config
use elem04, only : stre04
implicit none

integer :: iel, i, j

lpm = 0.d0; stp = 0.d0
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
  call stre04(ndim,ndf,nen,nst,icn,vmatp,ul,xl, lpm,stp)
end do
do i = 1, numnp
  stp(:,i) = stp(:,i)/lpm(i)
end do
ioutput = ioutput + 1
call write_tecplot

return
end subroutine output