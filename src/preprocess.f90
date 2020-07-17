subroutine preprocess

use config
use elem04, only : getsurfinfo
implicit none

! problem initialization
call initialize

! get the surface element information, for surface integral based contact only
if(intmet == 2) call getsurfinfo

! mass matrix
call compute_mass

! output the initial configuration
ioutput = 0
call write_tecplot

return
end subroutine preprocess
