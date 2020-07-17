subroutine write_tecplot

use config
implicit none

integer :: ip, iel, i

10    FORMAT('Zone T="Load step 0",N=',I6,',E=',I6,',Datapacking=point,zonetype=fequadrilateral')
11    FORMAT('Zone T="',I4,'-th LOAD STEP"',',N=',I6,',E=',I6,', &
              Datapacking=point,zonetype=fequadrilateral,connectivitysharezone=1')

if(ioutput == 0) then
  write(ivs,'(a)') 'Title="Quasi-Static Contact Simulation"'
  write(ivs,'(a)') 'Variables="X","Y","UX","UY","S","I<sub>1<sub>","MAT"'
  write(ivs,10) numnp, numel
  do ip = 1, numnp
    write(ivs,'(1x,6e15.6,I8)') ( xr(i,ip),i=1,ndim ), ( disp( (ip-1)*ndf + i ), i = 1, ndf ), (  stp(i,ip), i = 1, 2 ), ntype(ip)
  end do
  do iel = 1, numel
    write(ivs,'(1x,4i8)') (cn(i,iel), i = 1, nen)
  end do
else
  write(ivs,11) ioutput, numnp, numel
  do ip = 1, numnp
    write(ivs,'(1x,6e15.6,I8)') xr(1,ip)+disp((ip-1)*ndf + 1),  xr(2,ip)+disp((ip-1)*ndf + 2), &
     ( disp( (ip-1)*ndf + i ), i = 1, ndf ), ( stp(i,ip), i = 1, 2 ), ntype(ip)
  end do
endif

return
end subroutine write_tecplot