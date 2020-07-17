subroutine write_vtu

use config
implicit none
integer, parameter :: plu = 2015
character :: lct*16

integer :: i,j,ii,jj,ip

lct = 'PARAVIEW_000.vtu'
if(ioutput < 10) then
  write(lct(12:12),'(i1)') ioutput
elseif(ioutput < 100) then
  write(lct(11:12),'(i2)') ioutput
elseif(ioutput < 1000) then
  write(lct(10:12),'(i3)') ioutput
else
  write(*,*) 'Too many frames in write_vtu!'
  stop
endif
write(*,*) 'Saving PARAVIEW data to ',lct
open(unit=plu,file=lct,access='sequential')

!Write out top header
write(plu,1000)

write(plu,1020) numnp,numel               ! Start Mesh/Piece Section
write(plu,1010) '<Points>'                ! Start Point/Node data

write(plu,1030) 'Float64','nodes', 3      ! ndm
do ip = 1,numnp
  write(plu,2000) (xr(j,ip)+disp((ip-1)*ndf + j), j = 1, ndf), (0.d0,ii=ndf+1,3) 
end do ! i
write(plu,1010) '</DataArray>'            ! Close Node data

write(plu,1010) '</Points>'               ! Close Points section

write(plu,1010) '<Cells>'                 ! Start Cell Section
write(plu,1030) 'Int32','connectivity',1  ! Start Elements

do i = 1,numel
  write(plu,'(4i8)') (cn(ii,i)-1,ii=1,nen)
end do ! i

write(plu,1010) '</DataArray>'           ! Close Elements

write(plu,1030) 'Int32','offsets',1      ! Start Offsets
write(plu,2010) (nen*i, i = 1,numel)
write(plu,1010) '</DataArray>'            ! Close Offsets

write(plu,1030) 'Int32','types',1         ! Start Element types

write(plu,2010) (9,jj=1,numel)

write(plu,1010) '</DataArray>'               ! Close Element types
write(plu,1010) '</Cells>'                   ! Close Cell Section
                                              
write(plu,1010) '<PointData>'                ! Start Point Data

write(plu,1030) 'Float64','Displ',ndf        ! Start Displacements
do ip = 1,numnp
  write(plu,2000) (disp((ip-1)*ndf + j), j = 1, ndf)
end do ! i
write(plu,1010) '</DataArray>'               ! Close Displacements

! if(np(42).ne.0) then
!   write(plu,1030) 'Float64','Velocity',ndf   ! Start Velocity
!   do i = 1,numnp
!     write(plu,2000) (hr(npud+(i-1)*ndf+ii),ii = 0,ndf-1)
!   end do ! i
!   write(plu,1010) '</DataArray>'             ! Close Velocity
! 
!   write(plu,1030) 'Float64','Acceleration',ndf ! Start Acceleration
!   do i = 1,numnp
!     write(plu,2000) (hr(npud+nneq+(i-1)*ndf+ii),ii = 0,ndf-1)
!   end do ! i
!   write(plu,1010) '</DataArray>'              ! Close Acceleration
! endif
! 
! if(abs(istv).gt.0) then
!   write(plu,1030) 'Float64','Stress',abs(istv)   ! Start Stresses
!   do i = 1,numnp
!     write(plu,2000) (hr(npnp+(i-1) + ii*numnp),ii=1,abs(istv))
!   end do
!   write(plu,1010) '</DataArray>'              ! Close Stresses
! 
!   write(plu,1030) 'Float64','PStress',7       ! Start Principal Stresses
!   do i = 1,numnp
!     write(plu,2000) (hr(nper+(i-1) + ii*numnp),ii=1,7)
!   end do ! i
!   write(plu,1010) '</DataArray>'              ! Close Stresses
! 
! else
!   write(*,*) ' No stresses output to Paraview file'
! endif

write(plu,1030) 'Float64','Stress',2   ! Start Stresses
do ip = 1,numnp
  write(plu,'(2e13.5)')  (stp(i,ip), i=1,2)
end do
write(plu,1010) '</DataArray>'         ! Close Stresses

write(plu,1010) '</PointData>'              ! Close Point Data Section

write(plu,1010) '</Piece>'                  ! Close Mesh/Piece

!Close the XML file

write(plu,1010) '</UnstructuredGrid> </VTKFile>'

close(plu, status = 'keep')
  
!     I/O Formats

1000  format('<VTKFile type="UnstructuredGrid" version="0.1">',/'<UnstructuredGrid>')

1010  format(a)

1020  format('<Piece NumberOfPoints="',i10,'" NumberOfCells="',i10,'">')

1030  format('<DataArray type="',a,'" Name="',a,'" NumberOfComponents="',i2,'" format="ascii">')

!2000  format(1p,6e13.5)
2000  format(1p,3e13.5)

2010  format(i8)

return

end subroutine write_vtu