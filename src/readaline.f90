subroutine readaline(ifr,ifw,linestr)
!
!********************************************************************
!
!   input parameters:
!   ifr:     the input file handle
!   ifw:     the output(.ECHO) file handle
!   output parameter:
!   linestr: one line proper input from the input file
!  
!   NOTE:
!    1. the linestr should be long enough
!       e.g. the lenght > 80 in most cases
!    2. the comments in the input file will be written
!       to the output file, but the empty line will be ignored.
!
!********************************************************************
!
      implicit none
      integer ifr,ifw
      character*(*) linestr
!
!      !** local var
!
       logical readOK
       integer length,i
       character*1 c, c_TAB
!
       c_TAB=char(9)
!
       readOK=.false.
!
       do while ( .not. readOK )
          read(ifr,'(A)') linestr
!
	  length=len(linestr)
!
          do 100 i=1,length
             c=linestr(i:i)
             if (    (c.eq.'#') .or. (c.eq.'!') &
               .or. (c.eq.'c') .or. (c.eq.'C')  &
               .or. (c.eq.'*')                  ) then
			  ! comment statement,
			  ! read another line from file
! 
		 write(ifw,'(1x,A80)')  linestr 
!                  exit 
                   goto 101
!
	     elseif ( (c.eq.' ').or. (c.eq.c_TAB) ) then
!                  cycle
		   goto 100         ! space, ignore it.
! 
	     else
		   readOK=.true.    ! this line is proper input  
		   write(ifw,'(1x,A80)')  linestr 
!                  exit
		   goto 101
         endif
 100      continue
 101      continue
!
        enddo ! endwhile
!
return
end subroutine readaline
