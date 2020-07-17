program main
!
!-------------------------------------------------------------------------------------------------!
! This program is developed to study the Coarse-Grained Contact Algorithm, in Quasi-Static Case
! Two features, the original body-body integral and the newly developed surface-surface integral
! are implemented. Updated Lagrange Formulations are adopted, for the sake of efficiency.
! October 1st, 2014, Civil Engineering, UC Berkeley, by Houfu Fan
!-------------------------------------------------------------------------------------------------!
!
implicit none

! preprocess
call preprocess
! slove
call dynamic_solve

stop
end program



