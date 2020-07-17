subroutine dynamic_solve

use config
implicit none

real(8) :: time

integer :: iter, i, j, iload
real(8) :: lambda, dlambda
real(8) :: rdu, prf

integer :: istep = 0

lambda = 0.d0; dlambda = 0.01d0
! open file channel for prescribed f and u
open(100, file = 'pus_slide.plt')
write(100,*) 'Variables = "u","P(u)"'
! provide the intital velocity
do i = 1, numnp
  if(ntype(i) == 1) then
   vel_half((i-1)*ndf+1) = 2.d0
   vel_half((i-1)*ndf+2) = 0.d0
   !vel_half((i-1)*ndf+1) = 2.0d0
  endif
end do
! calculate the force vector and hence the acceleration
call compute_force
do i = 1, numnp
  do j = 1, ndf
    acc((i-1)*ndf+j) = -ftot((i-1)*ndf+j)/mass(i)
  end do
end do
time = 0.d0
open(200, file = 'energy.dat')
!----------------------------------------------------------------------------------------
! start of the explicit time integration
!----------------------------------------------------------------------------------------
do iload = 1, nload
  ! current time
  time = time + dt; istep = istep + 1
  ! get the half step velocity
  do i = 1, numnp
    if(ntype(i) == 3) cycle
    do j = 1, ndf
      vel_half((i-1)*ndf+j) = vel_half((i-1)*ndf+j) + acc((i-1)*ndf+j)*dt
    end do
  end do
  ! get the full step displacement
  do i = 1, numnp
    if(ntype(i) == 3) cycle
    do j = 1, ndf
      disp((i-1)*ndf+j) = disp((i-1)*ndf+j) + vel_half((i-1)*ndf+j)*dt
    end do
  end do
  ! prescribe the boundary condition
  disp(idd) = dud(:)*dble(iload)/dble(nload)
  ! calculate the force vector and hence the acceleration
  call compute_force
  do i = 1, numnp
    do j = 1, ndf
      acc((i-1)*ndf+j) = -ftot((i-1)*ndf+j)/mass(i)
    end do
  end do

  ! get the full step velocity
  do i = 1, numnp
    if(ntype(i) == 3) cycle
    do j = 1, ndf
      vel((i-1)*ndf+j) = vel_half((i-1)*ndf+j) + acc((i-1)*ndf+j)*0.5d0*dt
    end do
  end do

  ! compute the kinetic energy
  energy_kin = 0.d0
  do i = 1, numnp
    do j = 1, ndf
      energy_kin = energy_kin + 0.5d0*vel((i-1)*ndf+j)**2*mass(i)
    end do
  end do

  ! output the intermediate result and specific intervals
  if(mod(istep,iprint) == 0) then
    call output; call flush(ivs)
	call write_vtu
    write(*,*) 'current time step', istep
    write(200,'(5e20.10)') time, energy_int, energy_cont, energy_kin, energy_int+energy_cont+energy_kin
    !write(200,'(i5, 4e20.10)') ioutput, energy_int, energy_cont, energy_kin, energy_int+energy_cont+energy_kin
  endif
end do

return
end subroutine dynamic_solve
