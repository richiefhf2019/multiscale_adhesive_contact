subroutine solve


use config
implicit none

integer :: iter, i
real(8) :: lambda, dlambda
real(8) :: rdu, prf

integer :: nou, niter, nloop, status, steperr
real(8) :: alpha,beta,eta,dzeta,r0rn, norm0

lambda = 0.d0; dlambda = 0.01d0
! open file channel for prescribed f and u
open(100, file = 'pus_slide.plt')
write(100,*) 'Variables = "u","P(u)"'
!----------------------------------------------------------------------------------------
! start of the pseudo time (load) loop, quasi static simulation
!----------------------------------------------------------------------------------------
do while (lambda <= 1.d0)
  iter = 0
  write(*,*) '----------------------------------------------------------------------'
  write(*,*) 'Current loading ratio', lambda
  write(*,*) 'Current loading increment', dlambda
  write(*,*) '----------------------------------------------------------------------'
  lambda = lambda + dlambda
  disp(idd) = dud(:)*lambda*3.d0
  if(abs(lambda - 0.15d0)<0.001) dlambda = dlambda/2.d0
  if(abs(lambda - 0.85d0)<0.001) dlambda = dlambda*2.d0
  !************************************
  do ! start of newton's iteration
    iter = iter + 1
    
    ! form the stifness matrix ktot and force vector ftot
    call formkf

    ! residual of the current iteration
    residual = fext*lambda - ftot

    ! stiffness related to the free dofs
    kff = ktot(idf,idf)

    ! get residual at the free dofs, R^{(i+1)}
    rf = residual(idf)

    ! solve for the displacement inrecment \Delta d^{(i)}
    ! pbicg, solve Ax = b
    dxf = 0.d0; dxfr = 0.d0; nou = 0
    do
      call pbicg(dxf,dxfr,rf,nequ,nequ,12,nou,wrk,nloop,100 ,1.d-6 &
    	  ,norm0,alpha,beta,eta,dzeta,r0rn,status,steperr)
      if (status.lt.0) then
         write(*,*) 'stop nstat',status,steperr
         stop
      endif
      dxfr = matmul(kff,dxf)
      if(status == 1) then
        !write(*,'(A,i3,A,i3)') 'for iload = ', iload, '# of bicg iteration =', nloop
        exit
      endif
    end do
    if (steperr == 0) then
       write(*,*) 'nloop has reached maxit',nloop
    endif

    disp(idf) = disp(idf) + dxf ! d^{(i+1)} = d^{(i)} + \Delta d^{(i)

    ! check the energy residual
    rdu = dot_product(rf,dxf)

    ! convergence dectection, to jump out the load loop
    if (rdu <= tolerance) then
      write(*,*) 'After Newton iterations :', iter
      write(*,*) 'The energy residual becomes', rdu
      exit
    endif
    
    if( mod(iter,20) == 0) then
      write(*,*) 'iter =', iter, 'not converged'
      write(*,'(A, e15.5)')'with energy-residual =', rdu
      write(*,'(A, e15.5)')'dispalcement increment', sqrt(dot_product(dxf,dxf))
      if(iter == 60) exit
    endif
  end do ! end of the newton iteration for load step iload

  ! output intermediate result for visualization
  call output
  ! output the reaction force and prescribed displacement at the top
  prf = 0.d0
  do i = 1, size(idd)
    if(mod(i,1) == 0) then
      prf = prf - residual(idd(i))
    endif
  end do
  write(100,'(2e20.10)') lambda*1.d0-0.2d0, prf; call flush(100)
!----------------------------------------------------------------------------------------
! end of the pseudo time (load) loop, quasi static simulation
!----------------------------------------------------------------------------------------
end do

close(100)

return
end subroutine solve