      subroutine pbicg(xi,xr,b,lda,ndim,nlar,nou,wrk,itno,maxit ,tol &
	  ,norm,alpha,beta,eta,dzeta,r0rn,status,steperr)
!     We want to solve Ax=b
!     mat:  matrix A
!     lda: Input: leading dimension array of the matrix
!     ndim: Input: dimension  of the matrix: ndim.le.lda
!     nlar: Input: size of the work vector: nlar.ge.12
!     MAXIT: Input: Maximum of iteration for the iterative solver
!     NLOOP: Output: number of iteration for the iterative solver Should
!     be initialize to 0.
!     ncompte: number of Ax products.
!     STATUS: Output: if STATUS.lt.0 a problem occured in GPBICG
!     STATUS=1 the requested tolerance has been reached or the maximum number of
!     iterations has been reached
!     STEPERR: Output: if STEPERR.gt.0: indicates where the problem
!     occurs in PBICG. STEPERR=0 the maximum number of iterations has
!     been reached.  STEPERR=-1 routine completed without error
!     tol: Input: tolerance requested by the user. At the end we have:
!     r=||Ax-b||/||b||<tol
!     b: Input: right-hand side of Ax=b 
!     norm: Output: norm of b
!     xi: Input: initial guess; output:solution of the linear equation
!     xr: Input: xr = A xi
!     NOU: Local integer used by PBICG. Should be initialized to 0.
!     ALPHA,BETA,ETA,DZETA,R0RN: Local complex needs for GPBICG 
!     WRK: local array used by for PBICG
      implicit none

      integer :: ii,nou,ndim,itno,lda,maxit,status,steperr,nlar
!      double complex b(lda),wrk(lda,nlar),xi(lda),xr(lda)	  
      real(8) :: b(lda),wrk(lda,nlar),xi(lda),xr(lda)

!     .. local scalars ..
!      double complex alpha,beta,eta,dzeta,r0rn,ctmp,ctmp1,ctmp2,ctmp3,ctmp4,ctmp5
      real(8) :: alpha,beta,eta,dzeta,r0rn,ctmp,ctmp1,ctmp2,ctmp3,ctmp4,ctmp5	  
!      double precision tol,norm,residu
      real(8) :: tol,norm,residu	  

      if (nou.eq.1) goto 10
      if (nou.eq.2) goto 20
      if (nou.eq.3) goto 30
      if (nou.eq.4) goto 40

      status = 0
      steperr = -1

!     column index of the various variables stored in wrk array.
!     1:r0
!     2:p
!     3:r
!     4:y
!     5:t
!     6:ap
!     7:at
!     8:u
!     9:w
!     10:z
!     11:x
!     12:t old

!     compute norm and ax0 (x0 initial guess)
      norm=0.d0
      do ii=1,ndim
         wrk(ii,11)=xi(ii)
         norm=norm+abs(b(ii))**2.d0
      enddo
      norm=sqrt(norm)

      nou=1
      return

!     initialize r0=b-ax0,rot=ro,p0=v0=d0
 10   do ii=1,ndim
         wrk(ii,1)=b(ii)-xr(ii)
         wrk(ii,2)=0.d0
         wrk(ii,3)=wrk(ii,1)
         wrk(ii,5)=0.d0
         wrk(ii,9)=0.d0
         wrk(ii,8)=0.d0
         wrk(ii,10)=0.d0
      enddo
!     compute the initial residue
      ctmp=0.d0
      do ii=1,ndim
         ctmp=ctmp+wrk(ii,1)*(wrk(ii,1))
      enddo
      !write(*,*) 'residual initial',abs(sqrt(ctmp))/norm

      r0rn=ctmp

!     initialize rho,alpha,w=1,tau=norm,theta=0,eta=0
      beta=0.d0
      
!     begin the iteration sequence
      itno=-1
 100  itno=itno+1

!     compute p=r+beta*(p-u)
      do ii=1,ndim
         wrk(ii,2)=wrk(ii,3)+beta*(wrk(ii,2)-wrk(ii,8))
         xi(ii)=wrk(ii,2)
      enddo
      nou=2
      return

!     compute ap
 20   do ii=1,ndim
         wrk(ii,6)=xr(ii)
      enddo

!     compute alpha=r0r/r0ap
      ctmp=0.d0
      do ii=1,ndim
         ctmp=ctmp+(wrk(ii,1))*wrk(ii,6)
      enddo

      if (ctmp.eq.0.d0) then
         status=-1
         steperr=1
         return 
      endif

      alpha=r0rn/ctmp

!     compute y=t-r-alpha*w+alpha*ap et de t=r-alpha ap
      do ii=1,ndim
         wrk(ii,4)=wrk(ii,5)-wrk(ii,3)-alpha*wrk(ii,9)+alpha*wrk(ii,6)
         wrk(ii,12)=wrk(ii,5)
         wrk(ii,5)=wrk(ii,3)-alpha*wrk(ii,6)
         xi(ii)=wrk(ii,5)
      enddo
      nou=3
      return

!     compute at
 30   do ii=1,ndim
         wrk(ii,7)=xr(ii)
      enddo

!     compute dzeta and eta

      if (itno.eq.0) then

        eta=0.d0
        dzeta=0.d0
        ctmp=0.d0
        do ii=1,ndim
         dzeta=dzeta+(wrk(ii,7))*wrk(ii,5)
         ctmp=ctmp+(wrk(ii,7))*wrk(ii,7)
        enddo
        
        if (ctmp.eq.0.d0) then
           status=-1
           steperr=2
           return 
        endif
        dzeta=dzeta/ctmp

      else
         
         ctmp1=0.d0
         ctmp2=0.d0
         ctmp3=0.d0
         ctmp4=0.d0
         ctmp5=0.d0
         
         do ii=1,ndim
            ctmp1=ctmp1+(wrk(ii,7))*wrk(ii,7)
            ctmp2=ctmp2+(wrk(ii,4))*wrk(ii,4)
            ctmp3=ctmp3+(wrk(ii,7))*wrk(ii,4)
            ctmp4=ctmp4+(wrk(ii,7))*wrk(ii,5)
            ctmp5=ctmp5+(wrk(ii,4))*wrk(ii,5)
         enddo

         ctmp=ctmp1*ctmp2-ctmp3*(ctmp3)
         
         if (ctmp.eq.0.d0) then
            status=-1
            steperr=3
            return 
         endif

         dzeta=(ctmp2*ctmp4-ctmp5*ctmp3)/ctmp
         eta=(ctmp1*ctmp5-(ctmp3)*ctmp4)/ctmp

      endif

!     compute u
      do ii=1,ndim
         wrk(ii,8)=dzeta*wrk(ii,6)+eta*(wrk(ii,12)-wrk(ii,3)+beta*wrk(ii,8))
      enddo

!     compute z
      do ii=1,ndim
         wrk(ii,10)=dzeta*wrk(ii,3)+eta*wrk(ii,10)-alpha*wrk(ii,8)
      enddo

!     compute x and r
      residu=0.d0
      do ii=1,ndim
         wrk(ii,11)=wrk(ii,11)+alpha*wrk(ii,2)+wrk(ii,10)
         wrk(ii,3)=wrk(ii,5)-eta*wrk(ii,4)-dzeta*wrk(ii,7)
         residu=residu+abs(wrk(ii,3))**2.d0
      enddo
      residu=sqrt(residu)/norm

      if (residu.le.tol) then
         status=1
         do ii=1,ndim
            xi(ii)=wrk(ii,11)
         enddo
         nou=4
         return
 40   endif

!     compute beta
      ctmp=0.d0
      do ii=1,ndim
         ctmp=ctmp+(wrk(ii,1))*wrk(ii,3)
      enddo

      if (r0rn.eq.0.d0) then
         status=-1
         steperr=4
         return 
      endif

      beta=alpha*ctmp/dzeta/r0rn
      r0rn=ctmp

!     compute w
      do ii=1,ndim
         wrk(ii,9)=wrk(ii,7)+beta*wrk(ii,6)
      enddo

      if (itno.le.maxit) goto 100

      status=1
      steperr=0
      do ii=1,ndim
         xi(ii)=wrk(ii,11)
      enddo

      end subroutine pbicg