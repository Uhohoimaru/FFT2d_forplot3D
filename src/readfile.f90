subroutine readfile

use variables
!$ use omp_lib
implicit none

real(8) :: c


!------------program start

!----read grid

!open(10,file='grid/grid.000001',form='unformatted')
!read(10)jmax,kmax,lmax

!write(*,*)'grid point',jmax,kmax,lmax


!allocate(x(jmax,kmax,lmax))
!allocate(y(jmax,kmax,lmax))
!allocate(z(jmax,kmax,lmax))

!write(*,*)'ndim	mdim'
!write(*,*)ndim,mdim

!read(10)x,y,z

 
!write(*,*)jmax,kmax,lmax
!----------------------------
!-------read input
	a0 = 343.0d0
	theta = 0.0d0
	uj = a0 * 0.9d0
	rho0 = 1.0d0

!	pi = 3.14159265358979323846d0
	theta = theta * 2.0d0*pi/360.0d0
	write(*,*) cos(theta),sin(theta)
   
    ! eg for Tcos0 (change for Tcos1 ...)
    ! Instantaneous fields
 
    nb = 1 
    write(*,*)'flow import start'
    


	
    open(80,file=gridname,form="unformatted")
    read(80)jmax,kmax,lmax

	allocate(x(jmax,kmax,lmax))
	allocate(y(jmax,kmax,lmax))
	allocate(z(jmax,kmax,lmax))

        read(80) (((x(j,k,l),j=1,jmax),k=1,kmax),l=1,lmax),&
         (((y(j,k,l),j=1,jmax),k=1,kmax),l=1,lmax),&
         (((z(j,k,l),j=1,jmax),k=1,kmax),l=1,lmax)
    close(80)




	!data for cylindrical coordinate (instantaneous)

	!-------------------------------------------------
	!data for cartecian cordinate (instantaneous)
	allocate(rho(jmax,lmax,nbs))
	allocate(u(jmax,lmax,nbs))
	allocate(v(jmax,lmax,nbs))
	allocate(w(jmax,lmax,nbs))
	!allocate(e(jmax,lmax))
	allocate(p(jmax,lmax,nbs))
        allocate(buf(jmax,kmax,lmax,5))
	!-------------------------------------------------
        buf = 0.0
        p=0.0d0
	interval = 1



    dx = x(2,1,1) - x(1,1,1)
    dz = z(1,1,2) - z(1,1,1)
    dt = deltat
    do ii=1,nbs !nbinst
    nbg = ii*interval
    write(*,*)nbg,'/',nbs,'ok'
    

    write(filename,'(A,"/flow_z00001_",i5.5)')trim(adjustl(flowinst)),ii

    open(300,file=filename,form='unformatted')
    
    read(300)j,k,l

    if(j.ne.jmax.or.k.ne.kmax.or.l.ne.lmax) then
       write(*,*)'grid not match. end program'
       stop
    end if

    read(300)
    do n=1,5
       read(300) (((buf(j,k,l,n),j=1,jmax),k=1,kmax),l=1,lmax)
    end do

    do j=1,jmax
       do l=1,lmax

          rho(j,l,ii) = dble(buf(j,1,l,1))
          u(j,l,ii) = dble(buf(j,1,l,2))
          v(j,l,ii) = dble(buf(j,1,l,3))
          w(j,l,ii) = dble(buf(j,1,l,4))
          p(j,l,ii) = dble(buf(j,1,l,5))
          
       end do
    end do


!    write(*,*)'maxval(buf)',maxval(buf),'minval(buf)',minval(buf)
!    write(*,*)'maxval(buf)',maxval(p),'minval(buf)',minval(p)

    
 end do
 
 

end subroutine


















