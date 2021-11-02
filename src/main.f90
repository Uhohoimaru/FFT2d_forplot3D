program main
  
  use variables
  
  implicit none
  real(8) :: un_box_ave,dummy
!  complex(kind(0d0)),dimension(nx,ny) :: fft_rx
  
!  call readfile !input the file from plot 3D data


  allocate(coeff(5))

   open(25,file="input.dat")
   read(25,*)
   read(25,*)flowinst
   read(25,*)
 read(25,*)gridname
 read(25,*)
 read(25,*)outfile
 read(25,*)
 read(25,*)nst,nfi,intval !startnum,finishnum,interval
 read(25,*)
 read(25,*)coeff(1),coeff(2),coeff(3),coeff(4),coeff(5)!weight coefficient of rho,u,v,w,p
 read(25,*)
 read(25,*)nr !truncated matrix size
 read(25,*)
 read(25,*)deltat
 nbs = (nfi - nst + 1)/intval
write(*,*)'nbs=',nbs

 flowinst = trim(adjustl(flowinst))
 gridname = trim(adjustl(gridname))
 outfile =  trim(adjustl(outfile))

write(*,*)'flowdirectoryname', trim(adjustl(flowinst))
write(*,*)'gridname', trim(adjustl(gridname))
 write(*,*)'outdirectoryname',trim(adjustl(outfile))
  
 close(25)

!open(55,file=gridname,form='unformatted')
!read(55)jmax,kmax,lmax

!  allocate(x(jmax,kmax,lmax))
!  allocate(y(jmax,kmax,lmax))
!  allocate(z(jmax,kmax,lmax))

!  read(55) (((x(j,k,l),j=1,jmax),k=1,kmax),l=1,lmax),&
!       (((y(j,k,l),j=1,jmax),k=1,kmax),l=1,lmax),&
!       (((z(j,k,l),j=1,jmax),k=1,kmax),l=1,lmax)


!  close(55)





!--------------read file
 call readfile
!------------------------

!------input number of array to nx, ny----
!jmax = 4097
n = jmax ! number of points on z axis
ny = nbs !number of snapshots


! nx = 512 !number of size: spline interpolation
nx = 1024

hx = ( x(n,1,1) - x(1,1,1) )  / dble(nx-1)  ! width of x
hy = dt  ! width of t
write(*,*)"hx hy"
write(*,*)hx,hy



!do ii=1,nbs

!   do l=1,lmax
!      do j=1,jmax
!        p(j,l,ii) = dsin(4.0d0*x(j,1,l)) + dsin( 3.0*hy*dble(ii-1) )
!         p(j,l,ii) =  dsin( x(j,1,l) )  + dsin( hy*dble(ii-1) )
        !!!advection equation
      !   p(j,l,ii) = dsin(x(j,1,l) - hy*dble(ii-1))
       
! end do
!   end do
   
!end do


  allocate(un_box(nx,ny))
un_box = 0.0d0


allocate(f_sp(nx))
allocate(x_sp(nx))
allocate(fdata(n))
allocate(xdata(n))
allocate(fft_rx(nx,ny))
allocate(mf(0:nx-1,0:ny-1))
allocate(psd_n_dm(nx,ny))
f_sp = 0.0d0
x_sp = 0.0d0
fdata = 0.0d0
xdata = 0.0d0

do j=0,nx-1

 x_sp(j+1) = hx*dble(j)

end do



!-----------------------------------------
  


!dx = ( x(jmax,1,1) - x(1,1,1) ) / dble(4096)

!dx = x(2,1,1) - x(1,1,1)
!dt = deltat

 
!  allocate(un_box(nx,ny))
  allocate(kx(0:nx-1),mkx(0:nx-1))
   kx=0d0; mkx=0d0
  allocate(ky(0:Ny-1),mky(0:Ny-1))
   ky=0d0; mky=0d0
  allocate (mz(0:nx-1,0:ny-1))
!  z=dcmplx(0d0,0d0);
 mz=dcmplx(0d0,0d0)
 
    !--------FFT start
    !--------------filtering ------------------
    allocate(filter(nx,ny))



do lz=1,lmax



!!!-------interpolation to the x direction by using spline curve

do ii=1,nbs

 xdata(:) = x(:,1,1)
 fdata(:) = p(:,lz,ii)


! write(*,*)"maxval(fdata)",maxval(fdata)
!  write(*,*)"minval(fdata)",minval(fdata)    
 ! pause
 call spline_1d(nx-1,f_sp,x_sp,n-1,fdata,xdata)
 

!
!  write(*,*)maxval(f_sp)
!  write(*,*)minval(f_sp)    
!  pause
 
   do jx=1,nx

      un_box(jx,ii) = f_sp(jx) 
!* 0.1814 * 1.205 * (343 * 0.9) * (343 * 0.9)
 
   end do


end do


!!!------------------------------------------------------
!  write(*,*)maxval(un_box)
!  write(*,*)minval(un_box)    
!pause







!write(filename,'("mkdir ./result_FFT_mdim")')
!call system(filename)

!!write(filename,'("mkdir ./result_FFT_mdim_lambda")')
!call system(filename)



 !---------substract average of sample
!  do l=1,ny
! un_box_ave = 0.0d0

!     do j=1,nx
     
!        un_box_ave = un_box_ave + un_box(j,l) 

!     end do
!     un_box_ave = un_box_ave / dble(nx)
!     un_box(:,l) = un_box(:,l) - un_box_ave

!  end do
  
  
!  do j=1,nx
! un_box_ave = 0.0d0
!     do l=1,ny
    
!        un_box_ave = un_box_ave + un_box(j,l) 
!     end do
!     un_box_ave = un_box_ave / dble(ny)
!     un_box(j,:) = un_box(j,:) - un_box_ave


! end do




  write(*,*)"max(un_box) min(un_box)"
  write(*,*)maxval(un_box),minval(un_box)




!!!!contain values to un_box(fft container)


  
  !-----------------------------------------
  write(filename,'("./FFT2d_gnuplot/data/subaverage/data_subaverage_z_",f3.1,".txt")')z(1,1,lz)
      open(40,file=filename)
    
    do i=1,nx
       do j=1,ny
       write(40,*)hx*dble(i),hy*dble(j),un_box(i,j) 
       end do
       write(40,*)
    end do
    close(40)
  !------------make hx,hy,rf---------------------
    fft_rx(:,:) = un_box(:,:)

    call fft2d_sikinote(nx,ny,hx,hy,fft_rx,mf,mkx,mky,"forward")

!    call fft2d

    !----------output for gnuplot---------------
    

    !------------------------------------------------------------------------

    !make Power Spectra Density

!    fft_rx = mz
    psd_n_dm(:,:) = cdabs( mf(:,:) ) * cdabs( mf(:,:) ) 

    psd_n_dm(:,:) = ( psd_n_dm(:,:))  / ( hx * hy  )

    write(*,*)"psd_n_dm"
    write(*,*)maxval(psd_n_dm),minval(psd_n_dm)
    
!pause
    !!convert PSD(dB/Hz)
!    psd_n_dm(:,:) = 20d0*log10(psd_n_dm(:,:))
    !!
    
    !------------------------------------------------------------------------


    write(filename,'("./FFT2d_gnuplot/data/FFT_2d_pm3d/fft2d_pm3d_z_",f3.1,".txt")')z(1,1,lz)
    open(21,file=filename)
    do i=0,Nx-1
       do j=0,Ny-1
          
          write(21,'(4e20.10e3)') mkx(i),(-1.0d0)*mky(j),psd_n_dm(i+1,j+1) !attempt to make e^i(kx-wt), so need to convert mky -> -mky 
          !psd_n_dm(i+1,j+1) 
          !write(21,*)i,j,cdabs(mz(i,j))
          
       enddo
       write(21,*)
    enddo
    close(21)


    
    !----ranking-------------------------

    write(filename,'("./FFT2d_gnuplot/Allplot/ranking/z_",f3.1,".txt")')z(1,1,lz)

    open(150,file=filename)
    write(150,*)"wave number frequency     PSD"

    allocate(flg(0:nx-1,0:ny-1))
    flg = 0
   do n_itr=1,30
      max_psd = 0.0d0
    do i=0,Nx-1
       do j=0,Ny-1

          if(max_psd < psd_n_dm(i+1,j+1) .and. flg(i,j)==0 .and. -1.0d0 * mky(j) >= 0.0d0 ) then
         max_psd =  psd_n_dm(i+1,j+1)
         i_flg = i
         j_flg = j
          end if

       end do
    end do

    write(150,'(3e20.4)')mkx(i_flg),-1.0d0*mky(j_flg),max_psd
    flg(i_flg,j_flg) = 1 

    end do

    deallocate(flg)
    !--------------------------------------------



    !---filter only minus k space
    filter(:,:) = 1.0d0

    do l=0,ny-1
    do j=0,nx-1
       if( ( mkx(j) <= 0.0d0 .and. -1.0d0*mky(l) <= 0.0d0 ) .or. ( mkx(j) >= 0.0d0 .and. -1.0d0*mky(l) >= 0.0d0) )then!
!       if(  (mkx(j) >= 0.0d0) .or. -1.0d0*mky(l) <= 0.0d0  )then!

! .or. ( mkx(j) <= 0.0d0.and.mky(l) <= 0.0d0)) then

          fft_rx(j+1,l+1) = (0.0d0,0.0d0)
          filter(j+1,l+1) = (0.0d0,0.0d0)
             
          end if
          
    end do
    end do
        !---------------------------
    
       open(40,file='filter.txt')
    
    do i=1,nx
       do j=1,ny
       write(40,*)mkx(i-1),(-1.0d0)*mky(j-1),dble(filter(i,j)) !attempt to make e^i(kx-wt), so need to convert mky -> -mky 
       end do
       write(40,*)
    end do
    close(40)
    
    !-----------------------------------------
   ! call fft2d_inv

    call fft2d_sikinote(nx,ny,hx,hy,fft_rx,mf,mkx,mky,"backward")

    !---------------output the file with pm3d format (for gnuplot)------------
 write(filename,'("./FFT2d_gnuplot/data/FFT2d_recon/FFT2d_reconstruct_only_km_pm3d",f3.1,".txt")')z(1,1,lz)
 
    open(21,file=filename)
    do i=0,nx-1
       do j=0,ny-1
          write(21,'(4e20.10e3)')hx*dble(i),hy*dble(j),dble(fft_rx(i+1,j+1))
          !write(21,*)i,j,cdabs(mz(i,j))
          
       enddo
       write(21,*)
    enddo
    close(21)
    
    !------------------------------------------------------------------------


   write(filename,'("mkdir ./FFT_reconstruct/data_z_",f3.1)')z(1,1,lz)
   call system(filename)


    do j=0,ny-1

   
   write(filename,'("./FFT_reconstruct/data_z_",f3.1,"/z_0_",i5.5)')z(1,1,lz),j

   open(10,file=filename)

   

   do i=0,nx-1
            
      write(10,*)dble(i)*hx,dble(fft_rx(i+1,j+1))
   
   end do
   
      

close(10)

 
end do

  
!    call system("sh makeplot.sh")


end do





end program


function func(x,y)
  implicit none
  double precision::x,y
  double precision::func
 
  func=dsin(4d0*x) + dsin(3d0*y) 
 
  return
end function func














