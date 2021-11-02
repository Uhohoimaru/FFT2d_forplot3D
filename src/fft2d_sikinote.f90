include "/opt/intel/mkl/include/mkl_dfti.f90"


subroutine fft2d_sikinote(nx,ny,hx,hy,z,mz,mkx,mky,fft_status)

  
    use MKL_DFTI
!$use omp_lib
    implicit none


  integer::i,j,Nx,Ny
  double precision,allocatable::x(:),kx(:)
  double precision,allocatable::y(:),ky(:)
  double precision::hx,xmin,xmax
  double precision::hy,ymin,ymax
 
  double precision,dimension(0:nx-1) :: mkx  
  double precision,dimension(0:ny-1) :: mky  
  complex(kind(0d0)),dimension(0:nx-1,0:ny-1) :: z,mz
  
  double precision,parameter::pi=dacos(-1.d0)
  complex(kind(0d0))::func
  character(*) :: fft_status

  allocate(x(0:Nx-1),kx(0:Nx-1))
  x=0d0; kx=0d0; mkx=0d0
  allocate(y(0:Ny-1),ky(0:Ny-1))
  y=0d0; ky=0d0; mky=0d0
    mz=dcmplx(0d0,0d0)

  do i=0,Nx-1
     x(i)=dble(i)*hx
  enddo
  do j=0,Ny-1
     y(j)=dble(j)*hy
  enddo

  call dftf(Nx,kx,hx); call dftf(Ny,ky,hy) ! get frequency
 
  call dft2d(Nx,Ny,z,fft_status) ! dft fft_status = "forward" or "backward"

  call dfts2d(Nx,Ny,kx,ky,z,mkx,mky,mz) !sort frequency
 


end subroutine fft2d_sikinote



subroutine dftf(N,f,h)
  implicit none
  integer,intent(in)::N
  double precision,intent(in)::h
  double precision,intent(out)::f(0:N-1)

  integer::i
  double precision::mf(0:N-1)
 
  if(mod(N,2).eq.0)then
     do i=0,N-1
        mf(i)=(dble(i+1)-dble(N)*0.5d0)/(dble(N)*h)
     enddo

     do i=0,N-1
        if(i.le.N/2-2)then
           f(i+N/2+1)=mf(i)
        else
           f(i-N/2+1)=mf(i)
        endif
     enddo
  else
     do i=0,N-1
        mf(i)=(dble(i)-dble(N-1)*0.5d0)/(dble(N)*h)
     enddo

     do i=0,N-1
        if(i.le.(N-1)/2-1)then
           f(i+(N-1)/2+1)=mf(i)
        else
           f(i-(N-1)/2)=mf(i)
        endif
     enddo
  endif
 
  return
end subroutine dftf
!--------------------------------
subroutine dft2d(Nx,Ny,z,FB)
  !sikinote
  !date : 2015/12/29
  !developer : sikino
  !CC BY 4.0 (http://creativecommons.org/licenses/by/4.0/deed.ja)
  use MKL_DFTI
  implicit none
  integer,intent(in)::Nx,Ny
  complex(kind(0d0)),intent(inout)::z(0:Nx-1,0:Ny-1)
  character(*),intent(in)::FB

  complex(kind(0d0))::z1(0:Nx*Ny-1)
  integer::Status
  TYPE(DFTI_DESCRIPTOR),POINTER::hand
  integer::id(1:2)
 
  !DFT : Discrete Fourier Transform
  !
  !n    --> number of data.
  !z(i,j) --> value of data at x and y
  !FB   --> "forward"  : Forward DFT
  !     --> "backward" : Backward DFT
     
  call swap1d2d(nx,ny,z,z1,1)
  id(1)=nx; id(2)=ny
 
  Status = DftiCreateDescriptor(hand,DFTI_DOUBLE,DFTI_COMPLEX,2,id)
  Status = DftiCommitDescriptor(hand)

  if(trim(FB).eq."forward")then
     Status = DftiComputeForward(hand,z1)
  elseif(trim(FB).eq."backward")then
     Status = DftiComputeBackward(hand,z1)
  else
     write(6,*)"DFT string different"
     stop
  endif
  Status = DftiFreeDescriptor(hand)
 
  call swap1d2d(nx,ny,z,z1,-1)
 
  if(trim(FB).eq."backward")then
     z=z/dble(nx*ny)
  end if
   
  return
end subroutine dft2d
!---------------------------------------------
subroutine swap1d2d(nx,ny,z2,z1,isign)
  !change array matrix between 1D and 2D
  !  if isign ==  1  --> from 2D to 1D
  !  if isign == -1  --> from 1D to 2D
  implicit none
  integer,intent(in)::nx,ny,isign
  complex(kind(0d0)),dimension(0:nx-1,0:ny-1),intent(inout)::z2
  complex(kind(0d0)),dimension(0:nx*ny-1),intent(inout)::z1
  integer::j1,j2,k

  if(isign.eq.1)then
     do j2=0,ny-1
        do j1=0,nx-1
           k=j2*nx+j1
           z1(k)=z2(j1,j2)
        enddo
     enddo
  elseif(isign.eq.-1)then
     do j2=0,ny-1
        do j1=0,nx-1
           k=j2*nx+j1
           z2(j1,j2)=z1(k)
        enddo
     enddo
  endif
     
  return
end subroutine swap1d2d
!----------------------------
subroutine dfts2d(Nx,Ny,fx,fy,z,mfx,mfy,mz)
  implicit none
  integer,intent(in)::Nx,Ny
  double precision,intent(in)::fx(0:Nx-1),fy(0:Ny-1)
  complex(kind(0d0))::z(0:Nx-1,0:Ny-1)
  double precision,intent(out)::mfx(0:Nx-1),mfy(0:Ny-1)
  complex(kind(0d0)),intent(out)::mz(0:Nx-1,0:Ny-1)
  complex(kind(0d0))::mmz(0:Nx-1,0:Ny-1)
  integer::i,j,k,l  


  if(mod(Ny,2).eq.0)then
     do i=0,Ny-1
        if(i.le.Ny/2)then
           j=i+Ny/2-1
           mfy(j)=fy(i)
           mz(0:Nx-1,j)=z(0:Nx-1,i)
        else
           j=i-Ny/2-1
           mfy(j)=fy(i)
           mz(0:Nx-1,j)=z(0:Nx-1,i)
        endif
     enddo
  else
     do i=0,Ny-1
        if(i.le.(Ny-1)/2)then
           j=i+(Ny-1)/2
           mfy(j)=fy(i)
           mz(0:Nx-1,j)=z(0:Nx-1,i)
        else
           j=i-(Ny-1)/2-1
           mfy(j)=fy(i)
           mz(0:Nx-1,j)=z(0:Nx-1,i)
        endif
     enddo
  endif

  mmz=mz

  if(mod(Nx,2).eq.0)then
     do k=0,Nx-1
        if(k.le.Nx/2)then
           l=k+Nx/2-1
           mfx(l)=fx(k)
           mz(l,0:Ny-1)=mmz(k,0:Ny-1)
        else
           l=k-Nx/2-1
           mfx(l)=fx(k)
           mz(l,0:Ny-1)=mmz(k,0:Ny-1)
        endif
     enddo
  else
     do k=0,Nx-1
        if(k.le.(Nx-1)/2)then
           l=k+(Nx-1)/2
           mfx(l)=fx(k)
           mz(l,0:Ny-1)=mmz(k,0:Ny-1)
        else
           l=k-(Nx-1)/2-1
           mfx(l)=fx(k)
           mz(l,0:Ny-1)=mmz(k,0:Ny-1)
        endif
     enddo
  endif  

  return
end subroutine dfts2d
