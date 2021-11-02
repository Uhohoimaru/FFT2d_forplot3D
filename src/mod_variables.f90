module variables

    real(8) :: thetagT0, deltat, r0, rho0, p0, a0, deltat_nondim, theta,uj,dt,max_psd
    real(8) :: uu, vv, ww, uvw, e, ds, dx ,dz, len_x, len_y
    real(8), allocatable,dimension(:)    :: rg,zg
    real(8), allocatable,dimension(:,:)  ::rhoT0moy,urT0moy,utT0moy,uzT0moy,pT0moy
    real(8), allocatable,dimension(:,:)  :: vortT0moy,dilatT0moy
    real(8), allocatable,dimension(:,:)  :: rhoT0,urT0,utT0,uzT0,pT0,eT0
    real(8), allocatable,dimension(:,:)  :: vortT0,dilatT0

    real(8), allocatable ::  work(:),un_box(:,:),x(:,:,:),y(:,:,:),z(:,:,:)
 
    real(8),allocatable :: rho(:,:,:),u(:,:,:),v(:,:,:),w(:,:,:),p(:,:,:),coeff(:)
    
    real(8),allocatable :: r_FFT(:,:),f_sp(:),x_sp(:),fdata(:),xdata(:),psd_n_dm(:,:)


    integer,allocatable :: flg(:,:)
    
    real(8) xp,yp,eps,maxPSD
    real(8),allocatable :: buf(:,:,:,:)
    integer i,j,k,l,ii,n,nbs,m,mode,jkn,itr,nitr, nj, nl, jx, lz, i_flg, j_flg
    integer jmax,kmax,lmax,ary,nc,nst,nfi,intval, n_itr
    integer lda,ldvt,ldu,lwork,info,nmode,ndim,mdim
    integer jpa, lpa, nz, nr, nbinst, interval, nb, nbg,nx,ny
   
    complex(kind(0d0)),allocatable :: filter(:,:), mf(:,:), fft_rx(:,:),rx(:)
    complex(kind(0d0)) :: rx_ave
    
    
!--test
    
    double precision,allocatable::kx(:),mkx(:)
    double precision,allocatable:: ky(:),mky(:)
    double precision::hx,xmin,xmax
    double precision::hy,ymin,ymax
    complex(kind(0d0)),allocatable::mz(:,:)
    
  double precision,parameter::pi=dacos(-1.d0)
  ! complex(kind(0d0))::func
  double precision::func
  external::func
  
!----

character*1024 :: filename,flowinst,gridname,outfile
  

  data gamma,nbs,nmode,eps/1.4e0,2600,2600,1.0d-2/ !/number of samplings,heat capacity ratio,number of modes 




end module
