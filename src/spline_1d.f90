subroutine spline_1d(M,f,x,N,fdata,xdata)
  implicit none
  integer,intent(in)::N,M
  double precision,intent(in)::x(0:M),fdata(0:N),xdata(0:N)
  double precision,intent(out)::f(0:M) !,df(0:M),df2(0:M)
 
  integer::i,j
  double precision::t
  double precision,allocatable::a(:),b(:),c(:),d(:)
 
  allocate(a(0:N-1),b(0:N-1),c(0:N-1),d(0:N-1))
  call spline_abcd(a,b,c,d,N,xdata,fdata)
 
  do j=0,M
     if(x(j).lt.xdata(0))then
        t=x(j)-xdata(0)
        f(j)  = ((a(0)*t+b(0))*t+c(0))*t+d(0)
        !df(j) =3d0*a(0)*t*t+2d0*b(0)*t+c(0)
        !df2(j)=6d0*a(0)*t+2d0*b(0)
     elseif(x(j).gt.xdata(N))then
        t =x(j)-xdata(N-1)
        f(j)  = ((a(N-1)*t+b(N-1))*t+c(N-1))*t+d(N-1)
        !df(j) =3d0*a(N-1)*t*t+2d0*b(N-1)*t+c(N-1)
        !df2(j)=6d0*a(N-1)*t+2d0*b(N-1)
     else
        do i=0,N-1
           if(x(j).ge.xdata(i).and.x(j).le.xdata(i+1))then
              t=x(j)-xdata(i)
              f(j)  = ((a(i)*t+b(i))*t+c(i))*t+d(i)
              !df(j) =3d0*a(i)*t*t+2d0*b(i)*t+c(i)
              !df2(j)=6d0*a(i)*t+2d0*b(i)
              exit
           endif
        enddo
     endif
  enddo


!  open(550,file="test_spline.txt")
!  do i=0,m
!  write(550,*)x(i),f(i)
!!  end do
!  close(550)

!  open(550,file="test_original.txt")
!  do i=0,n
!  write(550,*)xdata(i),fdata(i)
!  end do
!  close(550)
!pause
  return
end subroutine spline_1d



subroutine spline_abcd(a,b,c,d,N,xdata,fdata)
  implicit none
  integer,intent(in)::N
  double precision,intent(out)::a(0:N-1),b(0:N-1),c(0:N-1),d(0:N-1)
  double precision,intent(in)::xdata(0:N),fdata(0:N)
 
  integer::i,j
  double precision,allocatable::l(:),mu(:),h(:),alpha(:),z(:)

  allocate(h(0:N-1),alpha(1:N-1))

  do i=0,N-1
     h(i)=xdata(i+1)-xdata(i)
  enddo
  do i=1,N-1
     alpha(i)=3d0*(fdata(i+1)-fdata(i))/h(i)-3d0*(fdata(i)-fdata(i-1))/h(i-1)
  enddo
   
  allocate(l(0:N-1),mu(0:N-1),z(0:N-1))
  l=0d0;  mu=0d0;  z=0d0
 
  l(0)=1d0
  do i=1,N-1
     l(i)=2d0*(xdata(i+1)-xdata(i-1))-h(i-1)*mu(i-1)
     mu(i)=h(i)/l(i)
     z(i)=(alpha(i)-h(i-1)*z(i-1))/l(i)
  enddo

  b(N-1)=z(N-1)
  c(N-1)=(fdata(N)-fdata(N-1))/h(N-1)-h(N-1)*(2d0*b(N-1))/3d0
  a(N-1)=-b(N-1)/(3d0*h(N-1))
  d(N-1)=fdata(N-1)
 
  do j=N-2,0,-1
     b(j)=z(j)-mu(j)*b(j+1)
     c(j)=(fdata(j+1)-fdata(j))/h(j)-h(j)*(b(j+1)+2d0*b(j))/3d0
     a(j)=(b(j+1)-b(j))/(3d0*h(j))
     d(j)=fdata(j)
  enddo
 
  return
end subroutine spline_abcd
