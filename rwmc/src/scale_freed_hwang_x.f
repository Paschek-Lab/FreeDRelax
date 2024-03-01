      program freedhwang
c
c
c

      implicit none

      integer nmax, nargmax
      parameter (nmax=100000, nargmax=20)
      
      integer i, n, narg
      double precision dd, D, R, b, xmin0, xmin 
      double precision pi, freed, xfreed
      double precision u, t, tmin, tmax, dt, norm0, norm
      double precision scale, g0, g1, g0x, t_trust
      double precision sfac, sx
      
      double precision time(nmax), g2(nmax)
      character*256  arg(nargmax)
      logical getline
      logical full
c-----

      sx=1.0d0
      full=.false.
c     sfac=dsqrt(2.0d0*pi)
      sfac=2.53d0

      narg = 0
      do i=1,nargmax
         call getarg(i,arg(i))
         if (arg(i)(1:1).ne.' ') narg = narg + 1
      enddo

      do i=1,narg
         if (arg(i).eq.'-D') then
            read(arg(i+1),*) D
         endif
         if (arg(i).eq.'-d') then
            read(arg(i+1),*) dd
         endif         
         if (arg(i).eq.'-R') then
            read(arg(i+1),*) R
         endif
         if (arg(i).eq.'-b') then
            read(arg(i+1),*) b
            R=b/2.0d0
         endif
         if (arg(i).eq.'-sfac') then
            read(arg(i+1),*) sfac
         endif
         if (arg(i).eq.'-sx') then
            read(arg(i+1),*) sx
         endif         
         if (arg(i).eq.'-full') then
            full=.true.
         endif    
      enddo

c       print*, '### ', D, dd, R
      
      n=0
      do while(getline(time(n+1),g2(n+1)).and.n.lt.nmax)
         n=n+1
      enddo

      
      pi=4.0d0*datan(1.0d0)

      xmin0=0.0d0
      norm0=(pi/54.0d0)
      
c     xmin=dsqrt(2.0d0*pi)*dd/R

      xmin=sfac*dd/R
      call calc_norm(norm, xmin)

      t_trust=R**2/(D*2.0d0*pi)
      
      print*, 0.0d0, 1.0d0, 1.0d0, 1.0d0, sx
      do i=1,n
         t=time(i)
         u=D*t/dd**2
         if (full) then
            if (t.gt.0.0d0) then
               g0x = xfreed(t,sx,D,dd,xmin0,norm0)               
               g0= freed(t,D,dd,xmin0,norm0)
               g1= freed(t,D,dd,xmin,norm)
               scale=g0/g1
               print*, t, g2(i)*scale, g0, scale, g0x
            endif
         else
            if (t.gt.0.0d0.and.t.lt.t_trust) then
               g0x = xfreed(t,sx,D,dd,xmin0,norm0)                              
               g0= freed(t,D,dd,xmin0,norm0)
               g1= freed(t,D,dd,xmin,norm)
               scale=g0/g1
               print*, t, g2(i)*scale, g0, scale, g0x
            endif
         endif
      enddo
      

      end

      double precision function freed (t, D, dd, xmin,norm)

      implicit none
      
      integer nmax, i
      double precision xmax, xmin, dx, x, x2, t
      double precision D, dd, u, b, sum, pi, norm

      nmax=100000

      pi=4.0d0*datan(1.0d0)
      
      u=D*t/dd**2
c      xmin=sqrt(2.0*pi)*1.0d0/8.0d0
c      xmin=sqrt(8.0)*1.0d0/4.0d0
c      xmin=2.75*1.0d0/4.0d0
c      xmin=0.0
      xmax=40.0d0/sqrt(u)
      dx=xmax/dble(nmax)

     
      sum=0.0d0
c      print*,'#', t, xmax
      do i=1,nmax
         x=xmin+dble(i)*dx
         x2=x*x
         b=x2/(81.0d0+9.0d0*x2-2.0d0*x2**2 + x2**3)

         sum=sum+b*dexp(-x2*u)
c         print*,'f ', x, b*dexp(-x2*u), sum
      enddo
c      freed=sum*dx*18.0/(dd**3)*
c     freed=sum*dx/(3.1415926/54.0)
      freed=sum*dx/norm
      end

     
      
      subroutine calc_norm(norm,xmin)
      
      implicit none

      integer nmax
      parameter (nmax=100000)
      
      integer i
      double precision xmin, xmax, dx, x, x2, sum, b
      double precision norm

      xmax=1000.0
      dx=xmax/dble(nmax)

     
      sum=0.0d0
      do i=1,nmax
         x=xmin+dble(i-1)*dx
         x2=x*x
         b=x2/(81.0d0+9.0d0*x2-2.0d0*x2**2 + x2**3)
         sum=sum+b
      enddo

      norm=sum*dx
      

      end

      logical function getline(xval,yval)
      implicit none
      character*120 input
      double precision xval,yval

      getline=.true.
 99   continue
      read(*,'(a)',end=100, err=100) input
      if (input(1:1).eq.'#'.or.input.eq.' ') goto 99
      read(input,*,err=100,end=100) xval, yval
      return
 100  continue
      getline=.false.
      end

      double precision function xfreed (t, sx, D, dd, xmin,norm)

      implicit none
      
      integer nmax, i
      double precision xmax, xmin, dx, x, x2, t, sx
      double precision D, dd, u, b, sum, pi, norm

      nmax=100000

      pi=4.0d0*datan(1.0d0)
      
      u=D*t/dd**2*sx**(2.0/3.0)
      xmax=40.0d0/sqrt(u)
      dx=xmax/dble(nmax)

     
      sum=0.0d0
      do i=1,nmax
         x=xmin+dble(i)*dx
         x2=x*x
         b=x2/(81.0d0+9.0d0*x2-2.0d0*x2**2 + x2**3)

         sum=sum+b*dexp(-x2*u)
      enddo

      xfreed=sum*dx/norm*sx
      end
