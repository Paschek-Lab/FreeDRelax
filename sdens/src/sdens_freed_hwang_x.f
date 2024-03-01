      program freedhwang
c
c
c

      implicit none

      integer nmax, nargmax
      parameter (nmax=100000, nargmax=100)
      
      integer i, it, n, nt, narg
      integer relax, domhz
      
      double precision dd, D, R, b, xmin0, xmin, sx
      double precision pi, freed, xfreed, dt
      double precision o, oout, ou, u, omegau, t, tmin, tmax, norm0, norm
      double precision scale, jomega, jomega2, o2, ou2, fac
      double precision j0, fm, tauadd, deltaj, deltaj2, rho
      double precision gamma
      
      double precision omega(nmax), time(nmax), deltag(nmax)
      character*256  arg(nargmax)
      logical getline
c-----

      sx=1.0d0
      rho=1.0d0
      scale=1.0d0
      relax=0
      domhz=0
      gamma=267.5153151
      
      
      open(10,file='diff.dat')
      do i=1,nmax
         nt=i
         read(10,*,end=100) time(i), deltag(i)
      enddo
 100  continue
      nt=nt-1
      close(10)
c      print*,'#',nt
      
      narg = 0
      do i=1,nargmax
         call getarg(i,arg(i))
         if (arg(i)(1:1).ne.' ') narg = narg + 1
      enddo

      do i=1,narg
         if (arg(i).eq.'-1H') then
            gamma=267.5153151
         endif
         if (arg(i).eq.'-19F') then
            gamma=251.815
         endif                  
         if (arg(i).eq.'-D') then
            read(arg(i+1),*) D
         endif
         if (arg(i).eq.'-d') then
            read(arg(i+1),*) dd
         endif         
         if (arg(i).eq.'-relax') then
            relax=1
         endif
         if (arg(i).eq.'-mhz') then
            domhz=1
         endif
         if (arg(i).eq.'-scale') then
            read(arg(i+1),*) scale
         endif
         if (arg(i).eq.'-sx') then
            read(arg(i+1),*) sx
         endif
         if (arg(i).eq.'-rho') then
            read(arg(i+1),*) rho
         endif
      enddo

c       print*, '### ', D, dd, R
      
      n=0
      do while(getline(omega(n+1)).and.n.lt.nmax)
         if (domhz.eq.1) then
            omega(n+1)=omega(n+1)*1.0d-6
         endif
         n=n+1
      enddo

      
      pi=4.0d0*datan(1.0d0)

      xmin0=0.0d0
      norm0=(pi/54.0d0)



      fac=2.0d0
c      fac=fac*(267.5153151e6)**4
      fac=fac*(gamma*1.0e6)**4      
      fac=fac*(6.62607015e-34/(2*pi))**2
      fac=fac*1.0d0/2.0d0*(3.0d0/2.0d0)
      fac=fac*(1.25663706212e-6/(4.0d0*pi))**2
c     fac=fac*2*4096/(4.96735e-9)**3
      fac=fac*rho*(1.0e9)**3
      fac=fac*4.0d0*pi/(3.0d0*(dd*1e-9)**3)
      fac=fac*1.0e-12


      fm=pi*sqrt(2.0d0)/15.0d0
c      fm=fm*(267.5153151e6)**4      
      fm=fm*(gamma*1.0e6)**4
      fm=fm*(6.62607015e-34/(2*pi))**2
      fm=fm*(1.25663706212e-6/(4.0d0*pi))**2
      fm=fm*(1+4.0d0*dsqrt(2.0d0))
c     fm=fm*2*4096/(4.96735e-9)**3
      fm=fm*rho*(1.0e9)**3
      fm=fm/dsqrt((D*1.0e-6)**3)


      dt=time(2)-time(1)
      j0=(4.0d0/9.0d0*dd**2/D)*sx**(1.0/3.0)

      deltaj=0.0d0
      do it=1,nt-1
         deltaj=deltaj+
     .        0.5d0*(deltag(it)+deltag(it+1))
     .        *(time(it+1)-time(it))
      enddo
      deltaj=deltaj*scale

      if (relax.eq.1) then
         print*, 0.0,
     .        (j0+deltaj)*fac,
     .        j0*fac,
     .        deltaj*fac,
     .        j0*fac         
      else
         print*, 0.0,
     .        j0+deltaj,
     .        j0,
     .        deltaj,
     .        j0
      endif
      
c      print*, 0.0, j0, j0
c      print*, 0.0, 4.0d0/9.0d0*dd**2/D


      
      do i=1,n
         o=omega(i)
         o2=omega(i)*2.0d0
         ou=dd**2*o/D
         ou2=dd**2*o2/D

         call sdensdelta
     .        (nmax, nt, time, deltag, o, dt, scale, deltaj, deltaj2)
         
         oout=o
         if (domhz.eq.1) then
            oout=o*1e6
         endif
         
         if (relax.eq.1) then
            jomega=xfreed(ou,sx,D,dd,xmin0,norm0)
            jomega2=xfreed(ou2,sx,D,dd,xmin0,norm0)

            print*, oout,
     .           (jomega+deltaj +4.0d0*(jomega2+deltaj2))*fac/5.0d0,
     .           (jomega+4.0d0*jomega2)*fac/5.0d0,            
     .           (deltaj+4.0d0*deltaj2)*fac/5.0d0,
     .           j0*fac-fm*sqrt(o*1.0d12)
            
         else
            jomega=xfreed(ou,sx,D,dd,xmin0,norm0)
            print*, oout, jomega+deltaj, jomega, deltaj,
     .           j0-sqrt(2.0d0)/6.0d0*dd**3/sqrt(D**3)*sqrt(o)
            
         endif
      enddo
      

      end

      subroutine sdensdelta
     .     (nmax, nt, time,deltag, o, dt, scale, deltaj, deltaj2)

      implicit none
      integer it, nt , nmax
      double precision time(nmax), deltag(nmax), dt, o, scale
      double precision deltaj, deltaj2

      deltaj=0.0d0
      deltaj2=0.0d0
      
      do it=1,nt-1
         deltaj=deltaj +
     .        0.5d0*(deltag(it)*cos(time(it)*o) +
     .               deltag(it+1)*cos(time(it+1)*o))
     .        * (time(it+1)-time(it))
         deltaj2=deltaj2 +
     .        0.5d0*(deltag(it)*cos(2.0d0*time(it)*o) +
     .        deltag(it+1)*cos(2.0d0*time(it+1)*o) )
     .        * (time(it+1)-time(it))
      enddo

      deltaj=deltaj*scale
      deltaj2=deltaj2*scale           

      end
      

      double precision function xfreed (ou, sx, D, dd, xmin,norm)

      implicit none
      
      integer nmax, i
      double precision xmax, xmin, dx, x, x2, t, sx, ous
      double precision D, dd, u, ou, b, c, sum, pi, norm

      nmax=1000000

      pi=4.0d0*datan(1.0d0)

      ous=ou/sx**(2.0/3.0)
      
      xmax=10.0d0/sqrt(ous)
      dx=xmax/dble(nmax)
     
      sum=0.0d0
      do i=1,nmax
         x=xmin+dble(i)*dx
         x2=x*x
         b=1.0d0/(81.0d0+9.0d0*x2-2.0d0*x2**2 + x2**3)
         c=1.0d0/(1.0d0+ous**2/(x2*x2))
         
         sum=sum+b*c
      enddo
      xfreed=sum*dx/norm*dd**2/D*sx**(1.0/3.0)
      end


    
      

      logical function getline(xval)
      implicit none
      character*120 input
      double precision xval,yval

      getline=.true.
 99   continue
      read(*,'(a)',end=100, err=100) input
      if (input(1:1).eq.'#'.or.input.eq.' ') goto 99
      read(input,*,err=100,end=100) xval
      return
 100  continue
      getline=.false.
      end
