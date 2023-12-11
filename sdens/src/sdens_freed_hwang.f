      program freedhwang
c
c
c

      implicit none

      integer nmax, nargmax
      parameter (nmax=100000, nargmax=10)
      
      integer i, it, n, nt, narg
      integer relax, domhz
      
      double precision dd, D, R, b, xmin0, xmin 
      double precision pi, freed, dt
      double precision o, oout, ou, u, omegau, t, tmin, tmax, norm0, norm
      double precision scale, jomega, jomega2, o2, ou2, fac
      double precision j0, fm, tauadd, deltaj, deltaj2, rho
      
      double precision omega(nmax), time(nmax), deltag(nmax)
      character*256  arg(nargmax)
      logical getline
c-----

      rho=1.0d0
      scale=1.0d0
      relax=0
      domhz=0
      
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
      fac=fac*(267.5153151e6)**4
      fac=fac*(6.62607015e-34/(2*pi))**2
      fac=fac*1.0d0/2.0d0*(3.0d0/2.0d0)
      fac=fac*(1.25663706212e-6/(4.0d0*pi))**2
c     fac=fac*2*4096/(4.96735e-9)**3
      fac=fac*rho*(1.0e9)**3
      fac=fac*4.0d0*pi/(3.0d0*(dd*1e-9)**3)
      fac=fac*1.0e-12


      fm=pi*sqrt(2.0d0)/15.0d0
      fm=fm*(267.5153151e6)**4
      fm=fm*(6.62607015e-34/(2*pi))**2
      fm=fm*(1.25663706212e-6/(4.0d0*pi))**2
      fm=fm*(1+4.0d0*dsqrt(2.0d0))
c     fm=fm*2*4096/(4.96735e-9)**3
      fm=fm*rho*(1.0e9)**3
      fm=fm/dsqrt((D*1.0e-6)**3)


      dt=time(2)-time(1)
      j0=(4.0d0/9.0d0*dd**2/D)
     
      deltaj=0.0d0
      deltaj=deltaj+0.5*deltag(1)
      deltaj=deltaj+0.5*deltag(nt)
      do it=2,nt-1
         deltaj=deltaj+deltag(it)
      enddo
      deltaj=deltaj*dt*scale

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

c         deltaj=0.0d0
c         deltaj=deltaj+0.5*(deltag(1)*cos(time(1)*o))
c         deltaj=deltaj+0.5*(deltag(nt)*cos(time(nt)*o))

c         deltaj2=0.0d0
c         deltaj2=deltaj2+0.5*(deltag(1)*cos(2.0d0*time(1)*o))
c         deltaj2=deltaj2+0.5*(deltag(nt)*cos(2.0d0*time(nt)*o))         

c         do it=2,nt-1
c            deltaj=deltaj+deltag(it)*cos(time(it)*o)
c            deltaj2=deltaj2+deltag(it)*cos(2.0d0*time(it)*o)
c         enddo

c         deltaj=deltaj*dt*scale
c         deltaj2=deltaj2*dt*scale     


         call sdensdelta
     .        (nmax, nt, time, deltag, o, dt, scale, deltaj, deltaj2)
         
         oout=o
         if (domhz.eq.1) then
            oout=o*1e6
         endif
         
         if (relax.eq.1) then
            jomega=freed(ou,D,dd,xmin0,norm0)
            jomega2=freed(ou2,D,dd,xmin0,norm0)

            print*, oout,
     .           (jomega+deltaj +4.0d0*(jomega2+deltaj2))*fac/5.0d0,
     .           (jomega+4.0d0*jomega2)*fac/5.0d0,            
     .           (deltaj+4.0d0*deltaj2)*fac/5.0d0,
     .           j0*fac-fm*sqrt(o*1.0d12)
            
         else
            jomega=freed(ou,D,dd,xmin0,norm0)
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
      deltaj=deltaj+0.5*(deltag(1)*cos(time(1)*o))
      deltaj=deltaj+0.5*(deltag(nt)*cos(time(nt)*o))

      deltaj2=0.0d0
      deltaj2=deltaj2+0.5*(deltag(1)*cos(2.0d0*time(1)*o))
      deltaj2=deltaj2+0.5*(deltag(nt)*cos(2.0d0*time(nt)*o))         

      do it=2,nt-1
         deltaj=deltaj+deltag(it)*cos(time(it)*o)
         deltaj2=deltaj2+deltag(it)*cos(2.0d0*time(it)*o)
      enddo
      
      deltaj=deltaj*dt*scale
      deltaj2=deltaj2*dt*scale     

      end
      

      double precision function freed (ou, D, dd, xmin,norm)

      implicit none
      
      integer nmax, i
      double precision xmax, xmin, dx, x, x2, t
      double precision D, dd, u, ou, b, c, sum, pi, norm

      nmax=1000000

      pi=4.0d0*datan(1.0d0)

c      u=1.0d0/ou
      
      xmax=10.0d0/sqrt(ou)
      dx=xmax/dble(nmax)

     
      sum=0.0d0
c      print*,'#', t, xmax
      do i=1,nmax
         x=xmin+dble(i)*dx
         x2=x*x
         b=1.0d0/(81.0d0+9.0d0*x2-2.0d0*x2**2 + x2**3)
         c=1.0d0/(1.0d0+ou**2/(x2*x2))
         
         sum=sum+b*c
c         print*,'f ', x, b*dexp(-x2*u), sum
      enddo
c      freed=sum*dx*18.0/(dd**3)*
c     freed=sum*dx/(3.1415926/54.0)
      freed=sum*dx/norm*dd**2/D
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
