      program freedhwang
c
c
c

      implicit none

      integer nmax, nargmax
      parameter (nmax=100000, nargmax=100)
      
      integer i, it, n, nt, narg, ntmax
      integer relax, domhz
      
      double precision Ak, tauk, betak
      double precision pi, dt, rhh
      double precision o, oout, u, omegau, t, tmin, tmax, norm0, norm
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
         if (arg(i).eq.'-A') then
            read(arg(i+1),*) Ak
         endif
         if (arg(i).eq.'-tau') then
            read(arg(i+1),*) tauk
         endif
         if (arg(i).eq.'-beta') then
            read(arg(i+1),*) betak
         endif     
         if (arg(i).eq.'-rhh') then
            read(arg(i+1),*) rhh
         endif  
         if (arg(i).eq.'-relax') then
            relax=1
         endif
         if (arg(i).eq.'-mhz') then
            domhz=1
         endif
         if (arg(i).eq.'-tmax') then
            read(arg(i+1),*) tmax
         endif          
      enddo

      n=0
      do while(getline(omega(n+1)).and.n.lt.nmax)
         if (domhz.eq.1) then
            omega(n+1)=omega(n+1)*1.0d-6
         endif
         n=n+1
      enddo

      
      pi=4.0d0*datan(1.0d0)



      fac=2.0d0
      fac=fac*(267.5153151e6)**4
      fac=fac*(6.62607015e-34/(2*pi))**2
      fac=fac*1.0d0/2.0d0*(3.0d0/2.0d0)
      fac=fac*(1.25663706212e-6/(4.0d0*pi))**2
      fac=fac*1.0d0/((rhh*1.0d-9)**6)
      fac=fac*1.0e-12



      dt=time(2)-time(1)
      ntmax=int(tmax/dt)

      j0=0.5
      do i=1,ntmax-1
         t=dble(i)*dt
         j0=j0+Ak*dexp(-(t/tauk)**betak)
      enddo
      j0=j0*dt

      deltaj=0.0d0
      deltaj=deltaj+0.5*deltag(1)
      deltaj=deltaj+0.5*deltag(nt)
      do it=2,nt-1
         deltaj=deltaj+deltag(it)
      enddo
      deltaj=deltaj*dt

      if (relax.eq.1) then
         print*, 0.0 ,
     .        (j0+deltaj)*fac,
     .        j0*fac,
     .        deltaj*fac
      else
         print*, 0.0 , j0+deltaj, j0, deltaj
      endif

      
      do i=1,n
         o=omega(i)
         o2=omega(i)*2.0d0

         call sdensdelta
     .        (nmax, nt, time, deltag, o, dt, deltaj, deltaj2)

         call sdenskohlrausch
     .        (ntmax, dt, o, ak, tauk, betak, jomega, jomega2)
         
         oout=o
         if (domhz.eq.1) then
            oout=o*1e6
         endif

         
         if (relax.eq.1) then
            print*, oout,
     .           (jomega+deltaj +4.0d0*(jomega2+deltaj2))*fac/5.0d0,
     .           (jomega+4.0d0*jomega2)*fac/5.0d0,            
     .           (deltaj+4.0d0*deltaj2)*fac/5.0d0    
            
         else
            print*, oout, jomega+deltaj, jomega, deltaj
         endif
      enddo
      
      end


      subroutine sdenskohlrausch
     .     (ntmax, dt, o, ak, tauk, betak, jomega, jomega2)

      implicit none

      integer ntmax, i
      double precision dt, o, ak, tauk, betak, jomega, jomega2
      double precision t

      jomega=0.5
      jomega2=0.5
      do i=1,ntmax-1
         t=dble(i)*dt
         jomega=jomega+Ak*dexp(-(t/tauk)**betak)*dcos(o*t)
         jomega2=jomega2+Ak*dexp(-(t/tauk)**betak)*dcos(2.0d0*o*t)
      enddo
      jomega=jomega*dt
      jomega2=jomega2*dt
      
      end
      
      subroutine sdensdelta
     .     (nmax, nt, time, deltag, o, dt, deltaj, deltaj2)

      implicit none
      integer it, nt , nmax
      double precision time(nmax), deltag(nmax), dt, o
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
      
      deltaj=deltaj*dt
      deltaj2=deltaj2*dt     

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
