	program diffu
********************************************************************
*
*       KMC of a Diffusing Particle  DP (C) 2023
*
*       Jump-length-unit: d=1
*       Time length unit: Delta t=1
*       Random orientation after each step
*       Sampling over nsample trajectories
*
********************************************************************
	
	implicit none

        integer*4 nmax
        parameter (nmax=40000)
	
	integer*4 i, is, nsteps, nsample, pbc
        double precision deltat, deltar, dself
	double precision xi0, yi0, zi0, r2i0, ri0
	double precision xi, yi, zi, r2i, ri
	double precision ux0, uy0, uz0, r30
	double precision ux, uy, uz, r3
	double precision dx, dy, dz
        double precision g, u, v, vv, p2i
 	double precision rlimit, r2limit, rmax, dd, weight
        double precision msd(0:nmax), g2(0:nmax), p2(0:nmax)
	double precision rr6(0:nmax)

********************************************************************
	

        read(*,*) deltat, rmax, rlimit, dself, nsample

	
        deltar=sqrt(dself*6.0*deltat)


	nsteps=nmax
 
        r2limit=rlimit**2

	do i=0,nsteps
           msd(i)=0.0d0
           p2(i)=0.0d0
           g2(i)=0.0d0
           rr6(i)=0.0d0
        enddo

	weight=0.0
	
        do is=1,nsample

	   xi0=0.0d0
	   yi0=0.0d0

	   call random_number(u)
	   call random_number(v)
	   zi0=rlimit + u*(rmax-rlimit)
	   vv=r2limit/zi0**2	   
	   do while (v.gt.vv)
	      call random_number(u)
	      call random_number(v)
	      zi0=rlimit + u*(rmax-rlimit)
	      vv=r2limit/zi0**2
c	      print*, v, vv, zi0
	   enddo

c	   print*, zi0
	   
	   r2i0=zi0*zi0
	   weight=weight+1.0d0/r2i0

c	   ri0=dsqrt(r2i0)
	   ri0=zi0
	   ux0=0.0d0
	   uy0=0.0d0
	   uz0=1.0d0

           xi=xi0
           yi=yi0
           zi=zi0

c	   p2(0)=p2(0)+1.0d0
c	   g2(0)=g2(0)+1.0d0/r2i0**3
c          rr6(0)=rr6(0)+1.0d0/r2i0**3

	   p2(0)=p2(0)   + 1.0/r2i0
	   g2(0)=g2(0)   + 1.0d0/r2i0
	   rr6(0)=rr6(0) + 1.0d0/r2i0	   
	   
	   do i=1,nsteps

c	      r2i=0.5*r2limit
c	      do while (r2i.lt.r2limit)
	      call random_vector(dx,dy,dz,deltar)
	      xi=xi+dx
	      yi=yi+dy
	      zi=zi+dz
	      r2i=xi*xi+yi*yi+zi*zi
c	      enddo
	      
	      if (r2i.lt.r2limit) then
		 ri=dsqrt(r2i)		 
		 xi=xi/ri*(2.0d0*rlimit-ri)
		 yi=yi/ri*(2.0d0*rlimit-ri)
		 zi=zi/ri*(2.0d0*rlimit-ri)
c		 xi=xi-dx
c		 yi=yi-dy
c		 zi=zi-dz
	      endif
	      r2i=xi*xi + yi*yi + zi*zi              
              msd(i)=msd(i)+r2i/r2i0
c	      msd(i)=msd(i)+r2i

	      ri=dsqrt(r2i)
	      ux=xi/ri
	      uy=yi/ri
	      uz=zi/ri
	      
	      p2i=1.5d0*uz*uz-0.5d0
	      p2(i)=p2(i)+p2i/(r2i0*r2i0)
 	      rr6(i)=rr6(i)+ri0/ri**3
	      g2(i)=g2(i)+p2i*ri0/(ri**3)

	   enddo	

        enddo

        open (10,file="p2_of_t.dat")
	do i=0,nsteps
	   write(10,*)  dble(i)*deltat,p2(i)/weight
        enddo
	close(10)


        open (15,file="g2_of_t.dat")
	do i=0,nsteps
	   write(15,*)  dble(i)*deltat,
     .  g2(i)/weight
        enddo
	close(15)       

	open (16,file="norm_g2_of_t.dat")
	do i=0,nsteps
	   write(16,*)  dble(i)*deltat,
     .     g2(i)/g2(0)
        enddo
	close(16)


c	open (17,file="rr6_of_t.dat")
c	do i=0,nsteps
c	   write(17,*)  dble(i)*deltat,
c     .     rr6(i)/weight
c        enddo
c	close(17)


c        open (20,file="msd_of_t.dat")
c        do i=1,nsteps
c           write(20,*)  dble(i)*deltat, msd(i)/weight
c        enddo
c        close(20)

	end	


	subroutine random_vector(dx,dy,dz,dr)

	implicit none	
	double precision dx, dy, dz, dr, u1, u2
	double precision ran1,ran2,ransq,ranh

        ransq=2.0
	do while (ransq.gt.1.0) 
	   call random_number(u1)
	   call random_number(u2)
	   ran1=1.0-2.0*u1
	   ran2=1.0-2.0*u2
	   ransq=ran1**2+ran2**2
	enddo
	ranh=2.0*sqrt(1.0-ransq)
	dx=ran1*ranh*dr
	dy=ran2*ranh*dr
	dz=(1.0-2.0*ransq)*dr
	   
	end
