	program diffu
********************************************************************
*
*       SIMPLE KMC of a Diffusing Particle  DP (C) 2021
*
*       Compute the <h(0)*h(t)> = <h(t)>, since h(0)=1
*       contact-pair population correlation
*       function of a 3D-Random Walker sampled over many 
*       trajectories starting at the origin
*       The contact-state h=1 is defined by minimum distance
*	r<rlimit to the origin. 
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
	double precision dx, dy, dz, box, hbox, vbox
        double precision g, u, p2i
	double precision rlimit, r2limit, dd, r2max
        double precision msd(0:nmax), g2(0:nmax), p2(0:nmax)
	double precision rr6(0:nmax)

********************************************************************
	
        nsample=1000000
        pbc=0

        read(*,*) pbc, hbox, deltat, rlimit, dself, nsample

c        deltat=0.05
c        dself=1.0/6.0
c        deltar=sqrt(0.395e-4)

c       rlimit=1.0
c        rlimit=0.9	

	r2max=hbox*hbox
	
        deltar=sqrt(dself*6.0*deltat)

c        hbox=5.0
        box=2.0*hbox

	vbox=box**3
c        vr=4.0/3.0*3.14159265*rlimit**3
c	vfac=vr/vbox

	nsteps=nmax
 
        r2limit=rlimit**2

	do i=0,nsteps
           msd(i)=0.0d0
           p2(i)=0.0d0
           g2(i)=0.0d0
           rr6(i)=0.0d0
        enddo

        do is=1,nsample

           r2i0=r2limit*0.9
           do while ((r2i0.lt.r2limit).or.(r2i0.gt.r2max))
c           do while ((r2i0.lt.r2limit))
              call random_number(u)
              xi0=box*u-hbox
              call random_number(u)
              yi0=box*u-hbox
              call random_number(u)
              zi0=box*u-hbox
              r2i0=xi0*xi0+yi0*yi0+zi0*zi0
           enddo
	   ri0=dsqrt(r2i0)
	   ux0=xi0/ri0
	   uy0=yi0/ri0
	   uz0=zi0/ri0

           xi=xi0
           yi=yi0
           zi=zi0

	   p2(0)=p2(0)+1.0d0
	   g2(0)=g2(0)+1.0d0/r2i0**3
	   rr6(0)=rr6(0)+1.0d0/r2i0**3
	   
	   do i=1,nsteps

	      call random_vector(dx,dy,dz,deltar)

	      xi=xi+dx
	      yi=yi+dy
	      zi=zi+dz
              r2i=xi*xi+yi*yi+zi*zi

	      if (r2i.lt.r2limit) then
		 ri=dsqrt(r2i)		 
		 xi=xi/ri*(2.0d0*rlimit-ri)
		 yi=yi/ri*(2.0d0*rlimit-ri)
		 zi=zi/ri*(2.0d0*rlimit-ri)
	      endif
	      
	      if (pbc.eq.1) then
              if (xi.gt.hbox) then
		 xi=xi-box
              endif
              if (xi.lt.-hbox) then
		 xi=xi+box
	      endif
              if (yi.gt.hbox) then
                yi=yi-box
              endif
              if (yi.lt.-hbox) then
		 yi=yi+box
	      endif
	      if (zi.gt.hbox) then
		 zi=zi-box
             endif
             if (zi.lt.-hbox) then
		 zi=zi+box
              endif
              endif

              r2i=xi*xi + yi*yi + zi*zi
              
              msd(i)=msd(i)+r2i

	      ri=dsqrt(r2i)
	      ux=xi/ri
	      uy=yi/ri
	      uz=zi/ri
	      
	      p2i=1.5d0*(ux0*ux+uy0*uy+uz0*uz)**2-0.5d0
	      p2(i)=p2(i)+p2i
	      rr6(i)=rr6(i)+1.0d0/(ri0**3*ri**3)
	      g2(i)=g2(i)+p2i/(ri0**3*ri**3)

	   enddo	

        enddo

        open (10,file="p2_of_t.dat")
	do i=0,nsteps
	   write(10,*)  dble(i)*deltat,p2(i)/dble(nsample)
        enddo
	close(10)


        open (15,file="g2_of_t.dat")
	do i=0,nsteps
	   write(15,*)  dble(i)*deltat,
     .  g2(i)/dble(nsample)
        enddo
	close(15)       

	open (16,file="norm_g2_of_t.dat")
	do i=0,nsteps
	   write(16,*)  dble(i)*deltat,
     .     g2(i)/g2(0)
        enddo
	close(16)


	open (17,file="rr6_of_t.dat")
	do i=0,nsteps
	   write(17,*)  dble(i)*deltat,
     .     rr6(i)/dble(nsample)
        enddo
	close(17)


        open (20,file="msd_of_t.dat")
        do i=1,nsteps
           write(20,*)  dble(i)*deltat, msd(i)/dble(nsample)
        enddo
        close(20)

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
