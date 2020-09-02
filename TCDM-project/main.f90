        program main
        implicit none
!================================================
        integer :: i,j,ierr
        real,allocatable,dimension(:,:) :: h1,h2,hM,h,hH,hhM
        real,allocatable,dimension(:,:) :: us,vs,u,v
!================================================
        integer :: nx,my,dt,time,runtime
        real :: dx,dy

        integer :: mxc,myc
        real :: hmax,sigmax,sigmay

        real :: Depth,u0,v0

        real :: rvm,b,Vmax,dr
        integer :: xc,yc

        real :: pi,d2r,g
        real :: omega,f
        
        common /DOMAIN2/ nx,my,dx,dy,dt,time,runtime
        common /TERRAIN2/ hmax,sigmax,sigmay,mxc,myc
        common /ENVIRONMENT2/ Depth,u0,v0
        common /TC2/ rvm,b,Vmax,dr,xc,yc
        common /GENERAL2/ pi,d2r,g,omega,f
        
        namelist /DOMAIN/ nx,my,dx,dy,dt,time,runtime
        namelist /TERRAIN/ hmax,sigmax,sigmay,mxc,myc
        namelist /ENVIRONMENT/ Depth,u0,v0
        namelist /TC/ rvm,b,Vmax,dr,xc,yc
        namelist /GENERAL/ pi,d2r,g,omega,f
!================================================

        !open(11,file='namelist.swm',status='old',iostat=ierr)
        open(11,file='namelist.swm',status='old')

        !if (ierr.ne.0) then

        !  call createNML()
        !  stop

        !else

          read(11,nml=DOMAIN)
          read(11,nml=TERRAIN)
          read(11,nml=ENVIRONMENT)
          read(11,nml=TC)
          read(11,nml=GENERAL)

          close(11)

        !endif
!=======================initial============================

        allocate(hM(nx,my),h1(nx,my),h2(nx,my))
        allocate(h(nx,my),hH(nx,my),hhM(nx,my))
        allocate(us(nx,my),vs(nx,my),u(nx,my),v(nx,my))
        
        do i = 1,nx
        do j = 1,my

          h1(i,j) = 0.
          h2(i,j) = 0.
          hM(i,j) = 0.
          us(i,j) = 0.
          vs(i,j) = 0.
           u(i,j) = 0.
           v(i,j) = 0.

        enddo
        enddo
        
        call calcuhM(hM)

        call calcuh1(h1,us,vs)
        
        call calcuh2(h2,u,v)

        h =  h1 + h2

        u = u + us
        v = v + vs
        
        do i = 1,nx
        do j = 1,my

          hH(i,j) = Depth + h(i,j)

        enddo
        enddo
        
!==========================================================
        hhM = hH - hM

        open(12,file='height.dat',form='unformatted',access='direct',recl=4*nx*my)
        
        print*,'min h',minval(hH)

        write(12,rec=1) hhM
        
!===================time integral==========================

        call integral(hH,u,v)
!==========================================================
        stop

        contains
          
        !  include 'createnml.f90'
          include 'calculateh1.f90'
          include 'calculateh2.f90'
          include 'calculatehM.f90'
          include 'sor2.f90'
          include 'advectionuv.f90'
          include 'integral.f90'

        end

!================================================
