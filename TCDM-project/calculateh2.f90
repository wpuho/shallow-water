        subroutine calcuh2(h2,u,v)
        implicit none

        real :: rr,Vtt
        real :: rx,ry,d,theta
        real,allocatable :: r(:),Vt(:)
        real,dimension(nx,my) :: u,v,h2
        integer :: n,mm,nn
        integer :: i,j
!================================================
        real :: rvm,b,Vmax
        integer :: xc,yc
        integer :: nx,my
!================================================
        real :: dr
        real :: dx,dy
        integer :: dt,time,runtime
!===============================================
        common /TC2/ rvm,b,Vmax,dr,xc,yc
        common /DOMAIN2/ nx,my,dx,dy,dt,time,runtime
!================================================
        rr = 0.
        Vtt = 1.
        n = 1
        
        do while (Vtt>=0.0001)

          n = n + 1
          rr =  rr + dr
          Vtt = Vmax*(rr/rvm)*exp((1-(rr/rvm)**b)/b)

        enddo

        allocate(r(n))
        allocate(Vt(n))

        r(1) = 0.
        Vt(1) = 1.
        n = 1
        
        do while (Vt(n)>=0.0001)

          n = n + 1
          r(n) =  r(n-1) + dr
          Vt(n) = Vmax*(r(n)/rvm)*exp((1-(r(n)/rvm)**b)/b)

        enddo

        Vt(1) = 0.
        mm = n
!================================================
        do 100 i = 1,nx
        do 100 j = 1,my

          rx = float(i-xc)*dx
          ry = float(j-yc)*dy
          d = (rx**2 + ry**2)**0.5

         do 200 nn = 1,mm-1
           
            if ((r(nn+1).gt.d).and.(r(nn).lt.d)) then

               Vtt = (Vt(nn+1)-Vt(nn))*((d-r(nn))/(r(nn+1)-r(nn))) + Vt(nn)
               theta = atan2(ry,rx)
               u(i,j) = -1*Vtt*sin(theta)
               v(i,j) = Vtt*cos(theta)

            else if (r(nn).eq.d) then

               Vtt = Vt(nn)
               theta = atan2(ry,rx)
               u(i,j) = -1*Vtt*sin(theta)
               v(i,j) = Vtt*cos(theta)

            else if (r(mm).eq.d) then

               Vtt = Vt(mm)
               theta = atan2(ry,rx)
               u(i,j) = -1*Vtt*sin(theta)
               v(i,j) = Vtt*cos(theta)

            else if (d.gt.r(mm)) then

                Vtt = Vmax*(d/rvm)*exp((1-(d/rvm)**b)/b)
                theta = atan2(ry,rx)
                u(i,j) = -1*Vtt*sin(theta)
                v(i,j) = Vtt*cos(theta)

            end if
           
          200 enddo
       
        100 enddo
!================================================
        do i = 1,nx
        do j = 1,my

          h2(i,j) = 0.

        enddo
        enddo
        
        call sor2(h2,u,v,1.)
        
        return
        end subroutine calcuh2

