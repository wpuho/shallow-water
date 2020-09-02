        subroutine sor2(h,u,v,omg)
        ! Using SOR to solve the Poisson equation: Laplace(h)= F(x,y)
        ! for the balanced height h: h(nx,ny)
        ! with the source term F(x,y): fxy(nx,ny) 
        ! from the velocity u,v: uu(nx,ny) and vv(nx,ny)
        ! and lateral b.c.: dh/dx= bc1= f*vv/g, dh/dy= bc2= -f*uu/g 
        !
        ! omega: 0 < omega < 2 for stability of the iteration (under testing...)
        ! stop iteration when diff is less than 0.001 or iterations > 10000
        implicit none

!        integer :: nx,my

        integer :: i,j,n
        real :: omg
        real,dimension(nx,my) :: h
        real,dimension(nx,my) :: h2
        real,dimension(nx,my) :: u,v,uu,vv

        real :: ux,uy,vx,vy
        real :: fvo,jco
        real,dimension(nx,my) :: fxy

        real :: bc1,bc2
        real :: hL,hR,hU,hD

        real :: ratio,resij
        real :: aveh,diff2(nx,my),diff
!================================================
        integer :: nx,my
        real :: dx,dy
        real :: g,f
!================================================
        real :: pi,d2r,omega
        integer :: dt,time,runtime
!================================================
        common /DOMAIN2/ nx,my,dx,dy,dt,time,runtime
        common /GENERAL2/ pi,d2r,g,omega,f
!================================================        
!        g = 9.8
        n = 0
        diff = 100.
        ! note: for fastering the convergence, first guess of h needs not to be set to zero,
        ! othereise, to zero at the first application...
        
        uu = u
        vv = v

        if(omg.eq.0.) then

          omg = 1.

        endif

        do while ((diff.gt.0.001).and.(n.lt.10000))

          n = n + 1
        
        ! compute fxy
          do j=1,my
          do i=1,nx
    
            if (i.eq.1) then

              ux = (uu(i+1,j)-uu(i,j))/dx
              vx = (vv(i+1,j)-vv(i,j))/dx

            else if (i.eq.nx) then

              ux = (uu(i,j)-uu(i-1,j))/dx
              vx = (vv(i,j)-vv(i-1,j))/dx

            else

              ux = (uu(i+1,j)-uu(i-1,j))/(2.*dx)
              vx = (vv(i+1,j)-vv(i-1,j))/(2.*dx)
                
            end if

            if (j.eq.1) then

              uy = (uu(i,j+1)-uu(i,j))/dy
              vy = (vv(i,j+1)-vv(i,j))/dy

            else if (j.eq.my) then

              uy = (uu(i,j)-uu(i,j-1))/dy
              vy = (vv(i,j)-vv(i,j-1))/dy
        
            else
                    
              uy = (uu(i,j+1)-uu(i,j-1))/(2.*dy)
              vy = (vv(i,j+1)-vv(i,j-1))/(2.*dy)

            end if

          fvo = f*(vx-uy)
          jco = ux*vy-uy*vx
          fxy(i,j) = (fvo+2.*jco)/g
        
          enddo
          enddo

          do j = 1,my
          do i = 1,nx

            bc1 = f*vv(i,j)/g
            bc2 = -f*uu(i,j)/g
            if (i.eq.1) then

              hL = h(i+1,j) - 2.*bc1*dx ! using lateral b.c.

            else

              hL = h(i-1,j)
       
            end if

            if (j.eq.1) then

              hD = h(i,j+1) - 2.*bc2*dy ! using lateral b.c.

            else

              hD = h(i,j-1)

            end if

            if (i.eq.nx) then

              hR = h(i-1,j) + 2.*bc1*dx ! using lateral b.c.

            else

              hR = h(i+1,j)

            end if

            if (j.eq.my) then

              hU = h(i,j-1) + 2.*bc2*dy ! using lateral b.c.

            else

              hU = h(i,j+1)

            end if

            if (dx.eq.dy) then

              h2(i,j) = h(i,j) + omg*(hU+hD+hR+hL-4.*h(i,j)-fxy(i,j)*dx*dx)/4.

            else

              ratio = (dx*dx)*(dy*dy)/(2.*(dx*dx)+2.*(dy*dy))
              resij = (hR+hL-2*h(i,j))/(dx*dx) + (hU+hD-2.*h(i,j))/(dy*dy) - fxy(i,j)
              h2(i,j) = h(i,j) + omg*ratio*resij

            end if
       
          enddo
          enddo


          aveh = sum(abs(h2))/size(h2)
          diff2 = abs((h2-h)/aveh)   ! using fractional difference
          diff = maxval(diff2)

          ! adjust the threshold for diff when getting slow convergence...
        
          h = h2    

          ! extract the constant value to reduce round-off error
  
          h = h - h(1,1)
          u = uu
          v = vv
        
          if (mod(n,100).eq.0) then
        
            print*,'Iteration n = ',n
            print*,'Iteration diff = ',diff
        
          end if

        enddo  

        return
        end subroutine sor2
