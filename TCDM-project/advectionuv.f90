        subroutine advectionuv(phi,u,v)
        implicit none

        integer :: i,j
        real :: dtt
        real,dimension(nx,my) :: u,v
        real,dimension(nx,my) :: uu,vv
        real,dimension(nx,my) :: phiu,phiuv
        real,dimension(nx,my) :: phi
!================================================
        real :: dx,dy
        integer :: nx,my,dt,time,runtime
 
        common /DOMAIN2/ nx,my,dx,dy,dt,time,runtime
!================================================
        dtt = float(dt)

        uu = u
        vv = v

        do i = 1,nx
        do j = 1,my

          phiu(i,j) = 0.

        enddo
        enddo

        do j = 1,my
        do i = 2,nx-1

        if (uu(i,j).ge.0.) then
      
          if ((i.eq.2).or.(i.eq.nx-1)) then

            phiu(i,j) = phi(i,j) - dtt*uu(i,j)*((phi(i,j)-phi(i-1,j))/(dx))

          else

            phiu(i,j) = phi(i,j) - dtt*uu(i,j)*((phi(i+1,j)-phi(i-1,j))/(2.*dx) -&
&(phi(i+1,j)-3.*phi(i,j)+3.*phi(i-1,j)-phi(i-2,j))/(6.*dx))

          endif
             
        elseif (uu(i,j).lt.0.) then

          if ((i.eq.nx-1).or.(i.eq.2)) then

            phiu(i,j) = phi(i,j) - dtt*uu(i,j)*((phi(i+1,j)-phi(i,j))/(dx))

          else

            phiu(i,j) = phi(i,j) - dtt*uu(i,j)*((phi(i+1,j)-phi(i-1,j))/(2*dx) -&
&(phi(i+2,j)-3.*phi(i+1,j)+3.*phi(i,j)-phi(i-1,j))/(6.*dx))
             
          endif

        endif
        
        enddo
        enddo
        
        phiu(1,:) = phiu(2,:)
        phiu(nx,:) = phiu(nx-1,:)
!        phiu(1,2:my-1) = phiu(2,2:my-1)
!        phiu(nx,2:my-1) = phiu(nx-1,2:my-1)
!        phiu(2:nx-1,1) = phiu(2:nx-1,2)
!        phiu(2:nx-1,my) = phiu(2:nx-1,my-1)

!        phiu(1,1) = (phiu(2,1) + phiu(1,2))/2.
!        phiu(nx,1) = (phiu(nx-1,1) + phiu(nx,2))/2.
!        phiu(1,my) = (phiu(1,my-1) + phiu(2,my))/2.
!        phiu(nx,my) = (phiu(nx-1,my) + phiu(nx,my-1))/2.

        do i = 1,nx
        do j = 1,my

          phiuv(i,j) = 0.

        enddo
        enddo

        do i = 1,nx
        do j = 2,my-1

          if (vv(i,j).ge.0.) then

            if ((j.eq.2).or.(j.eq.my-1)) then

              phiuv(i,j) = phiu(i,j) - dtt*vv(i,j)*((phiu(i,j)-phiu(i,j-1))/(dy))
                
            else

              phiuv(i,j) = phiu(i,j) - dtt*vv(i,j)*((phiu(i,j+1)-phiu(i,j-1))/(2*dy) -&
&(phiu(i,j+1)-3.*phiu(i,j)+3.*phiu(i,j-1)-phiu(i,j-2))/(6.*dy))

            endif

          elseif (vv(i,j).lt.0.) then

            if ((j.eq.my-1).or.(j.eq.2)) then

              phiuv(i,j) = phiu(i,j) - dtt*vv(i,j)*((phiu(i,j+1)-phiu(i,j))/(dy))

            else

              phiuv(i,j) = phiu(i,j) - dtt*vv(i,j)*((phiu(i,j+1)-phiu(i,j-1))/(2.*dy) -&
&(phiu(i,j+2)-3.*phiu(i,j+1)+3.*phiu(i,j)-phiu(i,j-1))/(6.*dy))
                
            endif

          endif
        
        enddo
        enddo

!        phiuv(1,2:my-1) = phiuv(2,2:my-1)
!        phiuv(nx,2:my-1) = phiuv(nx-1,2:my-1)
!        phiuv(2:nx-1,1) = phiuv(2:nx-1,2)
!        phiuv(2:nx-1,my) = phiuv(2:nx-1,my-1)
        phiuv(:,1) = phiuv(:,2)
        phiuv(:,my) = phiuv(:,my-1)

        phiuv(1,1) = (phiuv(2,1) + phiuv(1,2))/2.
        phiuv(nx,1) = (phiuv(nx-1,1) + phiuv(nx,2))/2.
        phiuv(1,my) = (phiuv(1,my-1) + phiuv(2,my))/2.
        phiuv(nx,my) = (phiuv(nx-1,my) + phiuv(nx,my-1))/2.

        phi = phiuv

        return
        end subroutine advectionuv

