        subroutine calcuh1(h1,us,vs)
        
        integer :: i,j
        real,dimension(nx,my) :: h1
        real,dimension(nx,my) :: us,vs
!================================================
        integer :: nx,my
        real :: dy,g,f
        real :: u0,v0 
!================================================
        integer :: dt,time,runtime
        real :: dx
        real :: Depth
        real :: pi,d2r,omega
!================================================
        common /DOMAIN2/ nx,my,dx,dy,dt,time,runtime
        common /ENVIRONMENT2/ Depth,u0,v0
        common /GENERAL2/ pi,d2r,g,omega,f
!================================================
        do i = 1,nx
        do j = 1,my

          us(i,j) = u0
          vs(i,j) = v0

        enddo
        enddo

        do i = 1,nx

          h1(i,1) = 0.

        do j = 1,my
        
          if (j.eq.1) then

            h1(i,2) = h1(i,1) - f*us(i,1)*dy/g
        
          else if(j.eq.my) then
                
            h1(i,my) = h1(i,my-1) - f*us(i,my-1)*dy/g
              
          else
            
            h1(i,j+1) = h1(i,j-1) - 2.*f*us(i,j)*dy/g
        
          end if
        
        enddo

        enddo
        
        return
        endsubroutine calcuh1
