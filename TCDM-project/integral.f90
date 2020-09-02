        subroutine integral(hH,u,v)

        integer :: i,j,t
        integer :: dd
        real :: dtt
        real,dimension(nx,my) :: u,v,uphi2,vphi2,u2,v2
        real,dimension(nx,my) :: hH,hphi2,h2
        real,dimension(nx,my) :: hux,huy
        real,dimension(nx,my) :: hM
!================================================
        real :: dx,dy
        integer :: nx,my,dt,time,runtime
        real :: pi,d2r,g,omega,f

        common /DOMAIN2/ nx,my,dx,dy,dt,time,runtime
        common /GENERAL2/ pi,d2r,g,omega,f
!================================================
        dd = 1
        dtt = float(dt)

        do i = 1,nx
        do j = 1,my

          hM(i,j) = 0.

        enddo
        enddo

        call calcuhM(hM)
        
        do t = 1,3600*time/dt

          uphi2 = u
          vphi2 = v

          call advectionuv(uphi2,u,v)

          call advectionuv(vphi2,u,v)
          
          do i = 2,nx-1
          do j = 2,my-1
        
            if (i.eq.2) then

              u2(i,j) = uphi2(i,j) - dtt*g*((hH(i+1,j)-hH(i,j))/(dx)) + dtt*f*vphi2(i,j)
            
            elseif (i.eq.nx-1) then

              u2(i,j) = uphi2(i,j) - dtt*g*((hH(i,j)-hH(i-1,j))/(dx)) + dtt*f*vphi2(i,j)

            else
              
              u2(i,j) = uphi2(i,j) - dtt*g*((hH(i+1,j)-hH(i-1,j))/(2.*dx)) + dtt*f*vphi2(i,j)

            end if

            if (j.eq.2) then

              v2(i,j) = vphi2(i,j) - dtt*g*((hH(i,j+1)-hH(i,j))/(dy)) - dtt*f*uphi2(i,j)
        
            elseif (j.eq.my-1) then

              v2(i,j) = vphi2(i,j) - dtt*g*((hH(i,j)-hH(i,j-1))/(dy)) - dtt*f*uphi2(i,j)

            else

              v2(i,j) = vphi2(i,j) - dtt*g*((hH(i,j+1)-hH(i,j-1))/(2.*dy)) - dtt*f*uphi2(i,j)

            end if
        
          enddo
          enddo

          u2(1,2:my-1) = u2(2,2:my-1)
          u2(nx,2:my-1) = u2(nx-1,2:my-1)
          u2(2:nx-1,1) = u2(2:nx-1,2)
          u2(2:nx-1,my) = u2(2:nx-1,my-1)

          u2(1,1) = (u2(2,1)+u2(1,2))/2.
          u2(nx,1) = (u2(nx-1,1)+u2(nx,2))/2.
          u2(1,my) = (u2(1,my-1)+u2(2,my))/2.
          u2(nx,my) = (u2(nx-1,my)+u2(nx,my-1))/2.

          v2(1,2:my-1) = v2(2,2:my-1)
          v2(nx,2:my-1) = v2(nx-1,2:my-1)
          v2(2:nx-1,1) = v2(2:nx-1,2)
          v2(2:nx-1,my) = v2(2:nx-1,my-1)

          v2(1,1) = (v2(2,1)+v2(1,2))/2.
          v2(nx,1) = (v2(nx-1,1)+v2(nx,2))/2.
          v2(1,my) = (v2(1,my-1)+v2(2,my))/2.
          v2(nx,my) = (v2(nx-1,my)+v2(nx,my-1))/2.

          hphi2 = hH

          call advectionuv(hphi2,u2,v2)

          do i = 1,nx
          do j = 1,my

            if (i.eq.1) then

              hux(i,j) = (u2(i+1,j)-u2(i,j))/(dx)

            elseif (i.eq.nx) then

              hux(i,j) = (u2(i,j)-u2(i-1,j))/(dx)

            else

              hux(i,j) = (u2(i+1,j)-u2(i-1,j))/(2.*dx)

            endif
            
            if (j.eq.1) then

              huy(i,j) = (v2(i,j+1)-v2(i,j))/(dy)

            elseif (j.eq.my) then

              huy(i,j) = (v2(i,j)-v2(i,j-1))/(dy)

            else

              huy(i,j) = (v2(i,j+1)-v2(i,j-1))/(2.*dy)

            endif
        
          enddo
          enddo

          do i = 1,nx
          do j = 1,my

            h2(i,j) = hphi2(i,j) - dtt*(hphi2(i,j)-hM(i,j))*(hux(i,j)+huy(i,j))

          enddo
          enddo

          u = u2
          v = v2

          hH = h2

          if (mod(t,3600*runtime/dt).eq.0) then

            dd = dd + 1
            hhM = hH - hM

            write(12,rec=dd) hhM

            print*,'min H',minval(hH)          
          
          end if

        enddo
       
        return
        end subroutine integral 
