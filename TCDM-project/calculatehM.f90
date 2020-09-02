        subroutine calcuhM(hM)

        integer :: i,j
        real :: mlon,mlat
        real,dimension(nx) :: lon
        real,dimension(my) :: lat
        real,dimension(nx,my) :: hM
!================================================        
        integer :: nx,my
        real :: dx,dy

        integer :: mxc,myc
        real :: hmax,sigmax,sigmay
!================================================
        integer :: dt,time,runtime
!================================================
        common /TERRAIN2/ hmax,sigmax,sigmay,mxc,myc
        common /DOMAIN2/ nx,my,dx,dy,dt,time,runtime
!================================================
        do i = 1,nx

          lon(i) = float(i)

        enddo

        do j = 1,my

          lat(j) = float(j)

        enddo

        mlon = lon(mxc)
        mlat = lat(myc)

        do i = 1,nx
        do j = 1,my

            hM(i,j) = hmax*exp(-((float(i)-mlon)*dx)**2./(sigmax**2)&
&-((float(j)-mlat)*dy)**2/(sigmay**2))
        
        enddo
        enddo

        return
        end subroutine calcuhM
