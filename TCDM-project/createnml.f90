        !subroutine createNML()
        program createNML
        implicit none
        
        integer :: nx,my,dt,time,runtime
        real :: dx,dy

        integer :: mxc,myc
        real :: hmax,sigmax,sigmay

        real :: Depth,u0,v0

        real :: rvm,b,Vmax,maxr,dr
        integer :: xc,yc

        real :: pi,d2r,g
        real :: omega,f

        namelist /DOMAIN/ nx,my,dx,dy,dt,time,runtime
        namelist /TERRAIN/ hmax,sigmax,sigmay,mxc,myc
        namelist /ENVIRONMENT/ Depth,u0,v0
        namelist /TC/ rvm,b,Vmax,dr,xc,yc
        namelist /GENERAL/ pi,d2r,g,omega,f

        open(11,file='namelist.swm')

        nx = 401
        my = 401
        dx = 5000.
        dy = 5000.
        dt = 5
        time = 72
        runtime = 6
        
        write(11,nml=DOMAIN)

        hmax = 3000. ! m
        sigmax = 80000. ! m
        sigmay = 160000. ! m
        mxc = 151
        myc = 201

        write(11,nml=TERRAIN)

        u0 = -4.
        v0 = 0.
        Depth = 5000.

        write(11,nml=ENVIRONMENT)

        rvm = 100000. 
        b = 2.
        Vmax = 20. ! m/s
        dr = 500. ! 同心圓間距 m
        xc = 261
        yc = 201

        write(11,nml=TC)

        pi = acos(-1.0) 
        d2r = pi/180.
        g = 9.8 ! m/s
        omega = (360./86400.)*d2r
        f = 2.*omega*sin(30.*d2r)
        
        write(11,nml=GENERAL)

        close(11)

        end
        !end subroutine createNML
