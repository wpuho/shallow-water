      SUBROUTINE createENV(U,V,H,Z,cor)
      REAL :: DX,DY
      REAL,DIMENSION(NX,NY) :: U,V,H,Z,cor,h1
      REAL :: U0,V0,H0,lat0
      REAL :: omg,pi,g
      COMMON /grid/ NX,NY
      COMMON /cnst/ omg,pi,g
      COMMON /res/ DX,DY
      NAMELIST /enviroment/ U0,V0,H0,lat0

      OPEN (12,FILE="namelist.SWM",STATUS="OLD")
      READ (12,NML=enviroment)
      CLOSE(12)

      WRITE(*,*) "Creating Environment..."
      WRITE(*,*) "  U0=",U0," V0=",V0
      WRITE(*,*) "  Depth=",H0," Cor(deg)=",lat0
      U = U0
      V = V0
      H = H0
      Z = 0
      IF (lat0.LT.0)THEN ! for beta plane
        CALL CalcBetaPlane(abs(lat0),cor)
      ELSE               ! f-plane
        cor = omg*sin(lat0*pi/180.)  
      ENDIF
     !geostrophic-wind
      h1 = 0
      DO i=1,NX ; DO j=1,NY
        IF (j.EQ.1)THEN
          h1(i,j+1) = h1(i,1  ) -2*cor(i,j)*U0*DY/g /2. ! j+1=2 : half gradient
        ELSEIF (j+1.LE.NY)THEN
          h1(i,j+1) = h1(i,j-1) -2*cor(i,j)*U0*DY/g     ! j+1=3 to NY
        ENDIF
      ENDDO ; ENDDO ! j , i
!      CALL calcDepth(U,V,cor,h1)
      CALL output_Var("SWM_input_cor.bin",cor)
      CALL output_Var("SWM_input_h1.bin",h1)
      H = H+h1
      RETURN
      ENDSUBROUTINE createENV
!     -----------------------------------------------------------------      
      SUBROUTINE createMNT(NMNT,Z)
      INTEGER :: NMNT
      REAL,DIMENSION(NX,NY) :: Z
      REAL    :: omg,pi,g
      REAL    :: DX,DY
      REAL    :: rx,ry
      REAL    :: ZM(10),Sx(10),Sy(10),ang(10)
      INTEGER :: Xc(10),Yc(10)
      COMMON /grid/ NX,NY
      COMMON /cnst/ omg,pi,g
      COMMON /res/ DX,DY
      NAMELIST /Mountains/ Xc,Yc,ZM,Sx,Sy,ang

      OPEN (13,FILE="namelist.SWM",STATUS="OLD")
      READ (13,NML=Mountains)
      CLOSE(13)
      ZM = ZM*1000. ! to m
      Sx = Sx*1000. ! to m
      Sy = Sy*1000. ! to m
      ang = ang/180.*pi

      DO k=1,NMNT
        WRITE(*,*) "Creating Mountain..."
        WRITE(*,*) "  X=",Xc(k)," Y=",Yc(k)," Max Depth=",ZM(k)
        WRITE(*,*) "  SigmaX=",Sx(k)," SigmaY=",Sy(k)
        DO i=1,NX ; DO j=1,NY
          rx     = ((i-Xc(k))*cos(ang(k))-(j-Yc(k))*sin(ang(k)))*DX 
          ry     = ((i-Xc(k))*sin(ang(k))+(j-Yc(k))*cos(ang(k)))*DY
          Z(i,j) = Z(i,j) + ZM(k)*exp(-(rx/Sx(k))**2-(ry/Sy(k))**2)
        ENDDO ; ENDDO ! j, i
      ENDDO !k
      CALL output_Var("SWM_input_Z.bin",Z)
      RETURN
      ENDSUBROUTINE createMNT
!     -----------------------------------------------------------------      
      SUBROUTINE createTC(NTC,U,V,H)
      INTEGER :: NTC
      REAL,DIMENSION(NX,NY) :: U,V,H,h2
      REAL,DIMENSION(NX,NY) :: Vtx,Vty 
      REAL :: DX,DY
      REAL    :: Vt
      REAL    :: r,rx,ry,dr,theta
      REAL    :: omg,pi,g
      INTEGER :: Xc(10),Yc(10)
      REAL    :: VM(10),RM(10),ga(10)
      COMMON /grid/ NX,NY
      COMMON /cnst/ omg,pi,g
      COMMON /res/ DX,DY
      NAMELIST /Cyclone/ Xc,Yc,VM,RM,ga

      OPEN (14,FILE="namelist.SWM",STATUS="OLD")
      READ (14,NML=Cyclone)
      CLOSE(14)
      RM = RM*1000. ! to m

      h2  = 0.
      Vtx = 0.
      Vty = 0.

      DO k=1,NTC
        WRITE(*,*) "Creating Tropical Cyclone..."
        WRITE(*,*) "  X=",Xc(k)," Y=",Yc(k)," gamma=",ga(k)
        WRITE(*,*) "  Max Speed=",VM(k)," Max Radius=",RM(k)

        DO i=1,NX ; DO j=1,NY
         ! wind
          rx  = (i-Xc(k))*DX ; ry = (j-Yc(k))*DY
          r   = sqrt( rx**2 + ry**2)
          Vt  = VM(k)*r/RM(k)*exp(1/ga(k)*(1-(r/RM(k))**ga(k)))
          theta    = atan2(ry,rx)
          Vtx(i,j) = -Vt*sin(theta) ! u-component of Vt
          Vty(i,j) =  Vt*cos(theta) ! v-component of Vt
        ENDDO ; ENDDO ! j, i
      ENDDO ! k
!     -----------------------------------------------------------------      
!     Calculate h2 according to nonlinear balance.
      IF (NTC.GT.0) CALL calcDepth(Vtx,Vty,cor,h2)
!     -----------------------------------------------------------------      
      U = U+Vtx
      V = V+Vty
      H = H+h2
      CALL output_Var("SWM_input_h2.bin",h2)
      RETURN
      END SUBROUTINE createTC
!     -----------------------------------------------------------------      
!======================================================================
      SUBROUTINE CalcBetaPlane(clat,cor)
      INTEGER :: NX,NY
      REAL,DIMENSION(NX,NY) :: cor
      REAL    :: DX,DY,clat,Re,dlat
      REAL    :: omg,pi,g
      COMMON /grid/ NX,NY
      COMMON /res/ DX,DY
      COMMON /cnst/ omg,pi,g

      Re   = 6371000. ! Earth Radius, in m
      dlat = DY/Re * 180./pi

      clat = clat - (NY/2)*dlat ! offset to south side
      DO j=1,NY
        clat = clat + dlat
        cor(:,j) = omg*sin(clat*pi/180.)
      ENDDO

      CALL output_Var("SWM_input_cor.bin",cor)
      RETURN
      END SUBROUTINE CalcBetaPlane
!======================================================================
! Using SOR to solve the Poisson equation: Laplace(h)= F(x,y)
! DO the balanced height h: h(nx,ny)
! with the source term F(x,y): fxy(nx,ny)
! from the velocity u,v: U(nx,ny) and V(nx,ny)
! and lateral b.c.: dh/DX= bc1= f*V/g, dh/DY= bc2= -f*U/g
!
! omega: 0 < omega < 2 DO stability of the iteration (under testing...)
! stop iteration when diff is less than 0.001 or iterations > 10000
!
!----------------------------------------------------------------------
      SUBROUTINE calcDepth(U,V,cor,Hin) ! (h,h2,nx,ny,DX,DY,U,V,fxy,f,omega)
      INTEGER :: NX,NY
      REAL,DIMENSION(NX,NY) :: U,V,cor,Hin,H,H2,fxy,diff2d
      INTEGER :: n
      REAL    :: DX,DY
      REAL    :: ux,vx,uy,vy
      REAL    :: fvo,jco,bc1,bc2
      REAL    :: ratio,resij,aveH
      REAL    :: g,omega,domg,diff,diff2
      REAL    :: hL,hR,hU,hD
      COMMON /grid/ NX,NY
      COMMON /res/ DX,DY
 109  FORMAT("    Iteration n=",i5," omega=",f6.3,&
                    " ave=",e10.3," diff=",e10.3)
      g = 9.8
      !
      n    = 0
      diff = 99999.
      diff2= 99999.
      ! note: DO fastering the convergence, first guess of h needs not to be set to zero,
      ! othereise, to zero at the first application...
      !
      omega = 1.9
      domg  = 0.1
      H = Hin

!      IF(omega.EQ.0) omega=1.
!     -----------------------------------------------------------------
      ! compute fxy
      DO i=1,NX ; DO j=1,NY
       ! - - - - - - - - - - - - - - - - - - -
        IF (i.EQ.1 .OR. i.EQ.NX)THEN
          IF (i.EQ.1)THEN
            ux = (U(i+1,j)-U(i,j))/DX
            vx = (V(i+1,j)-V(i,j))/DX
          ELSE
            ux = (U(i,j)-U(i-1,j))/DX
            vx = (V(i,j)-V(i-1,j))/DX
          ENDIF
        ELSE
          ux = (U(i+1,j)-U(i-1,j))/(2*DX)
          vx = (V(i+1,j)-V(i-1,j))/(2*DX)
        ENDIF
       ! - - - - - - - - - - - - - - - - - - -
        IF (j.EQ.1 .OR. j.EQ.NY)THEN
          IF(j.EQ.1) THEN
            uy = (U(i,j+1)-U(i,j))/DY
            vy = (V(i,j+1)-V(i,j))/DY
          ELSE
            uy = (U(i,j)-U(i,j-1))/DY
            vy = (V(i,j)-V(i,j-1))/DY
          ENDIF
        ELSE
          uy = (U(i,j+1)-U(i,j-1))/(2*DY)
          vy = (V(i,j+1)-V(i,j-1))/(2*DY)
        ENDIF
       ! - - - - - - - - - - - - - - - - - - -
        fvo = cor(i,j)*(vx-uy)
        jco = ux*vy-uy*vx
        fxy(i,j) = (fvo+2.*jco)/g
      ENDDO ; ENDDO ! j,i
!     -----------------------------------------------------------------
      DO WHILE( diff.GT.0.001 .AND. n.LT.10000)
        n  = n+1
        DO i=1,NX ; DO j=1,NY
          bc1 =  cor(i,j)*V(i,j)/g
          bc2 = -cor(i,j)*U(i,j)/g
          IF(i.EQ.1)THEN
            hL = H(i+1,j)-2*bc1*DX   ! using lateral b.c.
          ELSE
            hL = H(i-1,j)
          ENDIF
          IF(j.EQ.1) THEN
            hD = H(i,j+1)-2*bc2*DY   ! using lateral b.c.
          ELSE
            hD = H(i,j-1)
          ENDIF
          IF(i.EQ.NX) THEN
            hR = H(i-1,j)+2*bc1*DX   ! using lateral b.c.
          ELSE
            hR = H(i+1,j)
          ENDIF
          IF(j.EQ.NY) THEN
            hU = H(i,j-1)+2*bc2*DY   ! using lateral b.c.
          ELSE
            hU = H(i,j+1)
          ENDIF
          IF (DX.EQ.DY) THEN  !DX=DY
            H2(i,j)= H(i,j) &
                    +omega*(hU+hD+hR+hL-4.*H(i,j)-fxy(i,j)*DX*DX)/4.
          ELSE ! DX ~= DY
            ratio  = (DX*DX)*(DY*DY)/(2.*(DX*DX)+2.*(DY*DY))
            resij  = (hR+hL-2.*H(i,j))/(DX*DX) &
                    +(hU+hD-2.*H(i,j))/(DY*DY) &
                    -fxy(i,j)
            H2(i,j)= H(i,j)+omega*ratio*resij
          ENDIF
        ENDDO ; ENDDO ! j, i
          aveH   = avg2d(abs(H2),NX,NY)
          diff2d = abs((H2-H)/aveH)   ! using fractional difference
          diff   = max2d(diff2d,NX,NY)

          IF (diff.GT.diff2) THEN
            omega = omega-domg
            IF (omega .LE.0 )THEN
              omega = omega+domg
              domg  = domg/10.
              omega = omega-domg
            ENDIF
          ELSE
            diff2 = diff
            H = H2 ! adjust the threshold DO diff when getting slow convergence...
            H = H - H(1,1) ! extract the constant value to reduce round-off error
          ENDIF
        IF (mod(n,100) .EQ. 0) THEN
          WRITE(*,109) n,omega,aveH,diff
        ENDIF
        !
      ENDDO ! while
      !
      !
      WRITE(*,109) n,omega,aveH,diff
      Hin = Hin + H
      RETURN
      END SUBROUTINE calcDepth
!======================================================================
