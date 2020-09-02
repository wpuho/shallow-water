!====================================================================== 
! Caculate Vorticity
!----------------------------------------------------------------------
      SUBROUTINE calcVorticity(U,V,vor)
      INTEGER :: NX,NY
      REAL,DIMENSION(NX,NY) :: U,V,vor
      REAL    :: DX,DY
      REAL    :: dVDX,dUDY
      COMMON /grid/ NX,NY
      COMMON /res/ DX,DY

      DO i=1,NX ; DO j=1,NY
        IF (i.EQ.1 .OR. i.EQ.NX)THEN
          IF (j.EQ.1)THEN
            dVDX = SD(V(i,j),V(i+1,j),V(i+2,j),DX, 1.)
          ELSE                              
            dVDX = SD(V(i,j),V(i-1,j),V(i-2,j),DX,-1.)
          ENDIF
        ELSE
          dVDX = CFD(V(i-1,j),V(i+1,j),DX)
        ENDIF
       ! - - - - - - - - - - - - - - 
        IF (j.EQ.1 .OR. j.EQ.NX)THEN
          IF (j.EQ.1)THEN
            dUDY = SD(U(i,j),U(i,j+1),U(i,j+2),DY, 1.)
          ELSE
            dUDY = SD(U(i,j),U(i,j-1),U(i,j-2),DY,-1.)
          ENDIF
        ELSE
          dUDY = CFD(U(i,j-1),U(i,j+1),DY)
        ENDIF
        vor(i,j) = dVDX-dUDY
      ENDDO ; ENDDO ! j,i

      RETURN
      END SUBROUTINE calcVorticity
!====================================================================== 
! Caculate Potential Vorticity
!----------------------------------------------------------------------
      SUBROUTINE calcPV(vor,cor,H,Z,PV)
      INTEGER :: NX,NY
      REAL,DIMENSION(NX,NY) :: vor,cor,H,Z,PV
      COMMON /grid/ NX,NY
      COMMON /res/ DX,DY
      
      PV = (cor+vor) / (H-Z)

      RETURN
      END SUBROUTINE calcPV
!====================================================================== 
!     Using SOR to solve the Possion equation: Laplace(str)= vor
!     DO the streamfunction: str(nx,ny)
!     given the source term (vorticity): vor
!     and lateral b.c.: -d(str)/DY= U, d(str)/DX= V
!
!     omega: 0 < omega < 2 DO stability of the iteration (under testing...)
!     stop iteration when diff is less than 0.001 or iterations > 10000
!----------------------------------------------------------------------
      SUBROUTINE calcStream(U,V,vor,str)
      INTEGER :: NX,NY
      REAL,DIMENSION(NX,NY) :: U,V,str,str2,vor,diff2d
      INTEGER :: n
      REAL    :: DX,DY,domg
      REAL    :: ratio,resij,avestr
      REAL    :: strL,strD,strR,strU
      REAL    :: omega,diff,diff2
      COMMON /grid/ NX,NY
      COMMON /res/ DX,DY
      n     = 0
      str   = 0
      diff  = 99999 ! current diff
      diff2 = 99999 ! last diff
      omega = 1.9   ! DYnamic, decrease when diff > diff2 
      domg  = 0.1
 109  FORMAT("     Iteration n=",i5," omega=",f6.3,&
                    " ave=",e10.3," diff=",e10.3)

      DO WHILE( diff.GT.0.001 .AND. n.LT.10000)
        n=n+1
        DO i=1,NX ; DO j=1,NY
          IF(i.EQ.1)THEN
            strL=str(i+1,j)-2*V(i,j)*DX   ! using lateral b.c.
          ELSE
            strL=str(i-1,j)
          ENDIF
          IF(j.EQ.1)THEN
            strD=str(i,j+1)+2*U(i,j)*DY   ! using lateral b.c.
          ELSE
            strD=str(i,j-1)
          ENDIF
          IF(i.EQ.nx)THEN
            strR=str(i-1,j)+2*V(i,j)*DX   ! using lateral b.c.
          ELSE
            strR=str(i+1,j)
          ENDIF
          IF(j.EQ.ny)THEN
            strU=str(i,j-1)-2*U(i,j)*DY   ! using lateral b.c.
          ELSE
            strU=str(i,j+1)
          ENDIF
          IF (DX.EQ.DY)THEN !  DO DX=DY 
            str2(i,j)= str(i,j)+omega*(strU+strD+strR+strL &
                      -4*str(i,j)-vor(i,j)*DX*DX)/4
          ELSE  !  DO DX ~= DY
            ratio = (DX*DX)*(DY*DY)/(2*(DX*DX)+2*(DY*DY))
            resij = (strR+strL-2*str(i,j))/(DX*DX)  &
                   +(strU+strD-2*str(i,j))/(DY*DY) - vor(i,j)
            str2(i,j) = str(i,j)+omega*ratio*resij
          ENDIF
        ENDDO ; ENDDO ! i,j
        avestr = avg2d(abs(str2),NX,NY)
        diff2d = abs((str2-str)/avestr)   ! using fractional difference
        diff   = max2d(diff2d,NX,NY)
        
        IF (diff.GT.diff2 .OR. isnan(avestr)) THEN
          ! str reserved
          omega = omega-domg   !
          IF (omega .LE.0 )THEN
            omega = omega+domg
            domg  = domg/10.
            omega = omega-domg
          ENDIF
        ELSE
          diff2 = diff         ! update diff
          str = str2           ! adjust the threshold DO diff when getting slow convergence... 
          str = str - max2d(str,NX,NY) ! extract the constant value to reduce round-off error
        ENDIF
!        IF (mod(n,10) .EQ. 0)THEN
!        WRITE(*,109) n,omega,avestr,diff
!        ENDIF
      ENDDO
      WRITE( *,109) n,omega,avestr,diff
      IF (omega.LE.0) STOP
      RETURN
      END SUBROUTINE calcStream
!======================================================================      
      REAL FUNCTION avg2d(values,NX,NY)
      INTEGER :: NX,NY
      REAL    :: values(NX,NY)
      avg2d = 0.
      DO i=1,NX ; DO j=1,NY
      avg2d = avg2d+values(i,j)
      ENDDO ; ENDDO
      avg2d = avg2d/NX/NY
      RETURN
      END FUNCTION avg2d
!======================================================================      
      REAL FUNCTION max2d(values,NX,NY)
      INTEGER :: NX,NY
      REAL    :: values(NX,NY)
      max2d = -999999.
      DO i=1,NX ; DO j=1,NY
      max2d = max(max2d,values(i,j))
      ENDDO ; ENDDO
      RETURN
      END FUNCTION max2d
!======================================================================      
      REAL FUNCTION min2d(values,NX,NY)
      INTEGER :: NX,NY
      REAL    :: values(NX,NY)
      min2d =  999999.
      DO i=1,NX ; DO j=1,NY
      min2d = min(min2d,values(i,j))
      ENDDO ; ENDDO
      RETURN
      END FUNCTION min2d
!======================================================================      
      REAL FUNCTION max2d_ind(values,NX,NY,ix,iy)
      INTEGER :: NX,NY
      INTEGER :: ix,iy
      REAL    :: values(NX,NY)
      max2d_ind = -999999.
      DO i=1,NX ; DO j=1,NY
      IF (values(i,j).GT.max2d_ind)THEN
        ix        = i
        iy        = j
        max2d_ind = values(i,j)
      ENDIF
      ENDDO ; ENDDO
      RETURN
      END FUNCTION max2d_ind
!======================================================================      
      REAL FUNCTION min2d_ind(values,NX,NY,ix,iy)
      INTEGER :: NX,NY
      INTEGER :: ix,iy
      REAL    :: values(NX,NY)
      min2d_ind = 999999.
      DO i=1,NX ; DO j=1,NY
      IF (values(i,j).LT.min2d_ind)THEN
        ix        = i
        iy        = j
        min2d_ind = values(i,j)
      ENDIF
      ENDDO ; ENDDO
      RETURN
      END FUNCTION min2d_ind
!======================================================================      
! Calculating the Radial wind and tangential wind.
!   phi = wind direction in mathmatic angle
!----------------------------------------------------------------------
      SUBROUTINE calcWTC(U,V,spd,xc,yc,WT,WR)
      INTEGER :: NX,NY
      INTEGER :: xc,yc
      REAL    :: rx,ry,theta,phi,ws
      REAL,DIMENSION(NX,NY) :: U,V,spd,WT,WR
      COMMON /grid/ NX,NY

      DO i=1,NX ; DO j=1,NY
        rx    = i-xc
        ry    = j-yc
        theta = atan2(ry,rx)
        phi   = atan2(V(i,j),U(i,j))
        WT(i,j) = SPD(i,j)*sin(phi)*cos(theta) &
                 -SPD(i,j)*cos(phi)*sin(theta)
        WR(i,j) = SPD(i,j)*cos(phi)*cos(theta) &
                 +SPD(i,j)*sin(phi)*sin(theta) 
      ENDDO ; ENDDO ! j,i
      RETURN
      END SUBROUTINE calcWTC
!======================================================================
      SUBROUTINE calcAsymWind(WTA,U,V,xc,yc,UAsym,VAsym)
      INTEGER :: NX,NY,NR
      INTEGER :: xc,yc,mr,N
      REAL,DIMENSION(NX,NY) :: U,V
      REAL,DIMENSION(NR)    :: WTA
      REAL :: dis,UAsym,VAsym
      COMMON /grid/ NX,NY
      COMMON /gAM/  NR

      mr = max2d(WTA,1,NR)
!      DO i=mr,NR
!        IF (WTA(i).LT.0.5) THEN
!          mr=i-1
!          EXIT
!        ENDIF
!      ENDDO
      N = 0
      UAsym = 0.
      VAsym = 0.

      DO ir=0,mr
       ! left and right
        DO i=xc-ir,xc+ir,max(2*ir,1)
        DO j=yc-ir,yc+ir
          dis = sqrt((i-xc)**2.+(j-xc)**2.)
          IF (r.GT.mr) CYCLE
          IF (i.LT.1 .OR. i.GT.NX .OR. j.LT.1 .OR. j.GT.NY) CYCLE
          N = N+1
          UAsym = UAsym + U(i,j)
          VAsym = VAsym + V(i,j)
        ENDDO ! j
        ENDDO ! i
       ! bottom and top
        DO j=yc-ir,yc+ir,max(2*ir,1)
        DO i=xc-ir+1,xc+ir-1
          dis = sqrt((i-xc)**2.+(j-xc)**2.)
          IF (r.GT.mr) CYCLE
          IF (i.LT.1 .OR. i.GT.NX .OR. j.LT.1 .OR. j.GT.NY) CYCLE
          N = N+1
          UAsym = UAsym + U(i,j)
          VAsym = VAsym + V(i,j)
        ENDDO ! j
        ENDDO ! i
      ENDDO ! ir
      UAsym = UAsym/N
      VAsym = VAsym/N
      RETURN
      END SUBROUTINE calcAsymWind
!======================================================================

