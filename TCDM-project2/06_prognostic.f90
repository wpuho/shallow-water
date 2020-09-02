!     =================================================================
!     Update advection term of mometum.
!     =================================================================
      SUBROUTINE updateAdvWind(U,V,DT)
      INTEGER :: NX,NY
      REAL,DIMENSION(NX,NY) :: U,V,dU,dV
      REAL,DIMENSION(NX,NY) :: sU,sV
      REAL :: DX,DY,DT
      COMMON /grid/ NX,NY
      COMMON /res/ DX,DY

      sU = U
      sV = V
!     U adv of U- - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO i=1,NX ; DO j=1,NY
      IF (U(i,j).GE.0)THEN ! U>0
        IF (i-2.LT.1 .OR. i+1.GT.NX)THEN
          dU(i,j) = 0 ! No change
        ELSE
          dU(i,j) = -sU(i,j)*advCFD(U(i-1,j),U(i+1,j),U(i-2:i+1,j),DX)
        ENDIF
      ELSE ! U<0
        IF (i-1.LT.1 .OR. i+2.GT.NX)THEN
          dU(i,j) = 0 ! No change
        ELSE
          dU(i,j) = -sU(i,j)*advCFD(U(i-1,j),U(i+1,j),U(i-1:i+2,j),DX)
        ENDIF
      ENDIF ! U
      ENDDO ; ENDDO ! j,i
      U = U+DT*dU
      CALL ZeroGrad(U,.True.) ! replace outtest grids as zero gradient.
!     U adv of V- - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO i=1,NX ; DO j=1,NY
      IF (U(i,j).GE.0)THEN ! U>0
        IF (i-2.LT.1 .OR. i+1.GT.NX)THEN
          dV(i,j) = 0 ! No change
        ELSE
          dV(i,j) = -sU(i,j)*advCFD(V(i-1,j),V(i+1,j),V(i-2:i+1,j),DX)
        ENDIF
      ELSE ! U<0
        IF (i-1.LT.1 .OR. i+2.GT.NX)THEN
          dV(i,j) = 0 ! No change
        ELSE
          dV(i,j) = -sU(i,j)*advCFD(V(i-1,j),V(i+1,j),V(i-1:i+2,j),DX)
        ENDIF
      ENDIF ! U
      ENDDO ; ENDDO ! j,i
      V = V+DT*dV
      CALL ZeroGrad(U,.True.) ! replace outtest grids as zero gradient.
!     V adv of U- - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO i=1,NX ; DO j=1,NY
      IF (V(i,j).GE.0)THEN ! V>0
        IF (j-2.LT.1 .OR. j+1.GT.NY)THEN
          dU(i,j) = 0 ! No change 
        ELSE
          dU(i,j) = -sV(i,j)*advCFD(U(i,j-1),U(i,j+1),U(i,j-2:j+1),DY)
        ENDIF
      ELSE ! V<0 
        IF (j-1.LT.1 .OR. j+2.GT.NY)THEN
          dU(i,j) = 0 ! No change
        ELSE
          dU(i,j) = -sV(i,j)*advCFD(U(i,j-1),U(i,j+1),U(i,j-1:j+2),DY)
        ENDIF
      ENDIF
      ENDDO ; ENDDO ! j,i
      U = U+DT*dU
      CALL ZeroGrad(U,.True.) ! replace outtest grids as zero gradient.
!     V adv of V- - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO i=1,NX ; DO j=1,NY
      IF (V(i,j).GE.0)THEN ! V>0
        IF (j-2.LT.1 .OR. j+1.GT.NY)THEN
          dV(i,j) = 0 ! No change 
        ELSE
          dV(i,j) = -sV(i,j)*advCFD(V(i,j-1),V(i,j+1),V(i,j-2:j+1),DY)
        ENDIF
      ELSE ! V<0 
        IF (j-1.LT.1 .OR. j+2.GT.NY)THEN
          dV(i,j) = 0 ! No change
        ELSE
          dV(i,j) = -sV(i,j)*advCFD(V(i,j-1),V(i,j+1),V(i,j-1:j+2),DY)
        ENDIF
      ENDIF
      ENDDO ; ENDDO ! j,i
      V = V+DT*dV
      CALL ZeroGrad(V,.True.) ! replace outtest grids as zero gradient.
!     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      RETURN
      END SUBROUTINE updateAdvWind
!     =================================================================
!     Update pressure gradient forece and Coriolis force term of mometum.
!     =================================================================
      SUBROUTINE updateMometum(U,V,H,cor,DT)
      INTEGER :: NX,NY
      REAL,DIMENSION(NX,NY) :: U,V,H,cor
      REAL,DIMENSION(NX,NY) :: dU,dV
      REAL    :: dHdX,dHdY
      REAL    :: omg,pi,g
      REAL :: DX,DY,DT
      COMMON /grid/ NX,NY
      COMMON /res/ DX,DY
      COMMON /cnst/ omg,pi,g

!     geostrophic balance---------------------------------------------
      DO i=1,NX ; DO j=1,NY
        IF (i.EQ.1 .OR. i.EQ.NX)THEN !!!!!!!
          dU(i,j) = 0. ! geostrophic balance
        ELSE
          dHdX    = CFD(H(i-1,j),H(i+1,j),DX)
          dU(i,j) = -g*dHdX + cor(i,j)*V(i,j)
        ENDIF
       ! - - - - - - - - - - - - - - - - - - - - 
        IF (j.EQ.1 .OR. j.EQ.NY)THEN !!!!!
          dV(i,j) = 0. ! geostrophic balance
        ELSE
          dHdY    = CFD(H(i,j-1),H(i,j+1),DY)
          dV(i,j) = -g*dHdY - cor(i,j)*U(i,j) 
        ENDIF
      ENDDO ; ENDDO ! j,i
!     One side difference----------------------------------------------
!      DO i=1,NX ; DO j=1,NY
!        IF (i.EQ.1 .OR. i.EQ.NX)THEN !!!!!!!
!          IF(i.EQ.1) THEN
!            dHdX = SD(H(i,j),H(i+1,j),H(i+2,j),DX, 1.)
!          ELSE
!            dHdX = SD(H(i,j),H(i-1,j),H(i-2,j),DX,-1.)
!          ENDIF
!        ELSE
!          dHdX = CFD(H(i-1,j),H(i+1,j),DX)
!        ENDIF
!       ! - - - - - - - - - - - - - - - - - - - -
!        IF (j.EQ.1 .OR. j.EQ.NY)THEN !!!!!
!          IF(j.EQ.1) THEN
!            dHdY = SD(H(i,j),H(i,j+1),H(i,j+2),DY, 1.)
!          ELSE
!            dHdY = SD(H(i,j),H(i,j-1),H(i,j-2),DY,-1.)
!          ENDIF
!        ELSE
!          dHdY = CFD(H(i,j-1),H(i,j+1),DY)
!        ENDIF
!       ! - - - - - - - - - - - - - - - - - - - -
!        dU(i,j) = -g*dHdX + cor(i,j)*V(i,j)
!        dV(i,j) = -g*dHdY - cor(i,j)*U(i,j)
!      ENDDO ; ENDDO ! j,i
!     -----------------------------------------------------------------
      U = U+DT*dU
      V = V+DT*dV
      CALL ZeroGrad(U,.True.) ! replace outtest grids as zero gradient.
      CALL ZeroGrad(V,.True.) ! replace outtest grids as zero gradient.
      RETURN
      END SUBROUTINE updateMometum
!     =================================================================
!     Update advection term of depth(H).
!     =================================================================
      SUBROUTINE updateAdvDepth(U,V,H,Z,DT)
      INTEGER :: NX,NY
      REAL,DIMENSION(NX,NY) :: U,V,H,Z
      REAL,DIMENSION(NX,NY) :: dH,HH
      REAL :: DX,DY,DT
      COMMON /grid/ NX,NY
      COMMON /res/ DX,DY

!     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      HH = H-Z !!!!!!!!!!!!!!!!!!
      DO i=1,NX ; DO j=1,NY
      IF (U(i,j).GE.0)THEN ! U>=0
        IF (i-2.LT.1 .OR. i+1.GT.NX)THEN
          dH(i,j) = 0 ! No change
        ELSE
          dH(i,j) = -U(i,j)*advCFD(HH(i-1,j),HH(i+1,j),HH(i-2:i+1,j),DX)
        ENDIF
      ELSE ! U<0
        IF (i-1.LT.1 .OR. i+2.GT.NX)THEN
          dH(i,j) = 0 ! No change
        ELSE
          dH(i,j) = -U(i,j)*advCFD(HH(i-1,j),HH(i+1,j),HH(i-1:i+2,j),DX)
        ENDIF
      ENDIF
      ENDDO ; ENDDO ! j,i
      H = H+DT*dH
      CALL ZeroGrad(H,.False.) ! Replace outtest grids as zero gradient except the south & north side.
                               ! Remaining there have pressure gradient.
!     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      HH = H-Z !!!!!!!!!!!!!!!!!!
      DO i=1,NX ; DO j=1,NY
        IF (V(i,j).GE.0)THEN !V>=0
          IF (j-2.LT.1 .OR. j+1.GT.NY)THEN
            dH(i,j)= 0 ! No change
          ELSE
            dH(i,j)=-V(i,j)*advCFD(HH(i,j-1),HH(i,j+1),HH(i,j-2:j+1),DY)
          ENDIF
        ELSE ! V<0 
          IF (j-1.LT.1 .OR. j+2.GT.NY)THEN
            dH(i,j)= 0 ! No change 
          ELSE
            dH(i,j)=-V(i,j)*advCFD(HH(i,j-1),HH(i,j+1),HH(i,j-1:j+2),DY)
          ENDIF
        ENDIF
      ENDDO ; ENDDO ! j,i
      H = H+DT*dH
      CALL ZeroGrad(H,.False.) ! Replace outtest grids as zero gradient except the south & north side.
                               ! Remaining there have pressure gradient.
!     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      RETURN
      END SUBROUTINE updateAdvDepth
!     =================================================================
!     Update divergence term of depth(H).
!        sH is not change in subroutine.
!     =================================================================
      SUBROUTINE updateDepth(U,V,H,Z,DT,sH)
      INTEGER :: NX,NY
      REAL,DIMENSION(NX,NY) :: U,V,H,Z
      REAL,DIMENSION(NX,NY) :: dH,HH,sH
      REAL    :: dUdX,dVdY
      REAL :: DX,DY,DT
      COMMON /grid/ NX,NY
      COMMON /res/ DX,DY

      HH = sH-Z !!!!!!!!
!     geostrophic balance----------------------------------------------      
      DO i=1,NX ; DO j=1,NY
        IF (i.EQ.1 .OR. i.EQ.NX .OR. j.EQ.1 .OR. j.EQ.NY)THEN
          dH(i,j) = 0. ! geostrophic balance
        ELSE
          dUdX = CFD(U(i-1,j),U(i+1,j),DX)
          dVdY = CFD(V(i,j-1),V(i,j+1),DY)
          dH(i,j) = -HH(i,j)*(dUdX+dVdY)
        ENDIF
      ENDDO ; ENDDO ! j,i
!     One side difference----------------------------------------------      
!      DO i=1,NX ; DO j=1,NY
!       ! - - - - - - - - - - - - - - -
!        IF (i.EQ.1 .OR. i.EQ.NX)THEN
!          IF (i.EQ.1)THEN
!            dUdX = SD(U(i,j),U(i+1,j),U(i+2,j),DX, 1.)
!          ELSE
!            dUdX = SD(U(i,j),U(i-1,j),U(i-2,j),DX,-1.)
!          ENDIF
!        ELSE
!          dUdX = CFD(U(i-1,j),U(i+1,j),DX)
!        ENDIF
!       ! - - - - - - - - - - - - - - -
!        IF (j.EQ.1 .OR. j.EQ.NY)THEN
!          IF (j.EQ.1)THEN
!            dVdY = SD(V(i,j),V(i,j+1),V(i,j+2),DY, 1.)
!          ELSE
!            dVdY = SD(V(i,j),V(i,j-1),V(i,j-2),DY,-1.)
!          ENDIF
!        ELSE
!          dVdY = CFD(V(i,j-1),V(i,j+1),DY)
!        ENDIF
!       ! - - - - - - - - - - - - - - -
!        dH(i,j) = -HH(i,j)*(dUdX+dVdY)
!      ENDDO ; ENDDO ! j,i
!     -----------------------------------------------------------------      
      H = H+DT*dH
      CALL ZeroGrad(H,.False.) ! Replace outtest grids as zero gradient except the south & north side.
                               ! Remaining there have pressure gradient.
      RETURN
      END SUBROUTINE updateDepth
!     =================================================================
!     Replace outtest grids as zero gradient.
!     =================================================================
      SUBROUTINE ZeroGrad(Var,opt)
      INTEGER :: NX,NY
      LOGICAL :: opt
      REAL,DIMENSION(NX,NY) :: Var
      COMMON /grid/ NX,NY

      Var( 1,:) = Var(2   ,:)
      Var(NX,:) = Var(NX-1,:)

      IF (opt)THEN 
      Var(:, 1) = Var(:,2   )
      Var(:,NY) = Var(:,NY-1)
      ENDIF

      RETURN
      END SUBROUTINE ZeroGrad
!     =================================================================
