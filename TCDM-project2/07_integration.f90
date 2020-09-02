!----------------------------------------------------------------------
!  Euler   
!----------------------------------------------------------------------
      SUBROUTINE EL(U,V,H,Z,cor)
      INTEGER :: NX,NY
      REAL,DIMENSION(NX,NY) :: U,V,H,Z,cor,sH
      REAL :: DT
      COMMON /grid/ NX,NY
      COMMON /time/ DT

      CALL updateAdvWind (U,V,DT)        ! U=> U**
      CALL updateMometum (U,V,H,cor,DT)  ! U=> U^n+1

      sH = H
      CALL updateAdvDepth(U,V,H,Z,DT)    ! H=> H**
      CALL updateDepth   (U,V,H,Z,DT,sH) ! H=> H^n+1
      RETURN
      END SUBROUTINE EL
!----------------------------------------------------------------------
!  Heun    
!     h1(n+1) = h(n)+dt*F(n)
!
!     h (n+1) = h(n)+dt/2*( F(n)+F1(n+1) )
!----------------------------------------------------------------------
      SUBROUTINE HN(U,V,H,Z,cor)
      INTEGER :: NX,NY
      REAL,DIMENSION(NX,NY) :: U,V,H,Z,cor
      REAL,DIMENSION(NX,NY) :: U_tmp,V_tmp,H_tmp
      REAL,DIMENSION(NX,NY) :: U_udt,V_udt,H_udt
      REAL,DIMENSION(NX,NY) :: dU1,dV1,dH1,dU2,dV2,dH2
      REAL :: DT
      COMMON /grid/ NX,NY
      COMMON /time/ DT

      U_udt = U ; V_udt = V ; H_udt = H
      CALL updateAdvWind (U_udt,V_udt,DT)
      CALL updateMometum (U_udt,V_udt,H_udt,cor,DT)
      CALL updateAdvDepth(U_udt,V_udt,H_udt,Z ,DT)
      CALL updateDepth   (U_udt,V_udt,H_udt,Z ,DT,H_tmp)
      dU1 = U_udt-U ; dV1 = V_udt-V ; dH1 = H_udt-H

      U_tmp = U_udt ; V_tmp = V_udt ; H_tmp = H_udt
      CALL updateAdvWind (U_udt,V_udt,DT)
      CALL updateMometum (U_udt,V_udt,H_udt,cor,DT)
      CALL updateAdvDepth(U_udt,V_udt,H_udt,Z ,DT)
      CALL updateDepth   (U_udt,V_udt,H_udt,Z ,DT,H_tmp)
      dU2 = U_udt-U_tmp ; dV2 = V_udt-V_tmp ; dH2 = H_udt-H_tmp

     ! note : dT has included in F
      U = U + (dU1+dU2)/2.
      V = V + (dV1+dV2)/2.
      H = H + (dH1+dH2)/2.
      RETURN
      END SUBROUTINE HN
!----------------------------------------------------------------------
!  Runge-Kutta   
!     h2(n+0.5) = h(n) + dt/2*F (n)
!     h3(n+0.5) = h(n) + dt/2*F2(n+0.5)
!     h4(n+1)   = h(n) + dt  *F3(n+0.5)
!
!     h (n+1)   = h(n) + dt/6*(F(n)+2*F2(n+0.5)+2*F3(n+0.5)+F4(n+1))
!----------------------------------------------------------------------
      SUBROUTINE RK(U,V,H,Z,cor)
      INTEGER :: NX,NY
      REAL,DIMENSION(NX,NY) :: U,V,H,Z,cor
      REAL,DIMENSION(NX,NY) :: U_tmp,V_tmp,H_tmp
      REAL,DIMENSION(NX,NY) :: U_udt,V_udt,H_udt
      REAL,DIMENSION(NX,NY) :: dU1,dV1,dH1,dU2,dV2,dH2,&
                               dU3,dV3,dH3,dU4,dV4,dH4
      REAL :: DT
      COMMON /grid/ NX,NY
      COMMON /time/ DT

      U_udt = U     ; V_udt = V     ; H_udt = H                 ! h(n)
      U_tmp = U_udt ; V_tmp = V_udt ; H_tmp = H_udt
      CALL updateAdvWind (U_udt,V_udt,DT)
      CALL updateMometum (U_udt,V_udt,H_udt,cor,DT)
      CALL updateAdvDepth(U_udt,V_udt,H_udt, Z,DT)
      CALL updateDepth   (U_udt,V_udt,H_udt, Z,DT,H_tmp)
      dU1 = U_udt-U_tmp ; dV1 = V_udt-V_tmp ; dH1 = H_udt-H_tmp ! F(n)
      
      U_udt = U+dU1/2.; V_udt = V+dV1/2.; H_udt = H+dH1/2.      ! h2(n+0.5)
      U_tmp = U_udt ; V_tmp = V_udt ; H_tmp = H_udt
      CALL updateAdvWind (U_udt,V_udt,DT)
      CALL updateMometum (U_udt,V_udt,H_udt,cor,DT)
      CALL updateAdvDepth(U_udt,V_udt,H_udt,Z ,DT)
      CALL updateDepth   (U_udt,V_udt,H_udt,Z ,DT,H_tmp)
      dU2 = U_udt-U_tmp ; dV2 = V_udt-V_tmp ; dH2 = H_udt-H_tmp ! F2(n+0.5)

      U_udt = U+dU2/2.; V_udt = V+dV2/2.; H_udt = H+dH2/2.      ! h3(n+0.5)
      U_tmp = U_udt ; V_tmp = V_udt ; H_tmp = H_udt
      CALL updateAdvWind (U_udt,V_udt,DT)
      CALL updateMometum (U_udt,V_udt,H_udt,cor,DT)
      CALL updateAdvDepth(U_udt,V_udt,H_udt,Z ,DT)
      CALL updateDepth   (U_udt,V_udt,H_udt,Z ,DT,H_tmp)
      dU3 = U_udt-U_tmp ; dV3 = V_udt-V_tmp ; dH3 = H_udt-H_tmp ! F3(n+0.5)

      U_udt = U+dU3 ; V_udt = V+dV3 ; H_udt = H+dH3             ! h4(n+1)
      U_tmp = U_udt ; V_tmp = V_udt ; H_tmp = H_udt
      CALL updateAdvWind (U_udt,V_udt,DT)
      CALL updateMometum (U_udt,V_udt,H_udt,cor,DT)
      CALL updateAdvDepth(U_udt,V_udt,H_udt,Z ,DT)
      CALL updateDepth   (U_udt,V_udt,H_udt,Z ,DT,H_tmp)
      dU4 = U_udt-U_tmp ; dV4 = V_udt-V_tmp ; dH4 = H_udt-H_tmp ! F4(n+1)

     ! note : dT has included in F
      U = U + (dU1+2*dU2+2*dU3+dU4)/6.
      V = V + (dV1+2*dV2+2*dV3+dV4)/6.
      H = H + (dH1+2*dH2+2*dH3+dH4)/6.

      RETURN
      END SUBROUTINE RK
!----------------------------------------------------------------------
