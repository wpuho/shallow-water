      SUBROUTINE shapiro(Var,mu,times)
      INTEGER :: NX,NY
      INTEGER :: times
      REAL,DIMENSION(NX,NY) :: Var,oVar
      REAL    :: mu(times)
      COMMON /grid/ NX,NY

      DO k=1,times
        oVar = Var
        DO i=2,NX-1 ; DO j=2,NY-1
          Var(i,j) = (1-mu(k))*oVar(i,j) &
                    +   mu(k) *( oVar(i+1,j)+oVar(i-1,j) &
                                +oVar(i,j+1)+oVar(i,j-1))/4.
        ENDDO ; ENDDO ! j,i
      ENDDO !k

      RETURN
      END SUBROUTINE shapiro

