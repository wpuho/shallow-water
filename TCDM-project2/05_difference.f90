!======================================================================
!     Central finite difference:
!       next F (FN) minus last F (FL).
!----------------------------------------------------------------------
      REAL FUNCTION CFD(FL,FN,DS)
      REAL    :: FL,FN 
      REAL    :: DS
      CFD = (FN-FL)/2./DS
      RETURN
      END FUNCTION CFD
!======================================================================
!     Central finite difference for advection:
!       F1 is i-1
!       F2 is i+1
!       FF(1) is the left/bottom side.
!----------------------------------------------------------------------
      REAL FUNCTION advCFD(F1,F2,FF,DS)
      REAL    :: F1,F2,FF(4)
      REAL    :: DS
      advCFD = (F2-F1)/2./DS - (FF(4)-3*FF(3)+3*FF(2)-FF(1))/6./DS
      RETURN
      END FUNCTION advCFD
!======================================================================
!     One sided-difference
!       F0 is the outermost gridbox, F3 is the innermost.
!
!        |---+---+---+--- or ---+---+---+---|
!       F0  F1  F2                 F2  F1  F0 
!----------------------------------------------------------------------
      REAL FUNCTION SD(F0,F1,F2,DS,op)
      REAL    :: F0,F1,F2,op
      REAL    :: DS
      SD = op*(-3*F0+4*F1-F2)/2./DS
      RETURN
      END FUNCTION SD
!======================================================================
