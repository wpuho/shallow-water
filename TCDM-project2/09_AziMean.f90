!======================================================================
! 20161029 -- Shen-Cha Hsu, schsu81@gmail.com
!======================================================================
      SUBROUTINE CalcAziMean(var,xc,yc,outvar)
      INTEGER NY,NX,NR
      REAL  outvar(NR)
      REAL  var(NX,NY)
      INTEGER xc,yc
      REAL  dr,RR,DX,DY
      REAL  sum_var,sum_wgt,wgt
      INTEGER t,r,g,mg,NG,tNG
      LOGICAL gflg
      COMMON /grid/ NX,NY
      COMMON /gAM/  NR
      COMMON /res/ DX,DY


      dr = max(DX,DY)
      RR = dr/2.

      mg   = 1
      DO r=1,NR
!      print*,"-----------------------"
        tNG     = 0
        sum_var = 0.
        sum_wgt = 0.
        gflg = .FALSE.
        DO g=mg,max(NX,NY)
          NG = 0
         ! Left and Right Lines
          DO i=xc-(g-1),xc+(g-1),max(2*(g-1),1)
          DO j=yc-(g-1),yc+(g-1)
            IF (i.LT.1 .OR. NX.LT.i .OR. j.LT.1 .OR. NY.LT.j) CYCLE
            wgt = get_wgt(i,j,xc,yc,r-1,dr,RR)
            IF (wgt.NE.0) NG = NG+1
            sum_var = sum_var + wgt*var(i,j)
            sum_wgt = sum_wgt + wgt
          ENDDO ! i
          ENDDO ! j
         ! Top and Bottom Lines
          DO j=yc-(g-1)  ,yc+(g-1)  ,max(2*(g-1),1)
          DO i=xc-(g-1)+1,xc+(g-1)-1
            IF (i.LT.1 .OR. NX.LT.i .OR. j.LT.1 .OR. NY.LT.j) CYCLE
            wgt = get_wgt(i,j,xc,yc,r-1,dr,RR)
            IF (wgt.NE.0) NG = NG+1
            sum_var = sum_var + wgt*var(i,j)
            sum_wgt = sum_wgt + wgt
          ENDDO ! i
          ENDDO !j

          IF     (.NOT.gFLG .AND. NG.NE.0)THEN
            mg = g
            gFLG = .TRUE.
          ELSEIF (gFLG .AND. NG.EQ.0) THEN
            EXIT!g
          ENDIF
          tNG = tNG + NG
        ENDDO ! g
        IF (sum_wgt.NE.0)THEN
        outvar(r) = sum_var/sum_wgt
        ENDIF
      ENDDO ! r

      RETURN
      END SUBROUTINE CalcAziMean
!======================================================================
      REAL FUNCTION get_wgt(i,j,xc,yc,R,dr,RR)
      INTEGER R,xc,yc,i,j 
      REAL    DX,DY
      COMMON /res/ DX,DY
      REAL dr,dis,RR

      dis = abs(R*dr - sqrt( ((i-xc)*DX)**2. + ((j-yc)*DY)**2.))
      IF (dis.ge.RR)THEN
      get_wgt = 0.
      ELSE
      get_wgt = 1 - dis/RR
      ENDIF
!      print*,i,j,xc,yc
!      print*,dis,RR
!      print*,"wgt=",get_wgt
      RETURN
      END FUNCTION get_wgt
!======================================================================
