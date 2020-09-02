!     -----------------------------------------------------------------
      SUBROUTINE OUTPUT(U,V,H,Z,cor,str,time,ana)
      INTEGER :: NX,NY,NR
      REAL,DIMENSION(NX,NY) :: U,V,H,Z,cor,str
      REAL,DIMENSION(NX,NY) :: spd,vor,WT,WR,PV
      REAL,DIMENSION(NR)    :: UA,VA,HA,spdA,vorA,strA,WTA,WRA
      INTEGER :: xSpd,ySpd,xVor,yVor,xH,yH
      REAL    :: mSpd,mVor,mH
      REAL    :: UAsym,VAsym
      LOGICAL :: ana
      CHARACTER(LEN=12) time
      CHARACTER(LEN=50) fn
      COMMON /grid/ NX,NY
      COMMON /gAM/  NR

 920  FORMAT(A12,3(2I9,f15.7),2f15.7)

      fn="SWM_output_"//time//".bin"
      WRITE(*,*) "OUTPUT: "//TRIM(fn)
      OPEN(21,FILE=fn,STATUS='REPLACE', &
            ACCESS="DIRECT", FORM="UNFORMATTED",RECL=NX*NY)
      WRITE(21,REC=1) U
      WRITE(21,REC=2) V
      WRITE(21,REC=3) H
      IF (ana) THEN

        spd = sqrt(U*U+V*V)

        WRITE(*,*) "  Calculating Vorticity..."
        CALL calcVorticity(U,V,vor)
        WRITE(*,*) "  Calculating StreamFunction..."
        CALL calcStream   (U,V,vor,str)
        WRITE(*,*) "  Calculating PV..."
        CALL calcPV       (vor,cor,H,Z,PV)
        WRITE(21,REC=4) spd
        WRITE(21,REC=5) vor
        WRITE(21,REC=6) str
        WRITE(21,REC=7) PV 
     
      ! only work for 1 TC
        WRITE(*,*) "  Tracking..."
        mH   = min2d_ind(  H,NX,NY,xH  ,yH  )
        mSpd = max2d_ind(spd,NX,NY,xSpd,ySpd)
        mVor = max2d_ind(vor,NX,NY,xVor,yVor)
        WRITE(*,*) "  Calculating Tagential/Radial Wind..."
        CALL calcWTC(U,V,spd,xVor,yVor,WT,WR)
        WRITE(21,REC=8) WT
        WRITE(21,REC=9) WR

        WRITE(*,*) "  Calculating AziMean..."
        CALL CalcAziMean( WT,xVor,yVor, WTA)
        CALL CalcAziMean( WR,xVor,yVor, WRA)
        CALL CalcAziMean(  H,xVor,yVor,  HA)
        CALL CalcAziMean(vor,xVor,yVor,vorA)
        CALL OUTPUT_AM("SWM_AM_WT.txt" , WTA)
        CALL OUTPUT_AM("SWM_AM_WR.txt" , WRA)
        CALL OUTPUT_AM("SWM_AM_H.txt"  ,  HA)
        CALL OUTPUT_AM("SWM_AM_VOR.txt",vorA)

        WRITE(*,*) "  Calculating AsymWind..."
        CALL calcAsymWind(WTA,U,V,xVor,yVor,UAsym,VAsym)
        OPEN (101,FILE="SWM_TC.txt",STATUS='OLD',POSITION="APPEND")
        WRITE(101,920) time,xH,yH,mH,xSpd,ySpd,mSpd,xVor,yVor,mVor, &
                       UAsym,VAsym 
        CLOSE(101)
      ENDIF
      CLOSE(21)
      END SUBROUTINE OUTPUT
!     =================================================================
      SUBROUTINE OUTPUT_Z(Z)
      INTEGER :: NX,NY
      REAL,DIMENSION(NX,NY) :: Z
      COMMON /grid/ NX,NY

      open(22,file="SWM_output_Z.bin",status='replace', &
             access="direct", form="unformatted",RECL=NX*NY) 
      write(22,rec=1) Z
      close(22)

      END SUBROUTINE OUTPUT_Z
!     =================================================================
      SUBROUTINE OUTPUT_Var(fn,Z)
      INTEGER :: NX,NY
      REAL,DIMENSION(NX,NY) :: Z
      CHARACTER(LEN=*) fn
      COMMON /grid/ NX,NY

      WRITE(*,*) "OUTPUT: "//trim(fn)
      open(23,file=fn,status='replace', &
             access="direct", form="unformatted",RECL=NX*NY) 
      write(23,rec=1) Z
      close(23)

      END SUBROUTINE OUTPUT_Var
!     =================================================================
      SUBROUTINE OUTPUT_AM(fn,Var)
      INTEGER :: NR
      REAL,DIMENSION(NR) :: Var
      CHARACTER(LEN=*)  :: fn
      CHARACTER(LEN=20) :: cfmt
      COMMON /gAM/  NR

      WRITE(cfmt,'(A,I5,A)') "(",NR,"f15.7)"

      OPEN (24,FILE=fn,STATUS='OLD',POSITION="APPEND")
      WRITE(24,cfmt) Var
      CLOSE(24)
      
      END SUBROUTINE OUTPUT_AM
!     =================================================================
      SUBROUTINE ReplaceFile(fn)
      CHARACTER(LEN=*) fn
      OPEN (25,FILE=fn,STATUS='REPLACE')
      CLOSE(25)
      END SUBROUTINE ReplaceFile
!     =================================================================

