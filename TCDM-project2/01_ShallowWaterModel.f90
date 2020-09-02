      PROGRAM main
      INTEGER :: NX,NY,NR
      INTEGER :: NTC,NMNT
      INTEGER :: run_times,run_hours,run_days,out_times
      INTEGER :: itg_opt,smth_opt
      INTEGER :: ierr
      REAL    :: DX,DY,DT
      REAL,ALLOCATABLE,DIMENSION(:,:) :: U,V,H,Z,cor,str
      INTEGER :: FCST
      REAL    :: omg,pi,g
      LOGICAL :: ana
      CHARACTER(LEN=12) time
      COMMON /grid/ NX,NY
      COMMON /gAM/  NR
      COMMON /res/ DX,DY
      COMMON /time/ DT
      COMMON /cnst/ omg,pi,g
      NAMELIST /General/ NX,NY,DX,DY,DT,run_days,run_hours,out_times,&
                         itg_opt,smth_opt,NTC,NMNT,ana,NR 

      OPEN(11,FILE="namelist.SWM",STATUS="OLD",IOSTAT=ierr)
      IF (ierr.NE.0)THEN
        CALL createNML()
        STOP
      ELSE
        READ (11,NML=General)
        CLOSE(11)
        run_times = run_days*86400+run_hours*3600
      ENDIF
!     Parameter and initial values-------------------------------------
      omg  = 2*7.29E-5
      pi   = acos(-1.)
      g    = 9.8
      FCST = 0
      CALL fmtTime(FCST,time)
      ALLOCATE(U(NX,NY),V(NX,NY),H(NX,NY),Z(NX,NY),&
               cor(NX,NY),str(NX,NY))
      IF (ana) THEN
      CALL ReplaceFile("SWM_AM_WT.txt")
      CALL ReplaceFile("SWM_AM_WR.txt")
      CALL ReplaceFile("SWM_AM_H.txt")
      CALL ReplaceFile("SWM_AM_VOR.txt")
      OPEN(101,FILE="SWM_TC.txt"   ,STATUS="REPLACE") ! replace the file
      WRITE(101,'(A12,3(2A9,A15),2A15)') "Time","XminH","YminH","minH",&
                                       "XmaxWind","YminWind","maxWind",&
                                        "XmaxVor", "YmaxVor", "maxVor",&
                                        "UAsym","VAsym" ! write header
      CLOSE(101) 
      ENDIF
!     Initialization---------------------------------------------------
      CALL createENV(U,V,H,Z,cor)
      CALL createMNT(NMNT,Z) 
      CALL createTC (NTC,U,V,H) 
      WRITE(*,'(72("="))')
      CALL OUTPUT(U,V,H,Z,cor,str,time,ana) ! output initial condiction 
!     Integration -----------------------------------------------------
      DO WHILE (FCST .LT. run_times )
        FCST = FCST + DT 
        CALL fmtTime(FCST,time)
        WRITE(*,*) "Step=",FCST," Time=",time
!       Update- - - - - - - - - - - - - - - - - - -
        IF     (itg_opt.EQ.1)THEN
          CALL EL(U,V,H,Z,cor) ! Euler
        ELSEIF (itg_opt.EQ.2)THEN 
          CALL HN(U,V,H,Z,cor) ! Heun
        ELSEIF (itg_opt.EQ.3)THEN
          CALL RK(U,V,H,Z,cor) ! Runge-Kutta
        ELSE
          WRITE(*,*) "Uknown itg_opt"
          STOP
        ENDIF
!       Smooth- - - - - - - - - - - - - - - - - - -
        IF     (smth_opt.EQ.1)THEN
          CALL shapiro(U,(/0.125,-0.125/),2)
          CALL shapiro(V,(/0.125,-0.125/),2)
          CALL shapiro(H,(/0.125,-0.125/),2)
        ENDIF
!       - - - - - - - - - - - - - - - - - - - - - -
        IF (.NOT.chkStable(U,V,H)) THEN
          CALL OUTPUT(U,V,H,Z,cor,str,time,.False.)
          STOP
        ENDIF
!       Output- - - - - - - - - - - - - - - - - - -
        IF (mod(FCST,out_times).EQ.0) THEN
          CALL OUTPUT(U,V,H,Z,cor,str,time,ana) ! 1 hrs
        ENDIF
      ENDDO

      WRITE(*,*) "DONE"
      STOP
!     -----------------------------------------------------------------
      CONTAINS
        INCLUDE '02_init.f90'
        INCLUDE '03_output.f90'
        INCLUDE '04_analysis.f90'
        INCLUDE '05_difference.f90'
        INCLUDE '06_prognostic.f90'
        INCLUDE '07_integration.f90'
        INCLUDE '08_filter.f90'
        INCLUDE '09_AziMean.f90'
!       -----------------------------------------------------------------
!       Creating default namelist.
!       -----------------------------------------------------------------
        SUBROUTINE createNML()
        INTEGER :: NX,NY
        INTEGER :: NTC,NMNT
        INTEGER :: run_days,run_hours,out_times
        INTEGER :: itg_opt,smth_opt
        REAL    :: DX,DY,DT
        LOGICAL :: ana
        REAL    :: U0,V0,H0,lat0 ! env
        INTEGER :: Xc(2),Yc(2)
        REAL    :: ZM(2),Sx(2),Sy(2),ang(2) ! mnt
        REAL    :: VM(2),RM(2),ga(2)        ! tc
        NAMELIST /General/ NX,NY,DX,DY,DT,run_days,run_hours,out_times,&
                           itg_opt,smth_opt,NTC,NMNT,ana 
        NAMELIST /Enviroment/ U0,V0,H0,lat0
        NAMELIST /Mountains/ Xc,Yc,ZM,Sx,Sy,ang
        NAMELIST /Cyclone/ Xc,Yc,Vm,Rm,ga

        OPEN (10,FILE="namelist.SWM",STATUS="NEW",IOSTAT=ierr)

        NX =  401 ; NY =  401
        DX = 5000 ; DY = 5000
        DT = 5
        run_days  = 3
        run_hours = 0
        out_times = 3600
        itg_opt  = 1
        smth_opt = 1
        NTC  = 1
        NMNT = 1
        ana  = .True.
        WRITE(10,NML=General)

        U0 = -4.   ; V0 = 0.
        H0 = 5000. ; lat0 = 30
        WRITE(10,NML=Enviroment)

        Xc = (/151,201/)
        Yc = (/201,202/)
        ZM = (/3.   ,  2./)
        Sx = (/80.  , 85./)
        Sy = (/160. , 85./)
        ang= (/  0. , 15./)
        WRITE(10,NML=Mountains)

        Xc = (/261,202/)
        Yc = (/201,202/)
        VM = (/40,20/)
        RM = (/100. , 150./)
        ga = (/1.   ,  1.5/)
        WRITE(10,NML=Cyclone)
        STOP
        END SUBROUTINE createNML
!       -----------------------------------------------------------------
!       Converting simulation time to formated time.
!         Example : d00_00-58-10 for 1 day, 0 hour, 58 minite and 10 second.
!       -----------------------------------------------------------------
        SUBROUTINE fmtTime(FCST,time)
        INTEGER      :: FCST,dd,hh,mm,ss
        CHARACTER(LEN=12) time
 101    FORMAT("d",i2.2,"_",2(i2.2,"-"),i2.2)
        ss = mod(FCST     ,60)
        mm = mod(FCST/60  ,60)
        hh = mod(FCST/3600,24)
        dd =     FCST/86400
        WRITE(time,101) dd,hh,mm,ss
        RETURN 
        END SUBROUTINE fmtTime
!       -----------------------------------------------------------------
!       Whether the model is stable or not.
!         Accroding to CFL and nan.
!       -----------------------------------------------------------------
        LOGICAL FUNCTION chkStable(U,V,H)
        INTEGER :: NX,NY
        REAL    :: DX,DY,DT
        REAL    :: temp
        REAL,DIMENSION(NX,NY) :: U,V,H
        COMMON /grid/ NX,NY
        COMMON /res/ DX,DY
        COMMON /time/ DT

 103    FORMAT("X=",i3," Y=",i3," CFL=",f7.4," U=",f7.2," H=",f7.2)
        chkStable = .TRUE.
        DO i=1,NX
        DO j=1,NY
          temp = CFL(U(i,j),DX,DT)
          IF (temp .GT.1 .OR. isnan(H(i,j))) THEN
            chkStable = .FALSE.
            WRITE(*,103) i,j,temp,U(i,j),H(i,j)
          END IF
        ENDDO
        ENDDO

        RETURN
        END FUNCTION chkStable
!       -----------------------------------------------------------------
        REAL FUNCTION CFL(C,DX,DT)
        REAL    DX,DT
        REAL    C

        CFL = C*DT/DX

        RETURN
        END FUNCTION CFL
!       -----------------------------------------------------------------
!       Whether the value is infanity or not.
!         Note : invalid syntax in GNU fortran.
!       -----------------------------------------------------------------
!        LOGICAL FUNCTION isInf(input)
!        REAL :: inf
!        REAL :: input
!        inf=1./0.
!        IF (input.EQ.inf .OR. input.EQ.-inf)THEN
!         isInf = .TRUE. 
!        ELSE
!         isInf = .FALSE.
!        ENDIF
!        RETURN
!        END FUNCTION
!       -----------------------------------------------------------------
      END PROGRAM main
!     =================================================================

