C     $Id: init_ens_pdaf_read.F 18 2018-02-06 08:14:41Z lnerger $
#include "CPP_OPTIONS.h"
C
C     !ROUTINE: init_ens_pdaf_read --- Initialize ensemble
C
C     !INTERFACE:
      SUBROUTINE init_ens_pdaf_read(filtertype, dim_p, dim_ens, state_p,
     &     Uinv, ens_p, flag)

C     !DESCRIPTION:
C     User-supplied routine for PDAF.
C
C     The routine is called when the filter is
C     initialized in PDAF\_filter\_init.  It has
C     to initialize an ensemble of dim\_ens states.
C     Typically, the ensemble will be directly read from files.
C
C     The routine is called by all filter processes and 
C     initializes the ensemble for the PE-local domain.
C
C     !REVISION HISTORY:
C     2013-02 - Lars Nerger - Initial code
C     Later revisions - see svn log
C
C     !USES:
      USE mod_parallel_pdaf, 
     &     ONLY: mype_filter
      USE mod_assimilation,
     &     ONLY: myIter, myThid, read_ens

      IMPLICIT NONE

C     !ARGUMENTS:
      INTEGER, INTENT(in) :: filtertype                  ! Type of filter to initialize
      INTEGER, INTENT(in) :: dim_p                       ! PE-local state dimension
      INTEGER, INTENT(in) :: dim_ens                     ! Size of ensemble
      REAL*8, INTENT(inout) :: state_p(dim_p)            ! PE-local model state
C     It is not necessary to initialize the array 'state_p' for SEIK. 
C     It is available here only for convenience and can be used freely.
      REAL*8, INTENT(inout) :: Uinv(dim_ens-1,dim_ens-1) ! Array not referenced for SEIK
      REAL*8, INTENT(out)   :: ens_p(dim_p, dim_ens)     ! PE-local state ensemble
      INTEGER, INTENT(inout) :: flag                     ! PDAF status flag

C     !CALLING SEQUENCE:
C     Called by: PDAF_filter_init    (as U_ens_init)


C     *** local variables ***
#include "SIZE.h"
#include "EEPARAMS.h"
#include "DYNVARS.h"
#include "PARAMS.h"
#include "GRID.h"

      INTEGER :: i, j, k, bi, bj, member  ! Counters
      INTEGER :: localWarnings
      REAL*8 :: Tfreezing
      CHARACTER(len=4) :: ensstr          ! String for ensemble member
      CHARACTER*(MAX_LEN_MBUF) msgBuf

C     i,j,k  :: Loop counters
C     bi,bj  :: tile indices

C     **********************
C     *** INITIALIZATION ***
C     **********************

C     *** Generate full ensemble on filter-PE 0 ***
      IF (mype_filter==0) THEN
         WRITE (*, '(/a, 9x, a)') 'PDAFuser: ',
     &        'Initialize state ensemble'
         WRITE (*, '(a, 9x, a)') 'PDAFuser: ',
     &        '--- read ensemble from files'
         WRITE (*, '(a, 9x, a, i5)') 'PDAFuser: ',
     &        '--- Ensemble size:  ', dim_ens
      END IF

C     ********************************
C     *** Read ensemble from files ***
C     ********************************
      IF (read_ens) THEN
C     Set model variables to initial/restart conditions
       DO member = 1, dim_ens
        IF ( startTime .EQ. baseTime .AND.  nIter0 .EQ. 0
     &     .AND. pickupSuff .EQ. ' ' ) THEN
         WRITE(ensstr,'(A,I3.3)') '_',member
C     ---------
C     1. ETAN
C     ---------
C        Initialise surface position anomaly to zero
         DO bj = myByLo(myThid), myByHi(myThid)
          DO bi = myBxLo(myThid), myBxHi(myThid)
           DO j=1-OLy,sNy+OLy
            DO i=1-OLx,sNx+OLx
             etaN(i,j,bi,bj) = 0. _d 0
            ENDDO
           ENDDO
          ENDDO
         ENDDO
C        Read an initial state
         IF (pSurfInitFile .NE. ' ') THEN
          CALL READ_FLD_XY_RL( pSurfInitFile, ensstr, etaN, 0, myThid )
C         fill the overlap (+ BARRIER)
          _EXCH_XY_RL(etaN, myThid)
         ENDIF
C     ---------
C     2. THETA
C     ---------
C        Initialise temperature field to the vertical reference profile
         DO bj = myByLo(myThid), myByHi(myThid)
          DO bi = myBxLo(myThid), myBxHi(myThid)
           DO k=1,Nr
            DO j=1-OLy,sNy+OLy
             DO i=1-OLx,sNx+OLx
              theta(i,j,k,bi,bj) = tRef(k)
             ENDDO
            ENDDO
           ENDDO
          ENDDO
         ENDDO

         IF ( hydrogThetaFile .NE. ' ' ) THEN
          CALL READ_FLD_XYZ_RL( hydrogThetaFile, ensstr, theta,
     &                          0, myThid )
          _EXCH_XYZ_RL(theta,myThid)
         ENDIF

C--   Apply mask and test consistency
         localWarnings=0
         DO bj = myByLo(myThid), myByHi(myThid)
          DO bi = myBxLo(myThid), myBxHi(myThid)
           DO k=1,Nr
            IF ( maskIniTemp ) THEN
             DO j=1-OLy,sNy+OLy
              DO i=1-OLx,sNx+OLx
               IF (maskC(i,j,k,bi,bj).EQ.0.) theta(i,j,k,bi,bj) = 0.
              ENDDO
             ENDDO
            ENDIF
            IF ( tRef(k).NE.0. ) THEN
             DO j=1,sNy
              DO i=1,sNx
               IF (  maskC(i,j,k,bi,bj).NE.0.
     &         .AND. theta(i,j,k,bi,bj).EQ.0. ) THEN
                 localWarnings=localWarnings+1
               ENDIF
              ENDDO
             ENDDO
            ENDIF
           ENDDO
          ENDDO
         ENDDO
         IF (localWarnings.NE.0) THEN
          IF ( checkIniTemp ) THEN
           WRITE(msgBuf,'(A,I10,A)')
     &      ' INI_THETA: found', localWarnings,
     &      ' wet grid-pts with theta=0 identically.'
           CALL PRINT_ERROR( msgBuf , myThid)
           WRITE(msgBuf,'(A,A)')
     &     ' If this is intentional, you need to',
     &     ' set checkIniTemp=.FALSE. in "data", namelist PARM05'
           CALL PRINT_ERROR( msgBuf , myThid)
           STOP 'ABNORMAL END: S/R INI_THETA'
          ELSE
           WRITE(msgBuf,'(A,I10,A)')
     &      '** WARNINGS ** INI_THETA: found', localWarnings,
     &      ' wet grid-pts with theta=0 identically.'
           CALL PRINT_MESSAGE( msgBuf, errorMessageUnit,
     &                         SQUEEZE_RIGHT, myThid )
          ENDIF
         ENDIF
   
C--   Check that there are no values of temperature below freezing
Cpoint.
         IF ( checkIniTemp .AND. allowFreezing ) THEN
          Tfreezing=-1.9 _d 0
          DO bj = myByLo(myThid), myByHi(myThid)
           DO bi = myBxLo(myThid), myBxHi(myThid)
            DO k=1,Nr
             DO j=1-OLy,sNy+OLy
              DO i=1-OLx,sNx+OLx
               IF (theta(i,j,k,bi,bj) .LT. Tfreezing) THEN
                  theta(i,j,k,bi,bj) = Tfreezing
               ENDIF
              ENDDO
             ENDDO
            ENDDO
           ENDDO
          ENDDO
         ENDIF
C     ---------
C     3. SALT
C     ---------
C--   Initialise salinity field to the vertical reference profile
         DO bj = myByLo(myThid), myByHi(myThid)
          DO bi = myBxLo(myThid), myBxHi(myThid)
           DO k=1,Nr
            DO j=1-OLy,sNy+OLy
             DO i=1-OLx,sNx+OLx
              salt(i,j,k,bi,bj) = sRef(k)
             ENDDO
            ENDDO
           ENDDO
          ENDDO
         ENDDO
   
         IF ( hydrogSaltFile .NE. ' ' ) THEN
          CALL READ_FLD_XYZ_RL( hydrogSaltFile, ensstr, salt, 0, myThid)
          _EXCH_XYZ_RL( salt, myThid )
         ENDIF
   
C--   Apply mask and test consistency
         localWarnings=0
         DO bj = myByLo(myThid), myByHi(myThid)
          DO bi = myBxLo(myThid), myBxHi(myThid)
           DO k=1,Nr
            IF ( maskIniSalt ) THEN
             DO j=1-OLy,sNy+OLy
              DO i=1-OLx,sNx+OLx
               IF (maskC(i,j,k,bi,bj).EQ.0.) salt(i,j,k,bi,bj) = 0.
              ENDDO
             ENDDO
            ENDIF
            IF ( sRef(k).NE.0. ) THEN
             DO j=1,sNy
              DO i=1,sNx
               IF ( maskC(i,j,k,bi,bj).NE.0.
     &         .AND. salt(i,j,k,bi,bj).EQ.0. ) THEN
                 localWarnings=localWarnings+1
               ENDIF
              ENDDO
             ENDDO
            ENDIF
           ENDDO
          ENDDO
         ENDDO
         IF ( localWarnings.NE.0 ) THEN
          IF ( checkIniSalt ) THEN
           WRITE(msgBuf,'(A,I10,A)')
     &      ' INI_SALT: found', localWarnings,
     &      ' wet grid-pts with salt=0 identically.'
           CALL PRINT_ERROR( msgBuf , myThid)
           WRITE(msgBuf,'(A,A)')
     &     ' If this is intentional, you need to',
     &     ' set checkIniSalt=.false. in "data", namelist PARM05'
           CALL PRINT_ERROR( msgBuf , myThid)
           STOP 'ABNORMAL END: S/R INI_SALT'
          ELSE
           WRITE(msgBuf,'(A,I10,A)')
     &      '** WARNINGS ** INI_SALT: found', localWarnings,
     &      ' wet grid-pts with salt=0 identically.'
           CALL PRINT_MESSAGE( msgBuf, errorMessageUnit,
     &                         SQUEEZE_RIGHT, myThid )
          ENDIF
         ENDIF
C     ---------
C     4. VELOCITY
C     ---------
         DO bj = myByLo(myThid), myByHi(myThid)
          DO bi = myBxLo(myThid), myBxHi(myThid)
           DO k=1,Nr
            DO j=1-Oly,sNy+Oly
             DO i=1-Olx,sNx+Olx
              uVel(i,j,k,bi,bj)=0. _d 0
              vVel(i,j,k,bi,bj)=0. _d 0
             ENDDO
            ENDDO
           ENDDO
          ENDDO
         ENDDO
   
         IF (uVelInitFile .NE. ' ' .OR. vVelInitFile .NE. ' ') THEN
C     Read an initial state for each component if required
          IF (uVelInitFile .NE. ' ')
     &    CALL READ_FLD_XYZ_RL( uVelInitFile, ensstr, uVel, 0, myThid )
   
          IF (vVelInitFile .NE. ' ')
     &    CALL READ_FLD_XYZ_RL( vVelInitFile, ensstr, vVel, 0, myThid )
   
          CALL EXCH_UV_XYZ_RL(uVel,vVel,.TRUE.,myThid)
         ENDIF
   
         DO bj = myByLo(myThid), myByHi(myThid)
          DO bi = myBxLo(myThid), myBxHi(myThid)
           DO k=1,Nr
            DO j=1-Oly,sNy+Oly
             DO i=1-Olx,sNx+Olx
              uVel(i,j,k,bi,bj)=uVel(i,j,k,bi,bj)*_maskW(i,j,k,bi,bj)
              vVel(i,j,k,bi,bj)=vVel(i,j,k,bi,bj)*_maskS(i,j,k,bi,bj)
             ENDDO
            ENDDO
           ENDDO
          ENDDO
         ENDDO

         CALL collect_state_pdaf(dim_p, state_p)
         ens_p(:, member)  = state_p(:)

C      Or Read pickup files
        ELSEIF ( .NOT.useOffLine .OR. nonlinFreeSurf.GT.0 ) THEN
         CALL READ_PICKUP_ENS( nIter0, myThid, ensstr )
        ENDIF
       ENDDO
C     ********************************
C     *** Simply create ensemble   ***
C     ********************************
      ELSE

       CALL collect_state_pdaf(dim_p, state_p)
         
       DO member = 1, dim_ens
          ens_p(:, member)  = state_p(:) + float(member-1)*1.d-2
       ENDDO

      ENDIF

CC     Write ensemble
C      IF (write_ens) THEN
C
C         WRITE(suff,'(I10.10)') myIter
C
C         DO i=1,dim_ens
CC           map back to geometric/geographic fields
C            CALL distribute_state_pdaf(dim_p, ens_p(:,i))
C
C            WRITE(ensstr,'(I3.3)') i
C            CALL WRITE_FLD_XY_RL('ETAN_'//ensstr//fcstanls,
C     &           suff,etaN,myIter,myThid)
C            CALL WRITE_FLD_XYZ_RL('THETA_'//ensstr//fcstanls,
C     &           suff,theta,myIter,myThid)
C            CALL WRITE_FLD_XYZ_RL('SALT_'//ensstr//fcstanls,
C     &           suff,salt,myIter,myThid)
C            CALL WRITE_FLD_XYZ_RL('UVEL_'//ensstr//fcstanls,
C     &           suff,uvel,myIter,myThid)
C            CALL WRITE_FLD_XYZ_RL('VVEL_'//ensstr//fcstanls,
C     &           suff,vvel,myIter,myThid)
C
C
C         ENDDO
C      ENDIF




      
C     ****************
C     *** clean up ***
C     ****************

      END SUBROUTINE init_ens_pdaf_read
