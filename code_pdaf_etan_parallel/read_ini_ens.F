#include "PACKAGES_CONFIG.h"
#include "CPP_OPTIONS.h"

C--  File read_ini_ens.F: Routines to handle 
C                         reading ensemble initial files
C--   Contents
C--   o INI_PSURF_ENS
C--   o INI_THETA_ENS
C--   o INI_SALT_ENS
C--   o INI_VEL_ENS

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C     !ROUTINE: INI_PSURF_ENS
C     !INTERFACE:
      SUBROUTINE INI_PSURF_ENS( myThid, ensstr)

C     !DESCRIPTION: 
C     *==========================================================*
C     | SUBROUTINE INI_PSURF_ENS                                 |
C     | o Set initial free-surface height/pressure for ensemble. |
C     *==========================================================*

C     !USES:
      IMPLICIT NONE
C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "SURFACE.h"
#ifdef ALLOW_CD_CODE
# include "CD_CODE_VARS.h"
#endif
#ifdef ALLOW_SHELFICE
# include "SHELFICE.h"
#endif

C     !INPUT/OUTPUT PARAMETERS:
C     == Routine arguments ==
C     myThid ::  Number of this instance of INI_PSURF
      INTEGER myThid
      CHARACTER(len=5) :: ensstr          ! String for ensemble member

C     !LOCAL VARIABLES:
C     == Local variables ==
C     bi,bj  :: tiles indices
C     I,J    :: Loop counters
      INTEGER bi, bj
      INTEGER  I,  J
CEOP

C     ==============================
C       ETAN and related variables
C     ==============================

C     Read an initial state
      CALL READ_FLD_XY_RL( pSurfInitFile, ensstr, etaN, 0, myThid )
C     fill the overlap (+ BARRIER)
      _EXCH_XY_RL(etaN, myThid)

#ifdef ALLOW_CD_CODE
C--   By default, initialize etaNm1 with etaN :
      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        DO J=1-OLy,sNy+OLy
         DO I=1-OLx,sNx+OLx
          etaNm1(I,J,bi,bj) = etaN(I,J,bi,bj)
         ENDDO
        ENDDO
       ENDDO
      ENDDO
C     _EXCH_XY_RL(etaNm1, myThid)
#endif

#ifdef EXACT_CONSERV
C--   By default, initialize etaH with etaN :
      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          etaH(i,j,bi,bj) = etaN(i,j,bi,bj)
          etaHnm1(i,j,bi,bj) = etaN(i,j,bi,bj)
          dEtaHdt(i,j,bi,bj) = 0. _d 0
         ENDDO
        ENDDO
       ENDDO
      ENDDO
#endif /* EXACT_CONSERV */

#ifdef ALLOW_SHELFICE
      IF ( useShelfIce .AND. usingZCoords ) THEN
        DO bj=myByLo(myThid),myByHi(myThid)
         DO bi=myBxLo(myThid),myBxHi(myThid)
           DO j=1-OLy,sNy+OLy
            DO i=1-OLx,sNx+OLx
              phi0surf(i,j,bi,bj) = phi0surf(i,j,bi,bj)
     &          + shelficeLoadAnomaly(i,j,bi,bj)*recip_rhoConst
            ENDDO
           ENDDO
         ENDDO
        ENDDO
      ENDIF
#endif /* ALLOW_SHELFICE */

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C     !ROUTINE: INI_THETA_ENS
C     !INTERFACE:
      SUBROUTINE INI_THETA_ENS( myThid, ensstr )

C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE INI_THETA
C     | o Set model initial temperature field.
C     *==========================================================*
C     | There are several options for setting the initial
C     | temperature file
C     |  1. Inline code
C     |  2. Vertical profile ( uniform T in X and Y )
C     |  3. Three-dimensional data from a file. For example from
C     |     Levitus or from a checkpoint file from a previous
C     |     integration.
C     | In addition to setting the temperature field we also
C     | set the initial temperature tendency term here.
C     *==========================================================*
C     \ev

C     !USES:
      IMPLICIT NONE
C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#ifdef ALLOW_MNC
#include "MNC_PARAMS.h"
#endif

C     !INPUT/OUTPUT PARAMETERS:
C     == Routine arguments ==
C     myThid :: Number of this instance of INI_THETA
      INTEGER myThid
      CHARACTER(len=5) :: ensstr          ! String for ensemble member

C     !LOCAL VARIABLES:
C     == Local variables ==
C     bi,bj  :: Tile indices
C     i,j,k  :: Loop counters
      INTEGER bi, bj
      INTEGER i, j, k, localWarnings
      _RL     Tfreezing
      CHARACTER*(MAX_LEN_MBUF) msgBuf
CEOP

      CALL READ_FLD_XYZ_RL( hydrogThetaFile, ensstr, theta, 0,
     &                      myThid )
      _EXCH_XYZ_RL(theta,myThid)

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
     &      .AND. theta(i,j,k,bi,bj).EQ.0. ) THEN
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
     &   ' INI_THETA: found', localWarnings,
     &   ' wet grid-pts with theta=0 identically.'
        CALL PRINT_ERROR( msgBuf , myThid)
        WRITE(msgBuf,'(A,A)')
     &  ' If this is intentional, you need to',
     &  ' set checkIniTemp=.FALSE. in "data", namelist PARM05'
        CALL PRINT_ERROR( msgBuf , myThid)
        STOP 'ABNORMAL END: S/R INI_THETA'
       ELSE
        WRITE(msgBuf,'(A,I10,A)')
     &   '** WARNINGS ** INI_THETA: found', localWarnings,
     &   ' wet grid-pts with theta=0 identically.'
        CALL PRINT_MESSAGE( msgBuf, errorMessageUnit,
     &                      SQUEEZE_RIGHT, myThid )
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
c     ELSEIF ( allowFreezing ) THEN
c      CALL FREEZE_SURFACE( startTime, nIter0, myThid )
      ENDIF

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C     !ROUTINE: INI_SALT_ENS
C     !INTERFACE:
      SUBROUTINE INI_SALT_ENS( myThid, ensstr )

C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE INI_SALT
C     | o Set model initial salinity field.
C     *==========================================================*
C     | There are several options for setting the initial
C     | temperature file
C     |  1. Inline code
C     |  2. Vertical profile ( uniform S in X and Y )
C     |  3. Three-dimensional data from a file. For example from
C     |     Levitus or from a checkpoint file from a previous
C     |     integration.
C     | In addition to setting the salinity field we also
C     | set the initial salinity tendency term here.
C     *==========================================================*
C     \ev

C     !USES:
      IMPLICIT NONE
C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#ifdef ALLOW_MNC
#include "MNC_PARAMS.h"
#endif

C     !INPUT/OUTPUT PARAMETERS:
C     == Routine arguments ==
C     myThid :: Number of this instance of INI_SALT
      INTEGER myThid
      CHARACTER(len=5) :: ensstr          ! String for ensemble member

C     !LOCAL VARIABLES:
C     == Local variables ==
C     bi,bj  :: Tile indices
C     i,j,k  :: Loop counters
      INTEGER bi, bj
      INTEGER i, j, k, localWarnings
      CHARACTER*(MAX_LEN_MBUF) msgBuf
CEOP

      CALL READ_FLD_XYZ_RL( hydrogSaltFile, ensstr, salt, 0,
     &                      myThid )
      _EXCH_XYZ_RL( salt, myThid )

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
     &      .AND. salt(i,j,k,bi,bj).EQ.0. ) THEN
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
     &   ' INI_SALT: found', localWarnings,
     &   ' wet grid-pts with salt=0 identically.'
        CALL PRINT_ERROR( msgBuf , myThid)
        WRITE(msgBuf,'(A,A)')
     &  ' If this is intentional, you need to',
     &  ' set checkIniSalt=.false. in "data", namelist PARM05'
        CALL PRINT_ERROR( msgBuf , myThid)
        STOP 'ABNORMAL END: S/R INI_SALT'
       ELSE
        WRITE(msgBuf,'(A,I10,A)')
     &   '** WARNINGS ** INI_SALT: found', localWarnings,
     &   ' wet grid-pts with salt=0 identically.'
        CALL PRINT_MESSAGE( msgBuf, errorMessageUnit,
     &                      SQUEEZE_RIGHT, myThid )
       ENDIF
      ENDIF

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C     !ROUTINE: INI_VEL_ENS
C     !INTERFACE:
      SUBROUTINE INI_VEL_ENS( myThid, ensstr )

C     !DESCRIPTION: \bv
C     *=================================================================
C     | SUBROUTINE INI_VEL
C     | o Initialize flow field (either to zero or from input files)
C     *=================================================================
C     \ev

C     !USES:
      IMPLICIT NONE
C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"

C     !INPUT/OUTPUT PARAMETERS:
C     == Routine arguments ==
C     myThid ::  Number of this instance of INI_VEL
      INTEGER myThid
      CHARACTER(len=5) :: ensstr          ! String for ensemble member

C     !LOCAL VARIABLES:
C     == Local variables ==
C     i,j,k  :: Loop counters
C     bi,bj  :: tile indices
      INTEGER bi,bj,i,j,k
CEOP
      DO bj = myByLo(myThid), myByHi(myThid)
       DO bi = myBxLo(myThid), myBxHi(myThid)
        DO k=1,Nr
         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx
           uVel(i,j,k,bi,bj) = 0. _d 0
           vVel(i,j,k,bi,bj) = 0. _d 0
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO

      IF (uVelInitFile .NE. ' ' .OR. vVelInitFile .NE. ' ') THEN
C     Read an initial state for each component if required
       IF (uVelInitFile .NE. ' ')
     & CALL READ_FLD_XYZ_RL( uVelInitFile, ensstr, uVel, 0, myThid )

       IF (vVelInitFile .NE. ' ')
     & CALL READ_FLD_XYZ_RL( vVelInitFile, ensstr, vVel, 0, myThid )

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

      RETURN
      END
