!-------------------------------------------------------------------------------------------------------------------------------
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!                       Copyright (C) 2002 by
!   National Institute of Public Health and Environment
!           Laboratory for Air Research (RIVM/LLO)
!                      The Netherlands
!
! SUBROUTINE
! NAME               : %M%
! SCCS (SOURCE)      : %P%
! RELEASE - LEVEL    : %R% - %L%
! BRANCH - SEQUENCE  : %B% - %S%
! DATE - TIME        : %E% - %U%
! WHAT               : %W%:%E%
! AUTHOR             : HvJ/Franka Loeve (Cap Volmac)
! FIRM/INSTITUTE     : RIVM/LLO
! LANGUAGE           : FORTRAN-77/90
! DESCRIPTION        : Compute scaling factors for printing of concentrations and depositions. A scaling factor is the ratio between the 
!                      computed concentration (or deposition) and the value to be printed. All ratio's are based on an input source strength
!                      in g/s and an output in ug/m3 for concentrations and mol/ha/y for depositions.
! EXIT CODES         :
! FILES AND OTHER    :
!    I/O DEVICES
! SYSTEM DEPENDENCIES: HP FORTRAN
! CALLED FUNCTIONS   :
! UPDATE HISTORY     :
!-------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ops_scalefac(nrrcp, cpri, csec, drydep, wetdep, scale_con, scale_sec, scale_dep, cseccor, scale_sec_cor)

USE m_commonconst                                                              ! EPS_DELTA only

IMPLICIT NONE

! CONSTANTS
CHARACTER*512                                    :: ROUTINENAAM                ! 
PARAMETER        (ROUTINENAAM = 'ops_scalefac')

! SUBROUTINE ARGUMENTS - INPUT
INTEGER*4, INTENT(IN)                            :: nrrcp                      ! number of receptor points
real,    INTENT(IN)                            :: cpri(nrrcp)                ! array van primaire concentraties
real,    INTENT(IN)                            :: csec(nrrcp)                ! array van secundaire concentraties
real,    INTENT(IN)                            :: drydep(nrrcp)              ! array van droge depositie
real,    INTENT(IN)                            :: wetdep(nrrcp)              ! array van natte depositie
real,    INTENT(IN), OPTIONAL                  :: cseccor(nrrcp)             ! concentration of second secondary substance

! SUBROUTINE ARGUMENTS - OUTPUT
real,    INTENT(OUT)                           :: scale_con                  ! schaal vergr. concentratie
real,    INTENT(OUT)                           :: scale_sec                  ! schaal vergr. secundaire concentratie
real,    INTENT(OUT)                           :: scale_dep                  ! schaal vergr. droge depositie
real,    INTENT(OUT), OPTIONAL                 :: scale_sec_cor              ! schaal vergr. secundaire concentratie in gasvorm

! LOCAL VARIABLES
INTEGER*4                                        :: i                          ! teller over schaalfactoren
real                                           :: cmax                       ! grootst voorkomende primaire concentratie
real                                           :: csmax                      ! grootst voorkomende secundaire concentratie
real                                           :: csgmax                     ! grootst voorkomende secundaire concentratie in gasvorm
real                                           :: ddepmax                    ! grootst voorkomende droge depositie
real                                           :: depntmax                   ! grootst voorkomende natte depositie
real                                           :: s                          ! schaalfactor
real                                           :: tc                         ! teller aantal te grote prim. conc.
real                                           :: td                         ! teller aantal te grote droge dep.
real                                           :: tn                         ! teller aantal te grote natte dep.
real                                           :: ts                         ! teller aantal te grote sec. conc.
real                                           :: tg                         ! teller aantal te grote sec. conc. gas
real                                           :: scale_dry                  ! schaal vergr. concentratie
real                                           :: scale_wet                  ! schaal vergr. concentratie

! SCCS-ID VARIABLES
CHARACTER*81                                     :: sccsida                    ! 
sccsida = '%W%:%E%'//char(0)
!-------------------------------------------------------------------------------------------------------------------------------
!
! Compute maxima of primary and secondary concentrations and dry and wet depositions
!
cmax     = MAXVAL(cpri(:))
csmax    = MAXVAL(csec(:))
IF (PRESENT(cseccor)) csgmax = MAXVAL(cseccor(:))
ddepmax  = MAXVAL(drydep(:))
depntmax = MAXVAL(wetdep(:))

!
! Initialise scaling factors for primary and secondary concentrations and dry and wet depositions
!
scale_con = 1.0e-10
scale_sec = 1.0e-10
IF (PRESENT(cseccor)) scale_sec_cor = 1.0e-10
scale_dry = 1.0e-10
scale_wet = 1.0e-10
! 
! Loop over exponent in scaling factor (scaling factors range from 1e-10 to 1e30);
! Set scaling factor s, such that for a parameter x with maximum xmax: s*xmax < 2000 (or (2000/xmax) > s)
!
DO i = -10, 30
  s = 10**(FLOAT(i))
  IF (cmax .GT. (0. + EPS_DELTA) .AND. (2000./cmax) .GT. (s + EPS_DELTA)) THEN
     scale_con = s
  ENDIF
  IF (csmax .GT. (0. + EPS_DELTA) .AND. (2000./csmax) .GT. (s + EPS_DELTA)) THEN
     scale_sec = s
  ENDIF
  IF (PRESENT(cseccor)) THEN
    IF (csgmax .GT. (0. + EPS_DELTA) .AND. (2000./csgmax) .GT. (s + EPS_DELTA)) THEN
       scale_sec_cor = s
    ENDIF
  ENDIF
  IF (ddepmax .GT. (0. + EPS_DELTA) .AND. (2000./ddepmax) .GT. (s + EPS_DELTA)) THEN
     scale_dry = s
  ENDIF
  IF (depntmax .GT. (0. + EPS_DELTA) .AND. (2000./depntmax) .GT. (s + EPS_DELTA)) THEN
     scale_wet = s
  ENDIF
ENDDO
!
! Count number of values > 999
!
tc = COUNT(cpri(:)*scale_con .GT. 999.+EPS_DELTA)
ts = COUNT(csec(:)*scale_sec .GT. 999.+EPS_DELTA)
IF (PRESENT(cseccor)) tg = COUNT(cseccor(:)*scale_sec_cor .GT. 999.+EPS_DELTA)
td = COUNT(drydep(:)*scale_dry .GT. 999.+EPS_DELTA)
tn = COUNT(wetdep(:)*scale_wet .GT. 999.+EPS_DELTA)
!
! If more than 5% of values are > 999, divide scaling factor by 10
!
IF (tc .GT. (0.05*nrrcp + EPS_DELTA)) THEN
  scale_con = scale_con/10.
ENDIF
IF (ts .GT. (0.05*nrrcp + EPS_DELTA)) THEN
  scale_sec = scale_sec/10.
ENDIF
IF (PRESENT(cseccor)) THEN
  IF (tg .GT. (0.05*nrrcp + EPS_DELTA)) THEN
    scale_sec_cor = scale_sec_cor/10.
  ENDIF
ENDIF
IF (td .GT. (0.05*nrrcp + EPS_DELTA)) THEN
  scale_dry = scale_dry/10.
ENDIF
IF (tn .GT. (0.05*nrrcp + EPS_DELTA)) THEN
  scale_wet = scale_wet/10.
ENDIF
scale_dep = AMIN1(scale_dry, scale_wet)

RETURN
END SUBROUTINE ops_scalefac
