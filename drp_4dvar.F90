! --------------------------------------------------------------------------
! Description:
!
!   the main program of drp-4dvar
!   input: 
!        xb:  xb
!        xg:  xg
!        <xb00>: xb00
!        <xg00>: xg00
!        px: full px; ENS/px.e001 ...
!        <px00>:  ENS/px00.e001 ...
!        py: full py; ENS/py00.e001 py02.e001 ...
!        ob: obs. ob00 ob02 ...
!        <yb>: yb00 yb02 ... (if not provided, use mean(py) when calcualte ptb_py & d)
!        <beta_g>: betag
!   output:
!        xa:  xa
!        <xa00>: xa00
!        <gx>: ENS/gx.e001 ...
!        <beta_g>: betag
!
! History:
!
!   2014-08-21:
!
!       [shensi]: First creation.
! --------------------------------------------------------------------------

#include "inc_common.fpp"
program drp_4dvar
use utility
use module_types
use module_config
use module_drpsolve
implicit none

    type(type_localdesc) ::ld
    type(type_drpconfig) :: dc
    TIME_START(timer_bgn)
    call config_get_namelist(dc, ld)
    call config_check_settings(dc, ld)
    TIME_CLICK(timer, timer_bgn, timer_last)
    
    call drp_run(dc, ld)

 end
    
