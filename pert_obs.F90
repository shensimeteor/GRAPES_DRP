! --------------------------------------------------------------------------
! Description:
!
!   ./pert_obs.exe <input_obs_files> -o <output_obs_file_prefix>
!
! History:
!
!   2014-09-07:
!
!       [shensi]: First creation.
! --------------------------------------------------------------------------
#include "inc_common.fpp"
program pert_obs
use utility
use module_types
use module_obs
use ArgsManager
use mkl_vsl
use mkl_vsl_type
use vsl_external
implicit none
    character(len=*),parameter :: PROCEDURE_NAME="pert_obs"
    character(len=50) :: input, output
    integer :: narg, n_input, i, state
    integer,parameter :: MAX_OBS=5000
    real :: y(MAX_OBS), ptb(MAX_OBS)
!    character(len=2),parameter :: vars_name(3)=(/"gd","u","v"/)
!    real,parameter :: vars_obserror(3)=(/50, 2, 2/)  !!
    real :: hgt_obserr(3)=(/16, 9, 5/), psi_obserr(3)
    type(type_obsdesc) :: od
    type(vsl_stream_state) :: stream
    integer :: seed

    psi_obserr=hgt_obserr*10**5

    call args_readargs()
    narg=args_get_narg()
    if(narg /= 2) then
        RP_ERROR("./pert_obs.exe <input_obs_file> <output_obs_file> ")
    end if
    call args_getarg(1, input)
    call args_getarg(2, output)

    call obs_read_txt(od, y, (/trim(input)/), (/0/))
    
    seed=getseed_byclock()
    VSL_CHECK(vslnewstream(stream, VSL_BRNG_MCG31, seed))
    do i=1, od%n_obs
        VSL_CHECK(vsrnggaussian(VSL_METHOD_SGAUSSIAN_BOXMULLER2, stream, 1, ptb(i), 0., od%error(i)))
    end do
    y(1:od%n_obs)=y(1:od%n_obs)+ptb
    call obs_write_txt(od, y(1:od%n_obs), "temp_", 1)
    call system("mv temp_00 "//trim(output))

end program pert_obs
