! --------------------------------------------------------------------------
! Description:
!
!   calculate d(innovation/residual), output it with ob/ysim together into ob_ysim_d.csv
!   and calculate Jo(0)=0.5*d'*d, print to screen
!   ./calculate_d_Jo.exe <len_window, e.g. 6>  <prefix of ya/yb, so files are prefix??> <prefix of ob, i.e. ob??>
!
! History:
!
!   2014-09-28:
!
!       [shensi]: First creation.
! --------------------------------------------------------------------------
#include "inc_common.fpp"
program calculate_d_Jo
use utility
use ArgsManager
use module_types
use module_obs
implicit none
    character(len=*),parameter :: PROCEDURE_NAME="calculate_d_Jo"
    character(len=10) :: prefix_ysim, prefix_yo
    integer :: len_window, i, j
    character(len=20),allocatable :: ysim_files(:), yo_files(:)
    integer,parameter :: MAX_OBS=5000
    real :: y_all(MAX_OBS)
    real,allocatable :: yo(:), ysim(:), d(:)
    integer,allocatable :: tidx(:)
    type(type_obsdesc) :: od


    call args_readargs()
    if(args_get_narg() /= 3) then
        RP_ERROR("3 arguments need: <len_window> <prefix of y_sim> <prefix of yo>")
    end if
    call args_getarg(1, len_window)
    call args_getarg(2, prefix_ysim)
    call args_getarg(3, prefix_yo)

    allocate(ysim_files(len_window+1))
    allocate(yo_files(len_window+1))
    call get_files_byprefix(ysim_files, prefix_ysim)
    call get_files_byprefix(yo_files, prefix_yo)
    call check_file_exist_consistence(ysim_files, yo_files)

    allocate(tidx(len_window+1))
    tidx=(/(i-1,i=1,len_window+1)/)
    call obs_read_txt(od, y_all, ysim_files, tidx)
    allocate(ysim(od%n_obs), yo(od%n_obs), d(od%n_obs))
    ysim=y_all(1:od%n_obs)
    call obs_read_txt(od, yo, yo_files, tidx)
    d=yo-ysim

    call obs_normalize(od, ysim)
    call obs_normalize(od, yo)
    call obs_normalize(od, d)
    call write_float_matrix2( reshape((/yo, ysim, d/), (/od%n_obs,3/)), &
                        fname="ob_ysim_d.csv", element_fmat="F14.4,','")
    
    write(*,*) "Jo=", 0.5*sum(d*d)
    

contains
    !! files=(/ prefix//"00", prefix//"01" ... /)
    subroutine get_files_byprefix(files, prefix)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="get_files_byprefix"
        character(len=*) :: prefix
        character(len=*) :: files(:)
        integer :: n, i
        n=size(files)
        do i=1, n
            write(files(i), "(A,I2.2)") trim(prefix), i-1
        end do
    end subroutine get_files_byprefix  !!}}}

    subroutine check_file_exist_consistence(files1, files2)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="check_file_exist_consistence"
        character(len=*) :: files1(:), files2(:)
        integer :: n, i
        logical :: l1, l2
        n=size(files1)
        do i=1, n
            inquire(file=trim(files1(i)), exist=l1)
            inquire(file=trim(files2(i)), exist=l2)
            if(l1 /= l2 ) then
                RP_ERROR("file exist conflict!")
            end if
        end do
    end subroutine check_file_exist_consistence  !!}}}


end program calculate_d_Jo
