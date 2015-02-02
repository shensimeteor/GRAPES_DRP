! --------------------------------------------------------------------------
! Description:
!
!   kick_ens_bias.exe <xb> <list_ens.txt> <outdir>
!
! History:
!
!   2014-9-14
!
!       [shensi]: First creation.
! --------------------------------------------------------------------------
#include "inc_common.fpp"
program kick_ens_bias
use module_modelstate
use module_ensemble
use utility
use ArgsManager
implicit none
    character(len=*),parameter :: PROCEDURE_NAME="kick_ens_bias"
    character(len=100) :: xbfile, listfile, out_dir
    character(len=100), allocatable :: ens_files(:), out_files(:)
    integer :: n_mem, i
    type(type_modeldesc) :: md
    type(type_ensPx) :: ensPx
    integer,parameter :: NS=64*32*3 !!n_state & n_grids
    real :: xb(NS), pxmean(NS)

    !! read args
    call args_readargs()
    if(args_get_narg() /= 3) then
        print*, "./kick_ens_bias.exe <xbfile> <list_ens> <out_dir>"
        stop 2
    end if
    call args_getarg(1, xbfile)
    call args_getarg(2, listfile)
    call args_getarg(3, out_dir)

    !! prepare & read: xb * snaps
    call modeldesc_init(md)
    call modeldesc_getcoord(md, trim(xbfile))
    call modelstate_input(md, trim(xbfile), xb)
    call allocread_line(trim(listfile), n_mem, ens_files)
    
    call Px_read(md, ens_files, ensPx)
    call Px_average(ensPx, pxmean)
    call Px_full2pert(ensPx, pxmean)

    do i=1, n_mem
        ensPx%Px(:,i)=xb+ensPx%Px(:,i)
    end do

    allocate(out_files(n_mem))
    do i=1, n_mem
        out_files(i)=kickout_filepath(trim(ens_files(i)))
        out_files(i)=trim(out_dir)//"/kb_"//trim(out_files(i))
        call system("cp  "//trim(xbfile)//" "//trim(out_files(i)))
    end do
    call Px_write(md, out_files, ensPx)

end program



    


