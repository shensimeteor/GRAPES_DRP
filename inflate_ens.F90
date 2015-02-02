! --------------------------------------------------------------------------
! Description:
!
!   ./inflate_ens.exe <ens_list> <inflator> <out_dir> 
!
! History:
!
!   2014-09-26:
!
!       [shensi]: First creation.
! --------------------------------------------------------------------------
#include "inc_common.fpp"
program inflate_ens
use utility
use ArgsManager
use module_types
use module_modelstate
use module_ensemble
implicit none
    character(len=*),parameter :: PROCEDURE_NAME="inflate_ens"
    type(type_modeldesc) :: md
    type(type_ensPx) :: ensPx
    character(len=200) :: list_px, outdir
    real :: inflator
    character(len=200), allocatable :: pxfiles(:), newfiles(:)
    integer :: n, i, t, npx
    real,allocatable :: x_mean(:)

    call args_readargs()
    n=args_get_narg()
    if(n/=3) then
        RP_ERROR("3 arguments need: <list_ens> <inflator> <outdir>")
    end if
    call args_getarg(1, list_px)
    call args_getarg(2, inflator)
    call args_getarg(3, outdir)
    TRACEX(inflator)

    call allocread_line(trim(list_px), npx, pxfiles)

    call modeldesc_init(md)
    call modeldesc_getcoord(md, trim(pxfiles(1)))
    call Px_read(md, pxfiles, ensPx)
    allocate(x_mean(md%n_total))
    call Px_average(ensPx, x_mean)
    call Px_full2pert(ensPx, x_mean)
!! inflation
    do i=1, npx
        ensPx%Px(:,i)=ensPx%Px(:,i)*sqrt(inflator)
    end do
    
    call Px_pert2full(ensPx, x_mean)

    allocate(newfiles(npx))
    do i=1, npx
        newfiles(i)=trim(outdir)//'/'//trim(kickout_filepath(trim(pxfiles(i))))
        call system("cp "//trim(pxfiles(i))//" "//trim(newfiles(i)))
    end do
    call Px_write(md, newfiles, ensPx)

end program inflate_ens
