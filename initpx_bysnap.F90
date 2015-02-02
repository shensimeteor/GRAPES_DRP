! --------------------------------------------------------------------------
! Description:
!
!   ./initpx_bysnap.exe <xbase_file> <N> <dir_to_place_output_ens> <logical, whether to minus the bias> <ratio>
!   ens = {xb+(xj-xi)*ratio-bias, i/j is randomly set in (b,1,2,..), }
!
! History:
!
!   2014-08-11:
!
!       [shensi]: First creation.
! --------------------------------------------------------------------------
#include "inc_common.fpp"
program initpx_bysnap
use module_modelstate
use utility
use ArgsManager
use mkl_vsl
use mkl_vsl_type
use vsl_external
implicit none
    character(len=*),parameter :: PROCEDURE_NAME="initpx_bysnap"
    character(len=100) :: xbfile, listfile, ens_file, dir
    character(len=100),allocatable :: snap_files(:)
    integer :: n_mem, i, n_snap, idx1, idx2
    type(type_modeldesc) :: md
    logical :: kickbias
    integer,parameter :: NS=64*32*3   !!n_state & n_grids
    real :: xb(NS), ptb(NS), xnew(NS), ftemp, ratio, pxmean(NS)
    real,allocatable :: snaps(:,:), px(:,:) 
    type(vsl_stream_state) :: stream
    integer :: seed

    !! read args
    call args_readargs()
    if(args_get_narg() /= 5) then
        print*, "./initpx_bysnap.exe <path_to_xb> <N> <dir_to_place_output_ens> <to_kickbias> <ratio>"
        stop 2
    end if
    listfile="list_files.txt"
    call args_getarg(1, xbfile)
    call args_getarg(2, n_mem)
    call args_getarg(3, dir)
    call args_getarg(4, kickbias)
    call args_getarg(5, ratio)

    !! prepare & read: xb * snaps
    call modeldesc_init(md)
    call modeldesc_getcoord(md, trim(xbfile))
    call modelstate_input(md, trim(xbfile), xb)
    call allocread_line(trim(listfile), n_snap, snap_files)
    n_snap=n_snap+1
    if(n_snap**2/2 < n_mem) then !! a enpirical method to jugde the n_mem's size
        RP_WARNING("n_mem is too large for the n_snap given")
    end if
    allocate(snaps(NS, n_snap))
    snaps(:,1)=xb
    do i=2, n_snap
        call modelstate_input(md, trim(snap_files(i-1)), snaps(:,i))
    end do

    !! generate ensemble
    allocate(px(NS,n_mem))
    call random_seed()    
    do i=1, n_mem
        call random_number(ftemp)
        ftemp=ftemp*(n_snap)
        idx1=ceiling(ftemp)
        call random_number(ftemp)
        ftemp=ftemp*(n_snap)
        idx2=ceiling(ftemp)
        px(:,i)=(snaps(:,idx1)-snaps(:,idx2))*ratio
    end do
    if(kickbias) then
        pxmean=sum(px, 2)/n_mem
    else
        pxmean=0
    end if
    do i=1, n_mem
        px(:,i)=px(:,i)-pxmean+xb
        write(ens_file,'(A,I3.3,A)') trim(dir)//"/xens_e",i,".nc"
        call system("cp "//trim(xbfile)//" "//trim(ens_file))
        call modelstate_output(md, trim(ens_file), px(:,i))
    end do
end program



    


