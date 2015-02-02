#include "inc_common.fpp"
program test_enslocalman
use utility
use module_types
use module_ensemble
use module_modelstate
implicit none
    character(len=*), parameter :: PROCEDURE_NAME="test_enslocalman"
    type(type_modeldesc) :: md
    type(type_localdesc) :: ld
    type(type_ensPx) :: ensPx
    real,allocatable :: x(:), y(:), stdv(:)
    real :: y_max(5000)
    character(len=*),parameter :: ncfile="input.nc"
  !  character(len=*),parameter :: obsfiles(7)=(/"qg_obs.00000.txt", &
  !                      "qg_obs.00001.txt", "qg_obs.00002.txt", "qg_obs.00003.txt", &
  !                      "qg_obs.00004.txt", "qg_obs.00005.txt", "qg_obs.00006.txt"/)
  !  integer,parameter :: obs_tidx(7)=(/1,2,3,4,5,6,7/)
  !  character(len=*),parameter :: xlmfile="global_xmodes.csv", &
  !                                ylmfile="global_ymodes.csv"
    integer,parameter :: nvar=1
    character(len=3) :: var_names(nvar)=(/"psi"/)
    integer,parameter :: npts=15
    integer :: ref_xid(npts)=1, &
               ref_yid(npts)=(/1, 9, 16, 24, 32, &
                               1, 9, 16, 24, 32, &
                               1, 9, 16, 24, 32/), &
               ref_zid(npts)=(/1,1,1,1,1, &
                               2,2,2,2,2, &
                               3,3,3,3,3/)
    character(len=7) :: point_name(npts)= &
                (/"spol200", "smid200", "equa200", "nmid200", "npol200", &
                  "spol500", "smid500", "equa500", "nmid500", "npol500", &
                  "spol800", "smid800", "equa800", "nmid800", "npol800" /)
    integer :: i, j, idx
    real, allocatable :: row(:), column(:)
    character(len=80) :: corr_file
    real :: temp
    real,allocatable :: Pxrow(:), Pxrow2(:)
    integer :: ix=15, iy=10, iz=1, igrid
    type(type_localmodes) :: xlm, ylm
    real,allocatable :: Cor(:,:)


    call modeldesc_init(md)
    call modeldesc_getcoord(md, ncfile)
    call Px_read(md, (/ncfile/), ensPx)
    ensPx%Px=1.
    ld%use_localman=.true.
    ld%ishorizon=.true.
    ld%isvertical = .true.
    ld%lman_xm_file="global_xm.nc"
    ld%lman_ym_file="global_ym.nc"
    ld%lman_zm_file="global_zm.nc"

    call Px_localize(ensPx, ld)
    print*, "after read Localize Px"
    print*, "n_state=",ensPx%n_state, "n_total_mem=",ensPx%n_total_mem
 !   !! test ROUx * ROUx', passed
 !   xlm=ensPx%lman%x3dm(iy,iz)
 !   allocate(Cor(ensPx%lman%nx, ensPx%lman%nx))
 !   Cor=matmul(xlm%modes, transpose(xlm%modes))
 !   print*,Cor(10:15,10:15)
 !   !! test ROUy*ROUy'
 !   ylm=ensPx%lman%y2dm(iz)
 !   deallocate(Cor)
 !   allocate(Cor(ensPx%lman%ny, ensPx%lman%ny))
 !   Cor=matmul(ylm%modes, transpose(ylm%modes))
 !   print*, Cor(10:15,10:15)
 !   stop

    !! test Px_get_row
!    allocate(row(ensPx%n_total_mem), Pxrow(ensPx%n_state), column(ensPx%n_state))
!    !! 1. Pxrow=ensPx^T * (0,0,...1,...0)^T
!    column=0
!    igrid=modelstate_get_idx(md, 1, ix, iy, iz)
!    column(igrid)=1
!    call Px_multiply_vector(md, ensPx, column, Pxrow)
!    !! 2. Pxrow2=Px_get_row
!    call Px_get_row(md, ensPx, 1, ix, iy, iz, Pxrow2)
 !   print*, sum(abs(Pxrow-Pxrow2))/size(Pxrow)
!    stop

!    do i=1, ensPx%lman%ny
!        temp=sum(ensPx%lman%x3dm(i,1)%modes(1,:)**2)
!!        print*,i,temp
!    end do
!    do i=1, ensPx%lman%ny
!        temp=sum(ensPx%lman%y2dm(1)%modes(i,:)**2)
!        print*, i, temp
!    end do

    allocate(stdv(ensPx%n_state))
!    !! test sqrt(diag(Sx*Sx'))=sqrt(diag(Px*Px')) (here , n_mem=1) = abs(Px) : passed
!    call Px_get_stdv_px(ensPx, stdv) 
!    TRACEX(ensPx%n_mem)
!    TRACEX(get_average(ensPx%Px(:,1)))
!    TRACEX(get_average(stdv))
!    print*, "after read get Px stdv"
!    call system("cp  "//trim(ncfile)//" "//trim("enslocalman_pxstdv.nc"))
!    call modelstate_output(md, "enslocalman_pxstdv.nc", stdv)
    call Px_get_stdv(md, ensPx, stdv)
    call system("cp  "//trim(ncfile)//" "//trim("enslocalman_stdv.nc"))
    call modelstate_output(md, "enslocalman_stdv.nc", stdv)


    !!test multiply_vector: passed for nonlocal, unknown for local 
!   allocate(row(ensPx%n_total_mem))
!   allocate(column(ensPx%n_state))
!   row=0.1
!   call Px_multiply_vector(md, ensPx, row, column)
!   call system("cp  "//trim(ncfile)//" "//trim("test_multiply.nc"))
!   call modelstate_output(md, "test_multiply.nc", column)
!   call Px_get_row(md, ensPx, 1, 32, 13, 2, row)
!   print*, row


    !! test localization correlation, passed
    allocate(row(ensPx%n_total_mem))
    allocate(column(ensPx%n_state))
    do i=1, nvar
        do j=1, npts
            idx=modelstate_get_idx(md, i, ref_xid(j), ref_yid(j), ref_zid(j)) 
            print*," -------POINTS ",j,"----------------"
            write(*,'(6I,1F)') i, j, idx, ref_xid(j), ref_yid(j), ref_zid(j), stdv(idx)
            call Px_get_row(md, ensPx, i, ref_xid(j), ref_yid(j), ref_zid(j), row)
!            open(11, file="row_pt"//trim(tostring(ival=j))//".txt")
!            write(11, '(F)') row
!            close(11)
            print*,"after Px get row"
            call Px_multiply_vector(md, ensPx, row, column)
            column=column/((ensPx%n_mem))
            column=column/(stdv(idx)*stdv) !!
   !         column=column
    !        print*,column(md%vars(i)%ibgn-1+1:md%vars(i)%ibgn-1+10)
            write(corr_file, *) "corr_", trim(var_names(i)), "_"//trim(point_name(j))//".nc"
            corr_file=trim(adjustl(corr_file))
            call system("cp  "//ncfile//" "//trim(corr_file))
            call modelstate_output(md, trim(corr_file), column)
        end do
    end do

end 
