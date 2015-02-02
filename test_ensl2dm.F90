! --------------------------------------------------------------------------
! Description:
!
!   need: ncfile="input.nc", obsfiles=(/.../)
!
! History:
!
!   2014-10-28:
!
!       [shensi]: First creation.
! --------------------------------------------------------------------------
#include "inc_common.fpp"
program test_ensl2dm
use utility
use module_types
use module_ensemble
use module_modelstate
use module_obs
use module_config
implicit none
    type(type_modeldesc) :: md
    type(type_obsdesc) :: od
    type(type_localdesc) :: ld
    type(type_ensPx) :: ensPx
    type(type_ensPy) :: ensPy
    real,allocatable :: x(:), y(:), stdv(:)
    real :: y_max(5000)
    character(len=*),parameter :: ncfile="input.nc"
    character(len=*),parameter :: obsfiles(7)=(/"qg_obs.00000.txt", &
                        "qg_obs.00001.txt", "qg_obs.00002.txt", "qg_obs.00003.txt", &
                        "qg_obs.00004.txt", "qg_obs.00005.txt", "qg_obs.00006.txt"/)
    integer,parameter :: obs_tidx(7)=(/1,2,3,4,5,6,7/)
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
                (/"spol200", "smid200", "equat200", "nmid200", "npol200", &
                  "spol500", "smid500", "equat500", "nmid500", "npol500", &
                  "spol800", "smid800", "equat800", "nmid800", "npol800" /)
    integer :: i, j, idx
    real, allocatable :: row(:), column(:)
    character(len=80) :: corr_file


    call modeldesc_init(md)
    call modeldesc_getcoord(md, ncfile)
    !!read/write x test, passed
    !allocate(x(md%n_total))
    !call modelstate_input(md, ncfile, x)
    !call modelstate_output(md, "test_output.nc", x)
    !!Px
    call Px_read(md, (/ncfile/), ensPx)
    ensPx%Px=1.

    !!read y test, passed (basic check)
    call obs_read_txt(od, y_max, obsfiles, obs_tidx)
    allocate(y(od%n_obs))
    y=y_max(1:od%n_obs)
    !print*, od%varname((/1,180,181,360,361,370,1260/))
    !print*, od%tidx((/1,180,181,360,361,370,1260/))
    !print*, y((/1,180,181,360,361,370,1260/))
    !!Py
    call Py_read(od, reshape(obsfiles, (/7,1/)), obs_tidx, ensPy)
    print*, "after read Px/Py"
    !!localize
    call localdesc_setup(ld)
!!    ld%ishorizon=.false.
    ld%use_l2dm=.true.
    ld%l2dm_method="readinterp"
    ld%l2dm_dir="2dmodes"
!    ld%xmode_file=trim(xlmfile)
!    ld%ymode_file=trim(ylmfile)
    call PxPy_localize(ensPx, ensPy, md, od, ld)
    print*, "after read Localize Px/Py"
    print*, "n_state=",ensPx%n_state, "n_total_mem=",ensPx%n_total_mem
    allocate(stdv(ensPx%n_state))
    !! test sqrt(diag(Sx*Sx'))=sqrt(diag(Px*Px')) (here , n_mem=1) = abs(Px)=1 : passed
    call Px_get_stdv(md, ensPx, stdv) 
    print*, "after read get Px stdv"
    call system("cp "//trim(ncfile)//" ensl2dm_stdv.nc")
    call modelstate_output(md, "ensl2dm_stdv.nc", stdv)
!
!
    !!test multiply_vector: passed
  !  allocate(row(ensPx%n_total_mem))
  !  allocate(column(ensPx%n_state))
  !  row=0.1
  !  call Px_multiply_vector(md, ensPx, row, column)
  !  call system("cp erai_input.nc test_multiply.nc")
 !   call modelstate_output(md, "test_multiply.nc", column)
    !stop


!    !! test localization correlation, passed
    allocate(row(ensPx%n_total_mem))
    allocate(column(ensPx%n_state))
    do i=1, nvar
        do j=1, npts
            idx=modelstate_get_idx(md, i, ref_xid(j), ref_yid(j), ref_zid(j))
          !  print*, i, j, idx, ref_xid(j), ref_yid(j), stdv(idx)
            call Px_get_row(md, ensPx, i, ref_xid(j), ref_yid(j), ref_zid(j), row)
            print*,"after Px get row"
            call Px_multiply_vector(md, ensPx, row, column)
            column=column/((ensPx%n_mem))
            column=column
            print*,column(md%vars(i)%ibgn-1+1:md%vars(i)%ibgn-1+10)
            write(corr_file, *) "corr_", trim(var_names(i)), "_"//trim(point_name(j))//".nc"
            corr_file=trim(adjustl(corr_file))
            call system("cp  "//ncfile//" "//trim(corr_file))
            call modelstate_output(md, trim(corr_file), column)
        end do
    end do



end 
