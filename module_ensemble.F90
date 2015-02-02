! --------------------------------------------------------------------------
! Description:
!
!   about Px/Py, including extended ensemble localization,
!   Px_read/Py_read/Px_write/Py_write : call md_/obs_ read & write, so model/obs dependent
!   Px_localize
!  ! should be model-independent
!   every about model/grid should be based on modeldesc/ obs based on obsdesc
!   the tolerent modification when model changes is
!   1. the interface of modelstate/obs may be change
!   2. the detailed modes to select may be change
!
! History:
!
!   2014-07-09:
!
!       [shensi]: First creation.
! --------------------------------------------------------------------------
#include "inc_common.fpp"
module module_ensemble
use utility
use interp
use module_types
use module_modelstate
use module_obs
use module_local
use module_localman
implicit none
private
    ! variables
    public :: Px_read, Py_read, Px_write, Px_localize, Py_localize, &
              Px_multiply_vector, Py_multiply_vector, Py_transpose_multiply_vector, &
              Px_get_row, Px_get_stdv, Py_get_row, Px_full2pert, Py_full2pert, &
              Px_average,  Px_delete, Py_delete, Px_get_stdv_px, Px_pert2full
    integer,private :: module_ens_cnt=0
contains

    ! methods
!!read Px/Py
    subroutine Px_read(md, files, ensPx)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="Px_read"
        type(type_modeldesc) :: md
        character(len=*) :: files(:)
        type(type_ensPx) :: ensPx
        integer :: n_state, n_mem, i
        n_mem=size(files)
        ensPx%n_mem=n_mem
        ensPx%n_state=md%n_total
        ensPx%n_total_mem=n_mem
        allocate(ensPx%Px(ensPx%n_state, ensPx%n_mem))
        do i=1, n_mem 
            call modelstate_input(md, trim(files(i)), ensPx%Px(:,i))
        end do
    end subroutine Px_read  !!}}}

    subroutine Px_write(md, files, ensPx)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="Px_write"
        type(type_modeldesc) :: md
        character(len=*) :: files(:)
        type(type_ensPx) :: ensPx
        integer :: i
        if(size(files) /= ensPx%n_mem) then
            RP_ERROR("size(files) /= ensPx%n_mem")
        end if
        do i=1, ensPx%n_mem 
            call modelstate_output(md, trim(files(i)), ensPx%Px(:,i))
        end do
    end subroutine Px_write  !!}}}

    !!filess is a two dim. array, (n_time, n_mem)
    subroutine Py_read(od, filess, tindices, ensPy)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="Py_read"
        type(type_obsdesc) :: od
        character(len=*) :: filess(:,:)
        type(type_ensPy) :: ensPy
        integer :: n_mem, n_time, tindices(:), i
        n_time=size(filess,1)
        n_mem=size(filess,2)
        ensPy%n_total_mem=n_mem
        ensPy%n_mem=n_mem
        ensPy%n_state=od%n_obs
        allocate(ensPy%Py(ensPy%n_state, ensPy%n_mem))
        do i=1, n_mem
            call obs_read_txt(od, ensPy%Py(:,i), filess(:,i), tindices)
        end do
    end subroutine Py_read  !!}}}


!!do localize( include none-localization)
    subroutine Px_localize(ensPx, ld)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="Px_localize"
        type(type_ensPx) :: ensPx
        type(type_localdesc) :: ld
        ensPx%n_total_mem=ensPx%n_mem
        if(ld%use_localman) then
            if(ld%ishorizon) then
                call localman_input_xm(ensPx%lman, trim(ld%lman_xm_file))
                call localman_input_ym(ensPx%lman, trim(ld%lman_ym_file))
               ! call localman_set_ym_nolocal(ensPx%lman)
                ensPx%n_total_mem = ensPx%n_total_mem * ensPx%lman%x3dm(1,1)%m * ensPx%lman%y2dm(1)%m
            else
                call localman_input_domain(ensPx%lman, trim(ld%lman_xm_file))
                call localman_set_xm_nolocal(ensPx%lman)
                call localman_set_ym_nolocal(ensPx%lman)
            end if
            if(ld%isvertical) then
                call localman_input_zm(ensPx%lman, trim(ld%lman_zm_file))
                ensPx%n_total_mem = ensPx%n_total_mem * ensPx%lman%zm%m
            else 
                call localman_set_zm_nolocal(ensPx%lman)
            end if
            call localman_rownormalize(ensPx%lman)
        else
            RP_ERROR("ld%use_localman should be True")
        end if
        ensPx%has_localize=.true.
    end subroutine Px_localize  !!}}}
    
    !! after ensPx localize
    subroutine Py_localize(ensPx, ensPy, od, ld)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="Py_localize"
        type(type_ensPx) :: ensPx
        type(type_ensPy) :: ensPy
        type(type_localdesc) :: ld
        type(type_obsdesc) :: od
        ASSUREX(ensPx%has_localize)
        if(ld%use_localman) then
            if(ld%ishorizon .or. ld%isvertical) then
                call localman_interp_toscatter(ensPx%lman, od%xpos, od%ypos, od%zpos, &
                                               ensPy%xlm, ensPy%ylm, ensPy%zlm)
                call localmodes_rownormalize(ensPy%xlm)
                call localmodes_rownormalize(ensPy%ylm)
                call localmodes_rownormalize(ensPy%zlm)
            else
                call localmodes_generate_nolocal(ensPy%xlm, od%n_obs)
                call localmodes_generate_nolocal(ensPy%ylm, od%n_obs)
                call localmodes_generate_nolocal(ensPy%zlm, od%n_obs)
            end if
        else
            RP_ERROR("ld%use_localman should be True")
        end if
        ensPy%n_total_mem = ensPx%n_total_mem
        ensPy%has_localize=.true.
    end subroutine Py_localize  !!}}}

!!after localization
    subroutine Py_multiply_vector(od, ensPy, vct, out_vct)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="Py_multiply_vector"
        type(type_obsdesc) :: od
        type(type_ensPy) :: ensPy
        real :: vct(:), out_vct(:)
        integer :: nvct, nout
        real,allocatable :: vct_mtrx(:,:,:,:)
        integer :: n1, n2, n3, n4, i, j, k, t, iobs
        nvct=size(vct)
        nout=size(out_vct)
        ASSUREX(ensPy%has_localize)
        ASSUREX(nvct == ensPy%n_total_mem)
        ASSUREX(nout == ensPy%n_state)
        n1=ensPy%n_mem
        n2=ensPy%xlm%m
        n3=ensPy%ylm%m
        n4=ensPy%zlm%m
        allocate(vct_mtrx(n4,n3,n2,n1))
        vct_mtrx=reshape(vct, (/n4,n3,n2,n1/))
        out_vct=0
        do iobs=1, od%n_obs
            do t=1, n1
                do i=1, n2
                    do j=1, n3
                        do k=1, n4
                            out_vct(iobs)=out_vct(iobs)+ensPy%Py(iobs,t) &
                                             *ensPy%xlm%modes(iobs,i) &
                                             *ensPy%ylm%modes(iobs,j) &
                                             *ensPy%zlm%modes(iobs,k) &
                                             *vct_mtrx(k,j,i,t)
                        end do
                    end do
                end do
            end do
        end do
        deallocate(vct_mtrx)
    end subroutine Py_multiply_vector  !!}}}
   
    subroutine Py_transpose_multiply_vector(od, ensPy, vct, out_vct)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="Py_transpose_multiply_vector"
        type(type_obsdesc) :: od
        type(type_ensPy) :: ensPy
        real :: vct(:), out_vct(:)
        integer :: nvct, nout
        real,allocatable :: vct_mtrx(:,:,:,:)
        integer :: n1, n2, n3, n4, i, j, k, t, iobs
        nvct=size(vct)
        nout=size(out_vct)
        ASSUREX(ensPy%has_localize)
        ASSUREX(ensPy%n_state == nvct)
        ASSUREX(ensPy%n_total_mem == nout)
        n1=ensPy%n_mem
        n2=ensPy%xlm%m
        n3=ensPy%ylm%m
        n4=ensPy%zlm%m
        allocate(vct_mtrx(n4,n3,n2,n1))
        vct_mtrx=0
        do t=1, n1
            do i=1, n2
                do j=1, n3
                    do k=1, n4
                        do iobs=1, od%n_obs
                            vct_mtrx(k,j,i,t)=vct_mtrx(k,j,i,t) + &
                                              ensPy%Py(iobs,t) * &
                                              ensPy%xlm%modes(iobs,i) * &
                                              ensPy%ylm%modes(iobs,j) * &
                                              ensPy%zlm%modes(iobs,k) * &
                                              vct(iobs)
                        end do
                    end do
                end do
            end do
        end do
        out_vct=reshape(vct_mtrx, (/n1*n2*n3*n4/))
        deallocate(vct_mtrx)
    end subroutine Py_transpose_multiply_vector  !!}}}

    subroutine Px_multiply_vector(md, ensPx, vct, out_vct)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="Px_multiply_vector"
        type(type_modeldesc) :: md
        type(type_ensPx) :: ensPx
        real :: vct(:), out_vct(:)
        integer :: nvct, nout, px, py, pz, n_ens, i, j, k, s
        integer :: ix, iy, iz, nx_hr, ny_hr, nz_hr, ivar
        real,allocatable :: f3d(:,:,:) !! (nx, my, mz): ROUx * vector
        real,allocatable :: f3d2(:,:,:) !! (nx, ny, mz): diag(ROUy_ipy) * vector
        real,allocatable :: f3d3(:,:,:) !! (nx, ny, nz): diag(ROUz_ipz) * vector
        real,allocatable :: f3d_sum(:,:,:) !! same with f3d3
        real,allocatable :: f3d_sum2(:,:,:) !! same with f3d2, test
        real,allocatable :: f3d_sum3(:,:,:) !! same with f3d2, test
        real,allocatable :: f3d_hr(:,:,:) !! High Resolution (i.e. on model grid)
        real,allocatable :: f3d_Px(:,:,:) !!  same with f3d_hr
        real,allocatable :: vct_mtrx(:,:,:,:) !! matrix form of vector (pz, py, px, nens)
        real,allocatable :: xweight(:,:,:), yweight(:,:,:), zweight(:,:,:) !! local-grid -> model-grid 
        integer,allocatable :: xneighbor(:,:,:), yneighbor(:,:,:), zneighbor(:,:,:) !! ...
        real,allocatable :: yweight_m2n(:,:), zweight_m2n(:,:) !! local-grid: my->ny
        integer,allocatable :: yneighbor_m2n(:,:), zneighbor_m2n(:,:) !! local-grid: mz->nz
        integer :: nx_big, ny_big, nz_big
        print*,"test: just get into Px_multiply_vector, before get neighbors & weights"
        nvct=size(vct)
        nout=size(out_vct)
        ASSUREX(ensPx%has_localize)
        ASSUREX(ensPx%n_total_mem == nvct)
        ASSUREX(ensPx%n_state == nout)
        allocate(f3d(ensPx%lman%nx, ensPx%lman%my, ensPx%lman%mz))
        allocate(f3d2(ensPx%lman%nx, ensPx%lman%ny, ensPx%lman%mz))
        allocate(f3d3(ensPx%lman%nx, ensPx%lman%ny, ensPx%lman%nz))
        allocate(f3d_sum(ensPx%lman%nx, ensPx%lman%ny, ensPx%lman%nz))
        allocate(f3d_sum2(ensPx%lman%nx, ensPx%lman%ny, ensPx%lman%mz))
        allocate(f3d_sum3(ensPx%lman%nx, ensPx%lman%ny, ensPx%lman%mz))
        n_ens=ensPx%n_mem
        px=ensPx%lman%x3dm(1,1)%m
        py=ensPx%lman%y2dm(1)%m
        pz=ensPx%lman%zm%m
        allocate(vct_mtrx(pz,py,px,n_ens))
        vct_mtrx=reshape(vct, (/pz,py,px,n_ens/))
        nx_big=maxval(md%xcoords(:)%n)
        ny_big=maxval(md%ycoords(:)%n)
        nz_big=maxval(md%zcoords(:)%n)
        allocate(xweight(nx_big,2,md%n_xcoord), &
                 yweight(ny_big,2,md%n_ycoord), &
                 zweight(nz_big,2,md%n_zcoord), &
                 xneighbor(nx_big,2,md%n_xcoord), &
                 yneighbor(ny_big,2,md%n_ycoord), &
                 zneighbor(nz_big,2,md%n_zcoord))
        do i=1, md%n_xcoord
            call interp_1d_linear_nmg2sct_getWeight(ensPx%lman%xcoord,md%xcoords(i)%val, &
                   xneighbor(1:md%xcoords(i)%n,:,i), xweight(1:md%xcoords(i)%n,:,i))
        end do
        do i=1, md%n_ycoord
            call interp_1d_linear_nmg2sct_getWeight(ensPx%lman%ycoord,md%ycoords(i)%val, &
                   yneighbor(1:md%ycoords(i)%n,:,i), yweight(1:md%ycoords(i)%n,:,i))
        end do
        do i=1, md%n_zcoord
            call interp_1d_linear_nmg2sct_getWeight(ensPx%lman%zcoord,md%zcoords(i)%val, &
                   zneighbor(1:md%zcoords(i)%n,:,i), zweight(1:md%zcoords(i)%n,:,i))
        end do
        allocate(yweight_m2n(ensPx%lman%ny,2), yneighbor_m2n(ensPx%lman%ny, 2))
        allocate(zweight_m2n(ensPx%lman%nz,2), zneighbor_m2n(ensPx%lman%nz, 2))
        call localman_get_my2ny_Weight(ensPx%lman, yneighbor_m2n, yweight_m2n)
        call localman_get_mz2nz_Weight(ensPx%lman, zneighbor_m2n, zweight_m2n)
        TIME_CLICK(timer, timer_bgn, timer_last) 
        !! test neighbor/weight: passed {{{
      !  call write_integer_matrix2(xneighbor(:,:,1), fname="xneighbor.csv", &
      !          element_fmat= "2(I,',')") 
      !  call write_float_matrix2(xweight(:,:,1), fname="xweight.csv", &
      !          element_fmat= "2(F,',')")
      !  call write_integer_matrix2(yneighbor(:,:,1), fname="yneighbor.csv", &
      !          element_fmat= "2(I,',')") 
      !  call write_float_matrix2(yweight(:,:,1), fname="yweight.csv", &
      !          element_fmat= "2(F,',')")
      !  call write_integer_matrix2(zneighbor(:,:,1), fname="zneighbor.csv", &
      !          element_fmat= "2(I,',')") 
      !  call write_float_matrix2(zweight(:,:,1), fname="zweight.csv", &
      !          element_fmat= "2(F,',')")
      !  stop
      !   call write_integer_matrix2(yneighbor_m2n, fname="yneighbor_m2n.csv", &
      !          element_fmat="2(I,',')")
      !   call write_float_matrix2(yweight_m2n, fname="yweight_m2n.csv", &
      !          element_fmat="2(F,',')")
      !   call write_integer_matrix2(zneighbor_m2n, fname="zneighbor_m2n.csv", &
      !          element_fmat="2(I,',')")
      !   call write_float_matrix2(zweight_m2n, fname="zweight_m2n.csv", &
      !          element_fmat="2(F,',')") }}}
        
        out_vct=0
        do s=1, n_ens
            print*, "test: cycle for ensemble member:", s
            print*, "-test: before diag(Rouz)diag(Rouy)Roux*beta(j,k,s)"
            f3d_sum=0
            f3d_sum2=0
            f3d_sum3=0
            do k=1, pz
                !print*,"k=",k
                do j=1, py
                    call localman_xm_multiply_vector(ensPx%lman, vct_mtrx(k,j,:,s), f3d)
                    call PxMV_interp1(f3d, f3d2, yneighbor_m2n, yweight_m2n)
                    f3d_sum2=f3d_sum2+f3d2
                    call localman_ymi_schur_f3d(ensPx%lman, j, f3d2)
                    f3d_sum3=f3d_sum3+f3d2
                    call PxMV_interp2(f3d2, f3d3, zneighbor_m2n, zweight_m2n)
                    call localman_zmi_schur_f3d(ensPx%lman, k, f3d3)
                    f3d_sum=f3d_sum+f3d3
                end do
            end do
            TIME_CLICK(timer, timer_bgn, timer_last) 
            !! output f3d_sum* 
      !      module_ens_cnt=module_ens_cnt+1
      !      call write_float_matrix2(f3d_sum(:,:,1), fname="f3d_sumA1_"//trim(tostring(ival=module_ens_cnt))//".csv", &
      !                 element_fmat=  "91(F,',')")
      !      call write_float_matrix2(f3d_sum(:,:,2), fname="f3d_sumA2_"//trim(tostring(ival=module_ens_cnt))//".csv", &
      !                 element_fmat=  "91(F,',')")
      !      call write_float_matrix2(f3d_sum2(:,:,1), fname="f3d_sumB"//trim(tostring(ival=module_ens_cnt))//".csv", &
      !                 element_fmat=  "91(F,',')")
      !      call write_float_matrix2(f3d_sum3(:,:,1), fname="f3d_sumC"//trim(tostring(ival=module_ens_cnt))//".csv", &
      !                 element_fmat=  "91(F,',')")  
            !!interp 
            print*, "-test: before interp from localization grid to model grid"
            do ivar=1, md%nvar
                ix=md%vars(ivar)%xcoord_id
                iy=md%vars(ivar)%ycoord_id
                iz=md%vars(ivar)%zcoord_id
                nx_hr=md%xcoords(ix)%n
                ny_hr=md%ycoords(iy)%n
                nz_hr=md%zcoords(iz)%n
                allocate(f3d_hr(nx_hr,ny_hr,nz_hr), f3d_Px(nx_hr,ny_hr,nz_hr))
                call interp_3d_withWeight(f3d_sum, f3d_hr, &
                   xneighbor(1:nx_hr,:,ix),xweight(1:nx_hr,:,ix),  &
                   yneighbor(1:ny_hr,:,iy),yweight(1:ny_hr,:,iy),  &
                   zneighbor(1:nz_hr,:,iz),zweight(1:nz_hr,:,iz))
                call modelstate_get3dfield(md, ivar, ensPx%Px(:,s), f3d_Px)
                f3d_Px=f3d_Px * f3d_hr
                out_vct(md%vars(ivar)%ibgn:md%vars(ivar)%iend)= &
                       out_vct(md%vars(ivar)%ibgn:md%vars(ivar)%iend) + &
                       reshape(f3d_Px,(/nx_hr * ny_hr * nz_hr/))
                deallocate(f3d_hr, f3d_Px)
            end do
            TIME_CLICK(timer, timer_bgn, timer_last) 
        end do
        deallocate(f3d, f3d2, f3d3, f3d_sum)
        deallocate(vct_mtrx)
        deallocate(xweight, xneighbor, yweight, yneighbor, zweight, zneighbor)
        deallocate(yweight_m2n, yneighbor_m2n, zweight_m2n, zneighbor_m2n)
    contains 
        !! f3d, from (nx,my,mz) -> (nx,ny,mz)
        subroutine PxMV_interp1(f3d_in, f3d_out, yneighbors, yweights)  !!{{{
        implicit none
            character(len=*),parameter :: PROCEDURE_NAME="PxMV_interp1"
            integer :: yneighbors(:,:)
            real :: yweights(:,:), f3d_in(:,:,:), f3d_out(:,:,:)
            integer :: nx, ny, mz, my, i, j, k
            nx=size(f3d_in,1)
            my=size(f3d_in,2)
            mz=size(f3d_in,3)
            ny=size(f3d_out,2)
            if(.not. all(shape(yneighbors) == (/ny,2/))) then
                RP_ERROR("input yneighbors size wrong!")
            end if
            if(.not. all(shape(yweights) == (/ny,2/))) then
                RP_ERROR("input yweights size wrong!")
            end if
            if(.not. all(shape(f3d_out) == (/nx,ny,mz/)))then
                RP_ERROR("input f3d_out shape wrong!")
            end if
            do j=1, ny
        !        print*,j,yneighbors(j,1),yweights(j,1),yneighbors(j,2),yweights(j,2)
                f3d_out(:,j,:)=f3d_in(:,yneighbors(j,1),:)*yweights(j,1) + &
                               f3d_in(:,yneighbors(j,2),:)*yweights(j,2)
        !        print*, f3d_in(100,yneighbors(j,1),1), f3d_in(100,yneighbors(j,2),1), &
        !                f3d_out(100,j,1)
            end do
        end subroutine PxMV_interp1  !!}}}
        !! f3d, from (nx,ny,mz) -> (nx,ny,nz)
        subroutine PxMV_interp2(f3d_in, f3d_out, zneighbors, zweights)  !!{{{
        implicit none
            character(len=*),parameter :: PROCEDURE_NAME="PxMV_interp2"
            integer :: zneighbors(:,:)
            real :: zweights(:,:), f3d_in(:,:,:), f3d_out(:,:,:)
            integer :: nx, ny, mz, nz, i, j, k
            nx=size(f3d_in,1)
            ny=size(f3d_in,2)
            mz=size(f3d_in,3)
            nz=size(f3d_out,3)
            if(.not. all(shape(zneighbors) == (/nz,2/))) then
                RP_ERROR("input yneighbors size wrong!")
            end if
            if(.not. all(shape(zweights) == (/nz,2/))) then
                RP_ERROR("input yweights size wrong!")
            end if
            if(.not. all(shape(f3d_out) == (/nx,ny,nz/)))then
                RP_ERROR("input f3d_out shape wrong!")
            end if
            do k=1, nz
          !      print*,k,zneighbors(k,1),zweights(k,1),zneighbors(k,2),zweights(k,2)
                f3d_out(:,:,k)=f3d_in(:,:,zneighbors(k,1))*zweights(k,1) + &
                               f3d_in(:,:,zneighbors(k,2))*zweights(k,2)
            end do
        end subroutine PxMV_interp2  !!}}}
    end subroutine Px_multiply_vector  !!}}}

    subroutine Px_get_row(md, ensPx, ivar, ix, iy, iz, row)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="Px_get_row"
        type(type_modeldesc) :: md
        type(type_ensPx) :: ensPx
        integer :: ivar, ix, iy, iz
        real :: row(:), x1, y1, z1
        integer :: nrow
        real,allocatable :: row_mtrx(:,:,:,:)
        type(type_localmodes) :: xlm, ylm, zlm
        integer :: n1, n2, n3, n4, xid, yid, zid, i, j, k, t
        integer :: igrid, nx, ny, nz, ibgn, xlm_idx
        nrow=size(row)
        ASSUREX(ensPx%has_localize)
        ASSUREX(size(row) == ensPx%n_total_mem)
        ASSUREX((ivar>=1) .and. (ivar<=md%nvar))
        xid=md%vars(ivar)%xcoord_id
        yid=md%vars(ivar)%ycoord_id
        zid=md%vars(ivar)%zcoord_id
        nx= md%xcoords(xid)%n
        ny= md%ycoords(yid)%n
        nz= md%zcoords(zid)%n
        ASSUREX((ix<=nx) .and. (ix>=1))
        ASSUREX((iy<=ny) .and. (iy>=1))
        ASSUREX((iz<=nz) .and. (iz>=1))
        x1=md%xcoords(xid)%val(ix)
        y1=md%ycoords(yid)%val(iy)
        z1=md%zcoords(zid)%val(iz)
        call localman_interp_toscatter(ensPx%lman, (/x1/), (/y1/), (/z1/), xlm, ylm, zlm)
      !  call localmodes_rownormalize(xlm)
      !  call localmodes_rownormalize(ylm)
      !  call localmodes_rownormalize(zlm)
        n1=ensPx%n_mem
        n2=ensPx%lman%x3dm(1,1)%m
        n3=ensPx%lman%y2dm(1)%m
        n4=ensPx%lman%zm%m
        allocate(row_mtrx(n4,n3,n2,n1))
        igrid=(iz-1)*(nx*ny)+(iy-1)*nx+ix
        ibgn=md%vars(ivar)%ibgn
        do t=1, n1
            do i=1, n2
                do j=1, n3
                    do k=1, n4
                        row_mtrx(k,j,i,t)=ensPx%Px(ibgn-1+igrid,t)* &
                                   xlm%modes(1,i)*ylm%modes(1,j)*zlm%modes(1,k)
                    end do
                end do
            end do
        end do
        row=reshape(row_mtrx, (/nrow/))
        deallocate(row_mtrx)
        call localmodes_delete(xlm)
        call localmodes_delete(ylm)
        call localmodes_delete(zlm)
    end subroutine Px_get_row  !!}}}

    subroutine Py_get_row(od, ensPy, iobs, row)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="Py_get_row"
        type(type_obsdesc) :: od
        type(type_ensPy) :: ensPy
        real :: row(:)
        integer :: nrow
        real,allocatable :: row_mtrx(:,:,:,:)
        integer :: n1, n2, n3, n4, i, j, k, t, iobs
        nrow=size(row)
        if(.not. ensPy%has_localize) then
            RP_ERROR( "ensPy not has localized")
        end if
        if(nrow /= ensPy%n_total_mem ) then
            RP_ERROR( "size(row) /= ensPy%n_total_mem" )
        end if
        if(iobs > od%n_obs) then
            RP_ERROR( "iobs > od%n_obs" )
        end if
        n1=ensPy%n_mem
        n2=ensPy%xlm%m
        n3=ensPy%ylm%m
        n4=ensPy%zlm%m
        allocate(row_mtrx(n4,n3,n2,n1))
        do t=1, n1
            do i=1, n2
                do j=1, n3
                    do k=1, n4
                        row_mtrx(k,j,i,t)=ensPy%Py(iobs,t) &
                                         *ensPy%xlm%modes(iobs,i) &
                                         *ensPy%ylm%modes(iobs,j) &
                                         *ensPy%zlm%modes(iobs,k)
                    end do
                end do
            end do
        end do
        row=reshape(row_mtrx, (/nrow/)) 
        deallocate(row_mtrx)
    end subroutine Py_get_row  !!}}}

    !! sqrt(diag(Px*Px'))
    subroutine Px_get_stdv(md, ensPx, stdv)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="Px_get_stdv"
        type(type_modeldesc) :: md
        type(type_ensPx) :: ensPx
        real :: stdv(:)
        real,allocatable :: row(:)
        integer :: n_stdv
        integer :: ivar, xid, yid, zid, igrid 
        integer :: xidx, yidx, zidx, nx, ny, nz, ibgn
        if(.not. ensPx%has_localize) then
            RP_ERROR("ensPx not localized")
        end if
        n_stdv=size(stdv)
        if(n_stdv /= ensPx%n_state) then
            RP_ERROR("size(stdv) /= ensPx%n_state")
        end if
        allocate(row(ensPx%n_total_mem))
        do ivar=1, md%nvar
            xid=md%vars(ivar)%xcoord_id
            yid=md%vars(ivar)%ycoord_id
            zid=md%vars(ivar)%zcoord_id
            nx= md%xcoords(xid)%n
            ny= md%ycoords(yid)%n
            nz= md%zcoords(zid)%n
            ibgn=md%vars(ivar)%ibgn
            do igrid=1, md%vars(ivar)%n_grids
                zidx=ceiling(igrid*1./(nx*ny))
                yidx=ceiling(1.*(igrid-(zidx-1)*(nx*ny))/nx)
                xidx=igrid-(zidx-1)*(nx*ny)-(yidx-1)*nx
                call Px_get_row(md, ensPx, ivar, xidx, yidx, zidx, row)
                stdv(ibgn-1+igrid)=sqrt(sum(row**2))
            end do
        end do
        deallocate(row)
    end subroutine Px_get_stdv  !!}}}

    
!!operate on ensPx%Px, ensPy%Py, i.e., the dynamic ensemble (not on localization modes)
    subroutine Px_full2pert(ensPx, xcenter)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="Px_full2pert"
        type(type_ensPx) :: ensPx
        real :: xcenter(:)
        integer :: i
        if(size(xcenter) /= ensPx%n_state) then
            RP_ERROR("size(xcenter) /= ensPx%n_state")
        end if
        do i=1, ensPx%n_mem
            ensPx%Px(:,i)=ensPx%Px(:,i)-xcenter
        end do
    end subroutine Px_full2pert  !!}}}

    subroutine Py_full2pert(ensPy, ycenter)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="Py_full2pert"
        type(type_ensPy) :: ensPy
        real :: ycenter(:)
        integer :: i
        if(size(ycenter) /= ensPy%n_state) then
            RP_ERROR("size(ycenter) /= ensPy%n_state")
        end if
        do i=1, ensPy%n_mem
            ensPy%Py(:,i)=ensPy%Py(:,i)-ycenter
        end do
    end subroutine Py_full2pert  !!}}}

    subroutine Px_pert2full(ensPx, xcenter)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="Px_pert2full"
        type(type_ensPx) :: ensPx
        real :: xcenter(:)
        integer :: i
        if(size(xcenter) /= ensPx%n_state) then
            RP_ERROR("size(xcenter) /= ensPx%n_state")
        end if
        do i=1, ensPx%n_mem
            ensPx%Px(:,i)=ensPx%Px(:,i)+xcenter
        end do
    end subroutine Px_pert2full  !!}}}

    subroutine Px_average(ensPx, xmean)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="Px_average"
        type(type_ensPx) :: ensPx
        real :: xmean(:)
        integer :: i
        xmean=0
        do i=1, ensPx%n_mem
            xmean=xmean+ensPx%Px(:,i)
        end do
        xmean=xmean/ensPx%n_mem
    end subroutine Px_average  !!}}}

    !!similar with Px_get_stdv, but only operates on ensPx%Px (the dynamic ensemble)
    !!it should be same with result from Px_get_stdv if localmodes have been rownormalized
    !!with much little computation
    subroutine Px_get_stdv_px(ensPx, stdv)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="Px_get_stdv_px"
        type(type_modeldesc) :: md
        type(type_ensPx)  :: ensPx
        real :: stdv(:)
        integer :: i
        if(size(stdv) /= ensPx%n_state) then
            RP_ERROR("size(stdv) /= ensPx%n_state")
        end if
        do i=1, ensPx%n_state
            stdv(i)=sqrt(sum(ensPx%Px(i,:)**2))
        end do
    end subroutine Px_get_stdv_px  !!}}}

!!delete
    subroutine Px_delete(ensPx)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="Px_delete"
        type(type_ensPx) :: ensPx
        deallocate(ensPx%Px)
        if(ensPx%has_localize) then
            call localman_dealloc(ensPx%lman)
            ensPx%has_localize=.false.
        end if
        ensPx%n_mem=0
        ensPx%n_state=0
        ensPx%n_total_mem=0
    end subroutine Px_delete  !!}}}

    subroutine Py_delete(ensPy)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="Py_delete(ensPy)"
        type(type_ensPy) :: ensPy
        deallocate(ensPy%Py)
        if(ensPy%has_localize) then
            call localmodes_delete(ensPy%xlm)
            call localmodes_delete(ensPy%ylm)
            call localmodes_delete(ensPy%zlm)
            ensPy%has_localize=.false.
        end if
        ensPy%n_mem=0
        ensPy%n_state=0
        ensPy%n_total_mem=0
    end subroutine Py_delete  !!}}}


end module module_ensemble
