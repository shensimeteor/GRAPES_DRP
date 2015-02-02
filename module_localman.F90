! --------------------------------------------------------------------------
! Description:
!   NOW ASSUME :: x/y must be uniform(dx/dy constant),increasing order; z in any monotonous order
!
!   localman_set_domain/localman_input_domain :: set or input domain (nx,ny,nz,xcoord,ycoord,zcoord) &
!                                             even all set_*_nolocal, domain still be set/read
!   localman_set_xm/localman_input_xm/localman_set_xm_nolocal  ::  Not all of these 3*3 implemented, &
!   localman_set_ym/localman_input_ym/localman_set_ym_nolocal  ::  XM/YM must have same mz, must have same is_*_loc value, &
!   localman_set_zm/localman_input_zm/localman_set_zm_nolocal  ::  will NOT do rownormalize; Must call set/read_xm before set/read_ym
!   localman_output_xm/localman_output_ym/localman_output_zm
!   localman_interp_toscatter :: from lman --> xlm/ylm/zlm for speficied points series [xpos,ypos,zpos]
!   localman_rownormalize :: do row_normalize for xm/ym/zm
!   localman_xm_multiply_vector :: ROUx * vector(px) -> f3d(nx, my, mz)
!   localman_ymi_schur_f3d :: diag(ROUy_ipy) * f3d(nx,ny,mz) -> f3d
!   localman_zmi_schur_f3d :: diag(ROUz_ipz) * f3d(nx,ny,nz) -> f3d
!   localman_dealloc
!
! History:
!
!   2015-01-06: 
!
! --------------------------------------------------------------------------
#include "inc_common.fpp"
module module_localman
use utility
use module_types
use module_local
use netcdf
implicit none

    ! variables

contains
    ! methods
    subroutine localman_set_domain(lman, xcoord, ycoord, zcoord, is_x_cycle)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="localman_set_domain"
        real :: xcoord(:), ycoord(:), zcoord(:)
        logical :: is_x_cycle
        type(type_localman) :: lman
        if(lman%has_setdomain) then
            RP_ERROR("lman has set domained already!")
        end if
        lman%nx=size(xcoord)
        lman%ny=size(ycoord)
        lman%nz=size(zcoord)
        lman%is_x_cycle=is_x_cycle
        allocate(lman%xcoord(lman%nx))
        lman%xcoord=xcoord
        allocate(lman%ycoord(lman%ny))
        lman%ycoord=ycoord
        allocate(lman%zcoord(lman%nz))
        lman%zcoord=zcoord
        lman%has_setdomain=.true.
    end subroutine localman_set_domain  !!}}}

    !!set x modes ; if xl0(i)<0, means the xl0(i)=infinite, i.e.: the pole; x muse be cycle
    subroutine localman_set_xm(lman, ybgn, yend, zbgn, zend,  xl0, n_mode)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="localman_set_xm"
        type(type_localman) :: lman
        integer :: ybgn(:), yend(:), zbgn(:), zend(:), n_mode
        integer :: my, mz
        real :: xl0(:,:), dx, w1, w2
        integer :: i,j,k,s
        if(.not. lman%has_setdomain) then
            RP_ERROR("lman has NOT set domain yet!")
        end if
        if( lman%has_set_xm) then
            RP_ERROR("lman has set xmodes already!")
        end if
        if( .not. all((/size(ybgn), size(yend)/)==size(xl0,1))) then
            RP_ERROR("input parameter size conflicts")
        end if
        if( .not. all((/size(zbgn), size(zend)/)==size(xl0,2))) then
            RP_ERROR("input parameter size conflicts")
        end if
        my=size(ybgn)
        mz=size(zbgn)
        lman%my=my
        lman%mz=mz
        allocate(lman%ybgn(my), lman%yend(my))
        allocate(lman%zbgn(mz), lman%zend(mz))
        lman%ybgn=ybgn
        lman%yend=yend
        lman%zbgn=zbgn
        lman%zend=zend
        allocate(lman%x3dm(my,mz))
        dx=lman%xcoord(2)-lman%xcoord(1)
        do k=1, mz
            do j=1, my
                if(xl0(j,k)<0) then
                    call localmodes_generate_Fourier_pole(lman%x3dm(j,k), lman%nx, n_mode)
                else
                    call localmodes_generate_Fourier(lman%x3dm(j,k), lman%nx, dx, xl0(j,k), 1, n_mode)
                end if
            end do
        end do
        lman%has_set_xm=.true.
        lman%is_x_loc=.true.
    end subroutine localman_set_xm  !!}}}

    subroutine localman_set_xm_nolocal(lman)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="localman_set_xm_nolocal"
        type(type_localman) :: lman
        lman%my=1
        lman%mz=1
        allocate(lman%ybgn(1),lman%yend(1))
        allocate(lman%zbgn(1),lman%zend(1))
        lman%ybgn=1
        lman%yend=lman%ny
        lman%zbgn=1
        lman%zend=lman%nz
        allocate(lman%x3dm(1,1))
        call localmodes_generate_nolocal(lman%x3dm(1,1), lman%nx)
        lman%has_set_xm=.true.
        lman%is_x_loc=.false.
    end subroutine localman_set_xm_nolocal  !!}}}

    !! by Fourier also : results wrong: Become cyclic Correlation Matrix
    subroutine localman_set_ym(lman, zbgn, zend, yl0, n_mode)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="localman_set_ym"
        type(type_localman) :: lman
        integer :: zbgn(:), zend(:)
        real :: yl0(:), dy, w1, w2
        integer ::  n_mode, k, mz, j
        if(.not. lman%has_setdomain) then
            RP_ERROR("should call localman_set_domain before!")
        end if
        if( lman%has_set_ym) then
            RP_ERROR("lman has set ymodes already!")
        end if
        if(.not. lman%has_set_xm) then
            RP_ERROR("lman should set xm first!")
        end if
        if(.not. all((/size(zbgn), size(zend), size(yl0)/) == lman%mz)) then
            RP_ERROR("input parameter size conflicts")
        end if
        if(.not. all( zbgn == lman%zbgn ) ) then
            RP_ERROR("input parameter zbgn conflicts with lman%zbgn")
        end if
        if(.not. all( zend == lman%zend ) ) then
            RP_ERROR("input parameter zend conflicts with lman%zbgn")
        end if
        mz=lman%mz
        dy=lman%ycoord(2)-lman%ycoord(1)
        allocate(lman%y2dm(mz))
        do k=1, mz
            call localmodes_generate_Fourier(lman%y2dm(k), lman%ny, dy, yl0(k), 1, n_mode) !! sin?
        end do
        lman%has_set_ym=.true.
        lman%is_y_loc=.true.
    end subroutine localman_set_ym  !!}}}

    !! only for vertically uniformed  (mz=1)
    subroutine localman_set_ym_byEIG(lman, yl0, n_mode)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="localman_set_ym_byEIG"
        type(type_localman) :: lman
        real :: yl0, dy, w1, w2
        integer ::  n_mode, k, mz, j
        if(.not. lman%has_setdomain) then
            RP_ERROR("should call localman_set_domain before!")
        end if
        if( lman%has_set_ym) then
            RP_ERROR("lman has set ymodes already!")
        end if
        if(.not. lman%has_set_xm) then
            RP_ERROR("lman should set xm first!")
        end if
        if(lman%mz /= 1) then
            RP_ERROR("lman.mz not equal to 1, can NOT use set_ym_byEIG")
        end if
        allocate(lman%y2dm(1))
        call localmodes_generate_decompC(lman%y2dm(1), lman%ycoord, yl0, 1, .false., 0., op_nmodes=n_mode)
        lman%has_set_ym=.true.
        lman%is_y_loc=.true.
    end subroutine localman_set_ym_byEIG  !!}}}

    subroutine localman_set_ym_nolocal(lman)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="localman_set_ym_nolocal"
        type(type_localman) :: lman
        if(.not. lman%has_set_xm) then
            RP_ERROR("lman should set xm first!")
        end if
        if(lman%is_x_loc) then
            RP_ERROR("lman must have is_x_loc/is_y_loc same logical value!")
        end if
        allocate(lman%y2dm(1))
        call localmodes_generate_nolocal(lman%y2dm(1), lman%ny)
        lman%has_set_ym=.true.
        lman%is_y_loc=.false.
    end subroutine localman_set_ym_nolocal  !!}}}

    subroutine localman_set_zm_nolocal(lman)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="localman_set_zm_nolocal"
        type(type_localman) :: lman
        if(.not. lman%has_setdomain) then
            RP_ERROR("should call localman_set_domain before!")
        end if
        if( lman%has_set_zm) then
            RP_ERROR("lman has set zmode already!")
        end if
        call localmodes_generate_nolocal(lman%zm, lman%nz)
        lman%has_set_zm=.true.
        lman%is_z_loc=.false.
    end subroutine localman_set_zm_nolocal  !!}}}

    subroutine localman_dealloc(lman)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="localman_dealloc"
        type(type_localman) :: lman
        if(lman%has_setdomain) then
            deallocate(lman%xcoord, lman%ycoord, lman%zcoord)
            lman%nx=0
            lman%ny=0
            lman%nz=0
        end if
        if(lman%has_set_xm) then
            deallocate(lman%ybgn, lman%yend)
            deallocate(lman%zbgn, lman%zend)
            deallocate(lman%x3dm)
            lman%has_set_xm=.false.
            lman%is_x_loc=.false.
        end if
        if(lman%has_set_ym) then
            deallocate(lman%y2dm)
            lman%has_set_ym=.false.
            lman%is_y_loc=.false.
        end if
        if(lman%has_set_zm) then
            call localmodes_delete(lman%zm)
            lman%has_set_zm=.false.
            lman%is_z_loc=.false.
        end if
        lman%has_setdomain=.false.
    end subroutine localman_dealloc  !!}}}

    !!use lman to generate observation localmodes xlm/ylm/zlm
    subroutine localman_interp_toscatter(lman, xpos, ypos, zpos, xlm, ylm, zlm)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="localman_interp_toscatter"
        type(type_localman) :: lman
        real :: xpos(:), ypos(:), zpos(:)
        type(type_localmodes) :: xlm, ylm, zlm, xm1, xm2, ym1, ym2, xm3, xm4
        integer :: nsct, j, k, y1, y2, z1, z2, isct
        real :: dy, yw1, yw2, zw1, zw2
        logical :: found
        ASSUREX(lman%has_set_xm)
        ASSUREX(lman%has_set_ym)
        ASSUREX(lman%has_set_zm)
        nsct=size(xpos)
        if(.not. all(nsct == (/size(ypos), size(zpos)/))) then
            RP_ERROR("input xpos/ypos/zpos size conflict")
        end if
        dy=lman%ycoord(2)-lman%ycoord(1)
        !! for zlm
        call localmodes_interp_linear(lman%zcoord, lman%zm, zpos, zlm, .false.)
        !! for xlm/ylm
        call localmodes_init(xlm, nsct, lman%x3dm(1,1)%m)
        call localmodes_init(ylm, nsct, lman%y2dm(1)%m)
        do isct=1, nsct
        !!for xm->xlm
            found=.false.
            !!to find z1/z2/zw1/zw2
            do k=1, lman%mz
                if( (zpos(isct)-lman%zcoord(lman%zbgn(k))) * &
                    (zpos(isct)-lman%zcoord(lman%zend(k))) <= 0) then
                    z1=k
                    z2=k
                    zw1=1.
                    zw2=0.
                    found=.true.
                    exit
                end if
            end do
            if(.not. found) then
                do k=1, lman%mz-1
                    if( (zpos(isct)-lman%zcoord(lman%zend(k))) * &
                        (zpos(isct)-lman%zcoord(lman%zbgn(k+1))) <0 ) then
                        z1=k
                        z2=k+1
                        zw1=(lman%zcoord(lman%zbgn(k+1))-zpos(isct)) / &
                            (lman%zcoord(lman%zbgn(k+1)) - lman%zcoord(lman%zend(k)))
                        zw2=1-zw1
                        ASSUREX(zw1>0)
                        ASSUREX(zw2>0)
                        found=.true.
                        exit
                    end if
                end do
            end if
            ASSUREX(found)
            !!to find y1/y2/yw1/yw2
            found=.false.
            do j=1, lman%my
                if( (ypos(isct)-lman%ycoord(lman%ybgn(k))) * &
                    (ypos(isct)-lman%ycoord(lman%yend(k))) <= 0) then
                    y1=j
                    y2=j
                    yw1=1.
                    yw2=0.
                    found=.true.
                    exit
                end if
            end do
            if(.not. found) then
                do j=1, lman%my-1
                    if( (ypos(isct)-lman%ycoord(lman%yend(j))) * &
                        (ypos(isct)-lman%ycoord(lman%ybgn(j+1))) <0 ) then
                        y1=j
                        y2=j+1
                        yw1=(lman%ycoord(lman%ybgn(j+1))-ypos(isct)) / &
                            (lman%ycoord(lman%ybgn(j+1))-lman%ycoord(lman%yend(j)))
                        yw2=1-yw1
                        ASSUREX(yw1>0)
                        ASSUREX(yw2>0)
                        found=.true.
                        exit
                    end if
                end do
            end if
            ASSUREX(found)
            !!combine linearly & interp
            if(yw2 == 0) then
                call localmodes_init(xm1, lman%nx, lman%x3dm(1,1)%m)
                call localmodes_init(xm2, lman%nx, lman%x3dm(1,1)%m)
                xm1=lman%x3dm(y1,z1)
                xm2=lman%x3dm(y1,z2)
            else
                call localmodes_linear_combine(lman%x3dm(y1,z1), yw1, lman%x3dm(y2,z1), yw2, xm1)
                call localmodes_linear_combine(lman%x3dm(y1,z2), yw1, lman%x3dm(y2,z2), yw2, xm2)
            end if
            if(zw2 == 0) then
                call localmodes_init(xm3, lman%nx, lman%x3dm(1,1)%m)
                xm3=xm1
            else
                call localmodes_linear_combine(xm1, zw1, xm2, zw2, xm3)
            end if !!combination finish
            call localmodes_interp_linear(lman%xcoord, xm3, (/xpos(isct)/), xm4, lman%is_x_cycle)
            xlm%modes(isct,:)=xm4%modes(1,:) !!done xlm for 3D models
            call localmodes_delete(xm1)
            call localmodes_delete(xm2)
            call localmodes_delete(xm3)
            call localmodes_delete(xm4)
        !! for ym->ylm, using z1/z2/zw1/zw2
            if(zw2 == 0) then
                call localmodes_init(ym1, lman%ny, lman%y2dm(z1)%m)
                ym1=lman%y2dm(z1)
            else
                call localmodes_linear_combine(lman%y2dm(z1), zw1, lman%y2dm(z2), zw2, ym1)
            end if
            call localmodes_interp_linear(lman%ycoord, ym1, (/ypos(isct)/), ym2, .false.)
            ylm%modes(isct,:)=ym2%modes(1,:)
            call localmodes_delete(ym1)
            call localmodes_delete(ym2)
        end do
    end subroutine localman_interp_toscatter  !!}}}

    !!model dependent: for nolocal: when call, file_xm can be substitued by file_ym/file_zm
    subroutine localman_input_domain(lman, file_xm)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="localman_input_domain"
        character(len=*) :: file_xm
        type(type_localman) :: lman
        integer :: ncid, lonid, latid, levid, pxid, xmid, &
                   lon_vid, lat_vid, lev_vid, px_vid, i, j, k, &
                   nx, ny, nz, px
        real,allocatable :: lon(:), lat(:), lev(:)
        integer :: int_cycle
        logical :: is_cycle
        NC_CHECK(nf90_open(trim(file_xm), NF90_NOWRITE, ncid))
        NC_CHECK(nf90_inq_dimid(ncid, "lon", lonid))
        NC_CHECK(nf90_inquire_dimension(ncid, lonid, len=nx)) !!nx
        NC_CHECK(nf90_inq_dimid(ncid, "lat", latid))
        NC_CHECK(nf90_inquire_dimension(ncid, latid, len=ny)) !!ny
        NC_CHECK(nf90_inq_dimid(ncid, "lev", levid))  
        NC_CHECK(nf90_inquire_dimension(ncid, levid, len=nz)) !!nz
        if(.not. lman%has_setdomain) then
            allocate(lon(nx), lat(ny), lev(nz))
            NC_CHECK(nf90_inq_varid(ncid, "lon", lon_vid))
            NC_CHECK(nf90_get_att(ncid, lon_vid, "is_cycle",  int_cycle))
            if(int_cycle == 1) then
                is_cycle=.true.
            else
                is_cycle=.false.
            end if
            NC_CHECK(nf90_get_var(ncid, lon_vid, lon))
            NC_CHECK(nf90_inq_varid(ncid, "lat", lat_vid))
            NC_CHECK(nf90_get_var(ncid, lat_vid, lat))
            NC_CHECK(nf90_inq_varid(ncid, "lev", lev_vid))
            NC_CHECK(nf90_get_var(ncid, lev_vid, lev))
            call localman_set_domain(lman, lon, lat, lev, is_cycle)
            deallocate(lon, lat, lev)
        end if
    end subroutine localman_input_domain  !!}}}

    !! model dependent
    subroutine localman_output_xm(lman, file_xm, string)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="localman_output_xm"
        character(len=*) :: file_xm
        type(type_localman) :: lman
        integer :: ncid, lonid, latid, levid, pxid, xmid, &
                   lon_vid, lat_vid, lev_vid, px_vid, i, j, k, int_cycle, &
                   myid, mzid, ybgn_vid, yend_vid, zbgn_vid, zend_vid
        integer,allocatable :: px_idx(:)
        real,allocatable :: xmodes(:,:,:,:)
        character(len=*) :: string
        ASSUREX(lman%has_setdomain)
        ASSUREX(lman%has_set_xm)
        NC_CHECK(nf90_create(trim(file_xm), NF90_CLOBBER, ncid))
        NC_CHECK(nf90_def_dim(ncid, "lon", lman%nx, lonid))
        NC_CHECK(nf90_def_dim(ncid, "lat", lman%ny, latid))
        NC_CHECK(nf90_def_dim(ncid, "lev", lman%nz, levid))
        NC_CHECK(nf90_def_dim(ncid, "px", lman%x3dm(1,1)%m, pxid))
        NC_CHECK(nf90_def_dim(ncid, "my", lman%my, myid))
        NC_CHECK(nf90_def_dim(ncid, "mz", lman%mz, mzid))
        NC_CHECK(nf90_def_var(ncid, "lon", NF90_FLOAT, (/lonid/), lon_vid))
        NC_CHECK(nf90_def_var(ncid, "lat", NF90_FLOAT, (/latid/), lat_vid))
        NC_CHECK(nf90_def_var(ncid, "lev", NF90_FLOAT, (/levid/), lev_vid))
        NC_CHECK(nf90_def_var(ncid, "px", NF90_INT,(/pxid/), px_vid))
        NC_CHECK(nf90_def_var(ncid, "ybgn", NF90_INT,(/myid/), ybgn_vid))
        NC_CHECK(nf90_def_var(ncid, "yend", NF90_INT,(/myid/), yend_vid))
        NC_CHECK(nf90_def_var(ncid, "zbgn", NF90_INT,(/mzid/), zbgn_vid))
        NC_CHECK(nf90_def_var(ncid, "zend", NF90_INT,(/mzid/), zend_vid))
        NC_CHECK(nf90_def_var(ncid, "xmode", NF90_FLOAT, (/lonid,myid,mzid,pxid/), xmid))
        NC_CHECK(nf90_put_att(ncid, lon_vid, "units", "degrees_east"))
        if(lman%is_x_cycle) then
            int_cycle=1
        else
            int_cycle=0
        end if
        NC_CHECK(nf90_put_att(ncid, lon_vid, "is_cycle", int_cycle))
        NC_CHECK(nf90_put_att(ncid, lat_vid, "units", "degrees_north"))
        NC_CHECK(nf90_put_att(ncid, lev_vid, "units", "hPa"))
        NC_CHECK(nf90_put_att(ncid, NF90_GLOBAL, "description", trim(string)))
        NC_CHECK(nf90_put_att(ncid, xmid, "long_name", "x-localization modes"))
        NC_CHECK(nf90_enddef(ncid))
        NC_CHECK(nf90_put_var(ncid, lon_vid, lman%xcoord))
        NC_CHECK(nf90_put_var(ncid, lat_vid, lman%ycoord))
        NC_CHECK(nf90_put_var(ncid, lev_vid, lman%zcoord))
        allocate(px_idx(lman%x3dm(1,1)%m))
        px_idx=(/(i, i=1, lman%x3dm(1,1)%m)/)
        NC_CHECK(nf90_put_var(ncid, px_vid, px_idx))
        NC_CHECK(nf90_put_var(ncid, ybgn_vid, lman%ybgn)) 
        NC_CHECK(nf90_put_var(ncid, yend_vid, lman%yend)) 
        NC_CHECK(nf90_put_var(ncid, zbgn_vid, lman%zbgn)) 
        NC_CHECK(nf90_put_var(ncid, zend_vid, lman%zend)) 
        allocate(xmodes(lman%nx, lman%my, lman%mz, lman%x3dm(1,1)%m))
        do j=1, lman%my
            do k=1, lman%mz
                xmodes(:,j,k,:)=lman%x3dm(j,k)%modes(:,:)
            end do
        end do
        NC_CHECK(nf90_put_var(ncid, xmid, xmodes))
        deallocate(px_idx, xmodes)
        NC_CHECK(nf90_close(ncid))
    end subroutine localman_output_xm  !!}}}

    !! model-dependent
    subroutine localman_input_xm(lman, file_xm)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="localman_input_xm"
        character(len=*) :: file_xm
        type(type_localman) :: lman
        integer :: ncid, lonid, latid, levid, pxid, xmid, &
                   lon_vid, lat_vid, lev_vid, px_vid, i, j, k, &
                   nx, ny, nz, px, myid, mzid, &
                   ybgn_vid, yend_vid, zbgn_vid, zend_vid, my, mz
        real,allocatable :: xmodes(:,:,:,:), lon(:), lat(:), lev(:)
        integer :: int_cycle
        logical :: is_cycle
        if(lman%has_set_xm) then
            RP_ERROR("lman has set xmodes already!")
        end if
        NC_CHECK(nf90_open(trim(file_xm), NF90_NOWRITE, ncid))
        NC_CHECK(nf90_inq_dimid(ncid, "lon", lonid))
        NC_CHECK(nf90_inquire_dimension(ncid, lonid, len=nx)) !!nx
        NC_CHECK(nf90_inq_dimid(ncid, "lat", latid))
        NC_CHECK(nf90_inquire_dimension(ncid, latid, len=ny)) !!ny
        NC_CHECK(nf90_inq_dimid(ncid, "lev", levid))  
        NC_CHECK(nf90_inquire_dimension(ncid, levid, len=nz)) !!nz
        NC_CHECK(nf90_inq_dimid(ncid, "px", pxid))
        NC_CHECK(nf90_inquire_dimension(ncid, pxid, len=px)) !!px
        NC_CHECK(nf90_inq_dimid(ncid, "my", myid))
        NC_CHECK(nf90_inquire_dimension(ncid, myid, len=my)) !!my
        NC_CHECK(nf90_inq_dimid(ncid, "mz", mzid))
        NC_CHECK(nf90_inquire_dimension(ncid, mzid, len=mz)) !!mz
        if(.not. lman%has_setdomain) then
            allocate(lon(nx), lat(ny), lev(nz))
            NC_CHECK(nf90_inq_varid(ncid, "lon", lon_vid))
            NC_CHECK(nf90_get_att(ncid, lon_vid, "is_cycle",  int_cycle))
            if(int_cycle == 1) then
                is_cycle=.true.
            else
                is_cycle=.false.
            end if
            NC_CHECK(nf90_get_var(ncid, lon_vid, lon))
            NC_CHECK(nf90_inq_varid(ncid, "lat", lat_vid))
            NC_CHECK(nf90_get_var(ncid, lat_vid, lat))
            NC_CHECK(nf90_inq_varid(ncid, "lev", lev_vid))
            NC_CHECK(nf90_get_var(ncid, lev_vid, lev))
            call localman_set_domain(lman, lon, lat, lev, is_cycle)
            deallocate(lon,lat,lev)
        else
            ASSUREX(nx == lman%nx)
            ASSUREX(ny == lman%ny)
            ASSUREX(nz == lman%nz)
        end if
        lman%my=my
        lman%mz=mz
        allocate(lman%ybgn(my), lman%yend(my))
        allocate(lman%zbgn(mz), lman%zend(mz))
        NC_CHECK(nf90_inq_varid(ncid, "ybgn", ybgn_vid))
        NC_CHECK(nf90_inq_varid(ncid, "yend", yend_vid))
        NC_CHECK(nf90_inq_varid(ncid, "zbgn", zbgn_vid))
        NC_CHECK(nf90_inq_varid(ncid, "zend", zend_vid))
        NC_CHECK(nf90_get_var(ncid, ybgn_vid, lman%ybgn))
        NC_CHECK(nf90_get_var(ncid, yend_vid, lman%yend))
        NC_CHECK(nf90_get_var(ncid, zbgn_vid, lman%zbgn))
        NC_CHECK(nf90_get_var(ncid, zend_vid, lman%zend))
        allocate(xmodes(nx, my, mz, px))
        NC_CHECK(nf90_inq_varid(ncid, "xmode", xmid))
        NC_CHECK(nf90_get_var(ncid, xmid, xmodes))
        allocate(lman%x3dm(my, mz))
        do j=1, my
            do k=1, mz
                call localmodes_init(lman%x3dm(j,k), nx, px)
                lman%x3dm(j,k)%modes=xmodes(:,j,k,:)
            end do
        end do
        NC_CHECK(nf90_close(ncid))
        deallocate(xmodes)
        lman%has_set_xm=.true.
    end subroutine localman_input_xm  !!}}}

    !! model dependent
    subroutine localman_output_ym(lman, file_ym, string)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="localman_output_ym"
        character(len=*) :: file_ym
        type(type_localman) :: lman
        integer :: ncid, lonid, latid, levid, pyid, ymid, &
                   lon_vid, lat_vid, lev_vid, py_vid, i, j, k, &
                   mzid, zbgn_vid, zend_vid
        integer,allocatable :: py_idx(:)
        real,allocatable :: ymodes(:,:,:) 
        character(len=*) :: string
        integer :: int_cycle
        ASSUREX(lman%has_setdomain)
        ASSUREX(lman%has_set_ym)
        NC_CHECK(nf90_create(trim(file_ym), NF90_CLOBBER, ncid))
        NC_CHECK(nf90_def_dim(ncid, "lon", lman%nx, lonid))
        NC_CHECK(nf90_def_dim(ncid, "lat", lman%ny, latid))
        NC_CHECK(nf90_def_dim(ncid, "lev", lman%nz, levid))
        NC_CHECK(nf90_def_dim(ncid, "py", lman%y2dm(1)%m, pyid))
        NC_CHECK(nf90_def_dim(ncid, "mz", lman%mz, mzid))
        NC_CHECK(nf90_def_var(ncid, "lon",  NF90_FLOAT, (/lonid/), lon_vid))
        NC_CHECK(nf90_def_var(ncid, "lat",  NF90_FLOAT, (/latid/), lat_vid))
        NC_CHECK(nf90_def_var(ncid, "lev",  NF90_FLOAT, (/levid/), lev_vid))
        NC_CHECK(nf90_def_var(ncid, "zbgn", NF90_INT,(/mzid/), zbgn_vid))
        NC_CHECK(nf90_def_var(ncid, "zend", NF90_INT,(/mzid/), zend_vid))
        NC_CHECK(nf90_def_var(ncid, "py",   NF90_FLOAT,(/pyid/), py_vid))
        NC_CHECK(nf90_def_var(ncid, "ymode", NF90_FLOAT, (/latid,mzid,pyid/), ymid))
        NC_CHECK(nf90_put_att(ncid, lon_vid, "units", "degrees_east"))
        if(lman%is_x_cycle) then
            int_cycle=1
        else
            int_cycle=0
        end if
        NC_CHECK(nf90_put_att(ncid, lon_vid, "is_cycle", int_cycle))
        NC_CHECK(nf90_put_att(ncid, lat_vid, "units", "degrees_north"))
        NC_CHECK(nf90_put_att(ncid, lev_vid, "units", "hPa"))
        NC_CHECK(nf90_put_att(ncid, NF90_GLOBAL, "description", trim(string)))
        NC_CHECK(nf90_put_att(ncid, ymid, "long_name", "y-localization modes"))
        NC_CHECK(nf90_enddef(ncid))
        NC_CHECK(nf90_put_var(ncid, lon_vid, lman%xcoord))
        NC_CHECK(nf90_put_var(ncid, lat_vid, lman%ycoord))
        NC_CHECK(nf90_put_var(ncid, lev_vid, lman%zcoord))
        allocate(py_idx(lman%y2dm(1)%m))
        py_idx=(/(i, i=1,lman%y2dm(1)%m)/)
        NC_CHECK(nf90_put_var(ncid, py_vid, py_idx))
        NC_CHECK(nf90_put_var(ncid, zbgn_vid, lman%zbgn)) 
        NC_CHECK(nf90_put_var(ncid, zend_vid, lman%zend)) 
        allocate(ymodes(lman%ny, lman%mz, lman%y2dm(1)%m))
        do k=1, lman%mz
            ymodes(:,k,:)=lman%y2dm(k)%modes(:,:)
        end do
        NC_CHECK(nf90_put_var(ncid, ymid, ymodes))
        deallocate(py_idx, ymodes)
        NC_CHECK(nf90_close(ncid))
    end subroutine localman_output_ym  !!}}}

    !! model-dependent
    subroutine localman_input_ym(lman, file_ym)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="localman_input_ym"
        character(len=*) :: file_ym
        type(type_localman) :: lman
        integer :: ncid, lonid, latid, levid, pyid, ymid, &
                   lon_vid, lat_vid, lev_vid, py_vid, i, j, k, &
                   nx, ny, nz, py, myid, mzid, zbgn_vid, zend_vid, mz
        integer,allocatable :: py_idx(:), zbgn(:), zend(:)
        real,allocatable :: ymodes(:,:,:), lon(:), lat(:), lev(:)
        integer :: int_cycle
        logical :: is_cycle
        if(lman%has_set_ym) then
            RP_ERROR("lman has set ymodes already!")
        end if
        ASSUREX(lman%has_set_xm)
        NC_CHECK(nf90_open(trim(file_ym), NF90_NOWRITE, ncid))
        NC_CHECK(nf90_inq_dimid(ncid, "lon", lonid))
        NC_CHECK(nf90_inquire_dimension(ncid, lonid, len=nx)) !!nx
        NC_CHECK(nf90_inq_dimid(ncid, "lat", latid))
        NC_CHECK(nf90_inquire_dimension(ncid, latid, len=ny)) !!ny
        NC_CHECK(nf90_inq_dimid(ncid, "lev", levid))  
        NC_CHECK(nf90_inquire_dimension(ncid, levid, len=nz)) !!nz
        NC_CHECK(nf90_inq_dimid(ncid, "py", pyid))
        NC_CHECK(nf90_inquire_dimension(ncid, pyid, len=py)) !!py
        NC_CHECK(nf90_inq_dimid(ncid, "mz", mzid))
        NC_CHECK(nf90_inquire_dimension(ncid, mzid, len=mz)) !!mz
        if(.not. lman%has_setdomain) then
            allocate(lon(nx), lat(ny), lev(nz))
            NC_CHECK(nf90_inq_varid(ncid, "lon", lon_vid))
            NC_CHECK(nf90_get_att(ncid, lon_vid, "is_cycle",  int_cycle))
            if(int_cycle == 1) then
                is_cycle=.true.
            else
                is_cycle=.false.
            end if
            NC_CHECK(nf90_get_var(ncid, lon_vid, lon))
            NC_CHECK(nf90_inq_varid(ncid, "lat", lat_vid))
            NC_CHECK(nf90_get_var(ncid, lat_vid, lat))
            NC_CHECK(nf90_inq_varid(ncid, "lev", lev_vid))
            NC_CHECK(nf90_get_var(ncid, lev_vid, lev))
            call localman_set_domain(lman, lon, lat, lev, is_cycle)
            deallocate(lon,lat,lev)
        else
            ASSUREX(nx == lman%nx)
            ASSUREX(ny == lman%ny)
            ASSUREX(nz == lman%nz)
        end if
        ASSUREX(mz == lman%mz)
        allocate(zbgn(mz), zend(mz))
        NC_CHECK(nf90_inq_varid(ncid, "zbgn", zbgn_vid))
        NC_CHECK(nf90_inq_varid(ncid, "zend", zend_vid))
        NC_CHECK(nf90_get_var(ncid, zbgn_vid, zbgn))
        NC_CHECK(nf90_get_var(ncid, zend_vid, zend))
        ASSUREX(all(zbgn == lman%zbgn))
        ASSUREX(all(zend == lman%zend))
        deallocate(zbgn, zend)
        allocate(ymodes(ny, mz, py))
        NC_CHECK(nf90_inq_varid(ncid, "ymode", ymid))
        NC_CHECK(nf90_get_var(ncid, ymid, ymodes))
        allocate(lman%y2dm(lman%mz))
        do k=1, mz
            call localmodes_init(lman%y2dm(k), ny, py)
            lman%y2dm(k)%modes=ymodes(:,k,:)
        end do
        deallocate(ymodes)
        NC_CHECK(nf90_close(ncid))
        lman%has_set_ym=.true.
    end subroutine localman_input_ym  !!}}}

    !! model dependent
    subroutine localman_output_zm(lman, file_zm, string)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="localman_output_zm"
        character(len=*) :: file_zm
        type(type_localman) :: lman
        integer :: ncid, lonid, latid, levid, pzid, zmid, &
                   lon_vid, lat_vid, lev_vid, pz_vid, i, j, k
        integer,allocatable :: pz_idx(:)
        character(len=*) :: string
        integer :: int_cycle
        ASSUREX(lman%has_setdomain)
        ASSUREX(lman%has_set_zm)
        NC_CHECK(nf90_create(trim(file_zm), NF90_CLOBBER, ncid))
        NC_CHECK(nf90_def_dim(ncid, "lon", lman%nx, lonid))
        NC_CHECK(nf90_def_dim(ncid, "lat", lman%ny, latid))
        NC_CHECK(nf90_def_dim(ncid, "lev", lman%nz, levid))
        NC_CHECK(nf90_def_dim(ncid, "pz", lman%zm%m, pzid))
        NC_CHECK(nf90_def_var(ncid, "lon", NF90_FLOAT, (/lonid/), lon_vid))
        NC_CHECK(nf90_def_var(ncid, "lat", NF90_FLOAT, (/latid/), lat_vid))
        NC_CHECK(nf90_def_var(ncid, "lev", NF90_FLOAT, (/levid/), lev_vid))
        NC_CHECK(nf90_def_var(ncid, "pz",  NF90_FLOAT,(/pzid/), pz_vid))
        NC_CHECK(nf90_def_var(ncid, "zmode", NF90_FLOAT, (/levid,pzid/), zmid))
        NC_CHECK(nf90_put_att(ncid, lon_vid, "units", "degrees_east"))
        if(lman%is_x_cycle) then
            int_cycle=1
        else
            int_cycle=0
        end if
        NC_CHECK(nf90_put_att(ncid, lon_vid, "is_cycle", int_cycle))
        NC_CHECK(nf90_put_att(ncid, lat_vid, "units", "degrees_north"))
        NC_CHECK(nf90_put_att(ncid, lev_vid, "units", "hPa"))
        NC_CHECK(nf90_put_att(ncid, NF90_GLOBAL, "description", trim(string)))
        NC_CHECK(nf90_put_att(ncid, zmid, "long_name", "z-localization modes"))
        NC_CHECK(nf90_enddef(ncid))
        NC_CHECK(nf90_put_var(ncid, lon_vid, lman%xcoord))
        NC_CHECK(nf90_put_var(ncid, lat_vid, lman%ycoord))
        NC_CHECK(nf90_put_var(ncid, lev_vid, lman%zcoord))
        allocate(pz_idx(lman%zm%m))
        pz_idx=(/(i, i=1,lman%zm%m)/)
        NC_CHECK(nf90_put_var(ncid, pz_vid, pz_idx))
        NC_CHECK(nf90_put_var(ncid, zmid, lman%zm%modes))
        deallocate(pz_idx)
        NC_CHECK(nf90_close(ncid))
    end subroutine localman_output_zm  !!}}}

    !! model-dependent
    subroutine localman_input_zm(lman, file_zm)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="localman_input_zm"
        character(len=*) :: file_zm
        type(type_localman) :: lman
        integer :: ncid, lonid, latid, levid, pzid, zmid, &
                   lon_vid, lat_vid, lev_vid, pz_vid, i, j, k, &
                   nx, ny, nz, pz
        integer,allocatable :: pz_idx(:)
        real,allocatable :: lon(:), lat(:), lev(:)
        integer :: int_cycle
        logical :: is_cycle
        if(lman%has_set_zm) then 
            RP_ERROR("lman has set zmode already!")
        end if
        NC_CHECK(nf90_open(trim(file_zm), NF90_NOWRITE, ncid))
        NC_CHECK(nf90_inq_dimid(ncid, "lon", lonid))
        NC_CHECK(nf90_inquire_dimension(ncid, lonid, len=nx)) !!nx
        NC_CHECK(nf90_inq_dimid(ncid, "lat", latid))
        NC_CHECK(nf90_inquire_dimension(ncid, latid, len=ny)) !!ny
        NC_CHECK(nf90_inq_dimid(ncid, "lev", levid))  
        NC_CHECK(nf90_inquire_dimension(ncid, levid, len=nz)) !!nz
        NC_CHECK(nf90_inq_dimid(ncid, "pz", pzid))
        NC_CHECK(nf90_inquire_dimension(ncid, pzid, len=pz)) !!pz
        if(.not. lman%has_setdomain) then
            allocate(lon(nx), lat(ny), lev(nz))
            NC_CHECK(nf90_inq_varid(ncid, "lon", lon_vid))
            NC_CHECK(nf90_get_att(ncid, lon_vid, "is_cycle",  int_cycle))
            if(int_cycle == 1) then
                is_cycle=.true.
            else
                is_cycle=.false.
            end if
            NC_CHECK(nf90_get_var(ncid, lon_vid, lon))
            NC_CHECK(nf90_inq_varid(ncid, "lat", lat_vid))
            NC_CHECK(nf90_get_var(ncid, lat_vid, lat))
            NC_CHECK(nf90_inq_varid(ncid, "lev", lev_vid))
            NC_CHECK(nf90_get_var(ncid, lev_vid, lev))
            call localman_set_domain(lman, lon, lat, lev, is_cycle)
            deallocate(lon, lat, lev)
        else
            ASSUREX(nx == lman%nx)
            ASSUREX(ny == lman%ny)
            ASSUREX(nz == lman%nz)
        end if
        NC_CHECK(nf90_inq_varid(ncid, "zmode", zmid))
        call localmodes_init(lman%zm, nz, pz)
        NC_CHECK(nf90_get_var(ncid, zmid, lman%zm%modes))
        NC_CHECK(nf90_close(ncid))
        lman%has_set_zm=.true.
    end subroutine localman_input_zm  !!}}}

    !! normalize for xm/ym/zm
    subroutine localman_rownormalize(lman)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="localman_rownormalize"
        type(type_localman) :: lman
        integer :: j, k
        ASSUREX(lman%has_set_xm)
        ASSUREX(lman%has_set_ym)
        ASSUREX(lman%has_set_zm)
        do j=1, lman%my
            do k=1, lman%mz
                call localmodes_rownormalize(lman%x3dm(j,k))
            end do
        end do
        do k=1, lman%mz
            call localmodes_rownormalize(lman%y2dm(k))
        end do
        call localmodes_rownormalize(lman%zm)
    end subroutine localman_rownormalize  !!}}}

    !! for Sx*vector, calculate matrix multiplication: outf3d(nx,my,mz) = xmode*vector
    subroutine localman_xm_multiply_vector(lman, vector, outf3d)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="localman_xm_multiply_vector"
        type(type_localman) :: lman
        real :: vector(:), outf3d(:,:,:)
        real,allocatable :: x3d(:,:,:) !!nx,my,mz
        integer :: j, k, i, n
        ASSUREX(lman%has_setdomain)
        ASSUREX(lman%has_set_xm)
        ASSUREX(lman%has_set_ym)
        ASSUREX(lman%has_set_zm)
        ASSUREX(size(vector) == lman%x3dm(1,1)%m) !!size(vector) == px
        ASSUREX(all(shape(outf3d)==(/lman%nx, lman%my, lman%mz/)))
        allocate(x3d(lman%nx, lman%my, lman%mz))
        outf3d = 0
        do i=1, size(vector)
            do j=1, lman%my
                do k=1, lman%mz
                    x3d(:,j,k)=lman%x3dm(j,k)%modes(:,i)
                end do
            end do 
            outf3d=outf3d+x3d*vector(i)
        end do
        deallocate(x3d)
    end subroutine localman_xm_multiply_vector  !!}}}

    !! ym(:,ipy) (schur) vector, or diag(ym(:,ipy)) * vector:  f3d(nx,ny,mz)->f3d
    subroutine localman_ymi_schur_f3d(lman, ipy, f3d)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="localman_ymi_schur_f3d"
        type(type_localman) :: lman
        real :: f3d(:,:,:)
        integer :: j, k, ipy
        ASSUREX(lman%has_setdomain)
        ASSUREX(lman%has_set_xm)
        ASSUREX(lman%has_set_ym)
        ASSUREX(lman%has_set_zm)
        ASSUREX(all(shape(f3d)==(/lman%nx, lman%ny, lman%mz/)))
        ASSUREX( ipy <= lman%y2dm(1)%m )
        ASSUREX( ipy >= 1 )
        do k=1, lman%mz
            do j=1, lman%ny
                f3d(:,j,k)=f3d(:,j,k)*lman%y2dm(k)%modes(j,ipy)
            end do
        end do
    end subroutine localman_ymi_schur_f3d  !!}}}

    !! zm(:,ipz) (schur) vector, or diag(zm(:,ipz)) * vector:  f3d(nx,ny,nz)->f3d
    subroutine localman_zmi_schur_f3d(lman, ipz, f3d)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="localman_zmi_schur_f3d"
        type(type_localman) :: lman
        real :: f3d(:,:,:)
        integer :: j, k, ipz
        ASSUREX(lman%has_setdomain)
        ASSUREX(lman%has_set_xm)
        ASSUREX(lman%has_set_ym)
        ASSUREX(lman%has_set_zm)
        ASSUREX(all(shape(f3d)==(/lman%nx, lman%ny, lman%nz/)))
        ASSUREX( ipz <= lman%zm%m )
        ASSUREX( ipz >= 1 )
        do k=1, lman%nz
            f3d(:,:,k)=f3d(:,:,k)*lman%zm%modes(k,ipz)
        end do
    end subroutine localman_zmi_schur_f3d  !!}}}

    !!return neighbors(ny,2), weights(ny,2), which specify the neighbors/weights for the ny grid
    subroutine localman_get_my2ny_Weight(lman, neighbors, weights)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="localman_get_my2ny_Weight"
        type(type_localman) :: lman
        integer :: neighbors(:,:)
        real :: weights(:,:)
        integer :: j,k
        ASSUREX(lman%has_setdomain)
        ASSUREX(lman%has_set_xm)
        if(.not. all((/size(neighbors,1), size(weights,1)/) == lman%ny)) then
            RP_ERROR("input parameter size conflicts")
        end if
        if(.not. all((/size(neighbors,2), size(weights,2)/) == 2)) then
            RP_ERROR("input neighbors/weights 2nd dimension size should equal to 2")
        end if
        do j=1, lman%my
            neighbors(lman%ybgn(j):lman%yend(j),:)=j
            weights(lman%ybgn(j):lman%yend(j),1)=1.
            weights(lman%ybgn(j):lman%yend(j),2)=0.
        end do
        do j=1, lman%my-1
            do k=lman%yend(j)+1, lman%ybgn(j+1)-1
                neighbors(k,1)=j
                neighbors(k,2)=j+1
                weights(k,1)=(lman%ybgn(j+1)-k)*1.0/(lman%ybgn(j+1)-lman%yend(j))
                weights(k,2)=1-weights(k,1)
                ASSUREX(weights(k,1)>0)
                ASSUREX(weights(k,2)>0)
            end do
        end do
    end subroutine localman_get_my2ny_Weight  !!}}}

    subroutine localman_get_mz2nz_Weight(lman, neighbors, weights)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="localman_get_mz2nz_Weight"
        type(type_localman) :: lman
        integer :: neighbors(:,:)
        real :: weights(:,:)
        integer :: j, k
        ASSUREX(lman%has_setdomain)
        ASSUREX(lman%has_set_xm)
        if(.not. all((/size(neighbors,1), size(weights,1)/) == lman%nz)) then
            RP_ERROR("input parameter size conflicts")
        end if
        if(.not. all((/size(neighbors,2), size(weights,2)/) == 2)) then
            RP_ERROR("input neighbors/weights 2nd dimension size should equal to 2")
        end if
        do j=1, lman%mz
            neighbors(lman%zbgn(j):lman%zend(j),:)=j
            weights(lman%zbgn(j):lman%zend(j),1)=1.
            weights(lman%zbgn(j):lman%zend(j),2)=0.
        end do
        do j=1, lman%mz-1
            do k=lman%zend(j)+1, lman%zbgn(j+1)-1
                neighbors(k,1)=j
                neighbors(k,2)=j+1
                weights(k,1)=(lman%zcoord(lman%zbgn(j+1))-lman%zcoord(k)) / &
                             (lman%zcoord(lman%zbgn(j+1))-lman%zcoord(lman%yend(j)))
                weights(k,2)=1-weights(k,1)
                ASSUREX(weights(k,1)>0)
                ASSUREX(weights(k,2)>0)
            end do
        end do
    end subroutine localman_get_mz2nz_Weight  !!}}}

end module module_localman
