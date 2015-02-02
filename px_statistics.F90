!! add local2dman
#include "inc_common.fpp"
program px_statistics
use utility
use module_types
use module_ensemble
use module_modelstate
use module_config
implicit none
    character(len=*),parameter :: PROCEDURE_NAME="px_statistics"
    type(type_modeldesc) :: md
    type(type_localdesc) :: ld
    type(type_ensPx) :: ensPx
  
    logical :: do_stdv, do_cova, do_corr, do_horizon, do_vertical, use_localman
    logical :: do_zone_average, do_zone_abs_average
    
    character(len=100) :: list_files, prefix
    integer :: npts, n_mem, i, idx, j, s, k
    integer,parameter :: max_npt=100
    integer :: lon_idx(max_npt), lat_idx(max_npt), lev_idx(max_npt)
    character(len=10) :: pts_name(max_npt)

    integer,parameter :: nvar=1, nx=64, ny=32, nz=3
    character(len=3) :: var_names(nvar)=(/"psi"/)
    integer :: var_idx(nvar)=(/1/)
    character(len=100),allocatable :: files(:)
    character(len=8),parameter :: stdvfile="stdv.nc", &
                              cova_prefix="cova_", corr_prefix="corr_"
    character(len=50) :: cova_file, corr_file
    
    real,allocatable :: stdv(:), corr(:), cova(:), row(:), column(:), xmean(:), &
                        corr_item(:), cova_item(:)
    real :: x3d(nx,ny,nz,nvar)

    namelist /namelist_px/ do_stdv, do_cova, do_corr, list_files, prefix, &
                             npts, lon_idx, lat_idx, lev_idx, pts_name, &
                             do_horizon, do_vertical, use_localman, &
                             do_zone_average, do_zone_abs_average

    lon_idx=-1
    lat_idx=-1
    pts_name=" "
    open(11, file='namelist.px')
    read(11, nml=namelist_px)
    close(11)
    if(npts>=0 .and. count(lon_idx > 0) /= npts) then
        RP_ERROR("lon_idx & npts are not compatible")
    end if
    if(npts>=0 .and. count(lat_idx > 0) /= npts) then
        RP_ERROR("lat_idx & npts are not compatible")
    end if
    if(npts>=0 .and. count(lev_idx > 0) /= npts) then
        RP_ERROR("lev_idx & npts are not compatible")
    end if
    if(npts>=0 .and. count(len_trim(pts_name) > 0) /= npts) then
        RP_ERROR("pts_name & npts are not compatible")
    end if

    if(do_zone_abs_average .and. (.not. do_zone_average)) then
        RP_ERROR("if you want do_zone_abs_average, you should set do_zone_average .true. as well")
    end if

    if(npts == -1) then !! all lat/lev, lon=1
        npts=ny*nz
        do i=1, nz
            lat_idx((i-1)*ny+1: i*ny)=(/(j,j=1,ny)/)
            lev_idx((i-1)*ny+1: i*ny)=i
        end do
        lon_idx(1:npts)=1
        do i=1, nz
            do j=1, ny
                idx=(i-1)*ny+j
                write(pts_name(idx),"(A, I2.2, A, I1.1)") "y",j,"z",i
            end do
        end do
    end if

    call allocread_line(trim(list_files), n_mem, files)
    call modeldesc_init(md)
    call modeldesc_getcoord(md, files(1))
    call Px_read(md, files, ensPx)
    allocate(xmean(md%n_total))
    call Px_average(ensPx, xmean)
    call Px_full2pert(ensPx, xmean)

    ld%ishorizon=do_horizon
    ld%isvertical=do_vertical
    ld%use_localman=use_localman
    ld%lman_xm_file="global_xm.nc"
    ld%lman_ym_file="global_ym.nc"
    ld%lman_zm_file="global_zm.nc"

    call Px_localize(ensPx, ld)
    allocate(stdv(md%n_total), cova(md%n_total), corr(md%n_total))
    allocate(cova_item(md%n_total), corr_item(md%n_total))
    call Px_get_stdv_px(ensPx, stdv) 
    stdv=stdv/sqrt(1.0*(ensPx%n_mem-1))

    TRACEX(ensPx%n_mem)
    if(do_stdv) then
        call system("cp "//trim(files(1))//" "//trim(prefix)//trim(stdvfile))
        call modelstate_output(md, trim(prefix)//trim(stdvfile), stdv)
        RP_LOG("stdv file output --")
    end if

    allocate(row(ensPx%n_total_mem))
    allocate(column(ensPx%n_state))
    do j=1, npts
        RP_LOG("-process point:"//trim(pts_name(j)))
        if(do_zone_average) then
            cova=0
            corr=0
            do s=1, nx
                idx=modelstate_get_idx(md, 1, s, lat_idx(j), lev_idx(j))
                call Px_get_row(md, ensPx, 1, s, lat_idx(j), lev_idx(j), row)
                call Px_multiply_vector(md, ensPx, row, column)
                cova_item=column/((ensPx%n_mem-1))
                x3d=reshape(cova_item, (/nx,ny,nz,nvar/))
                do k=1, nvar 
                    call xshift_cova_corr(nx, ny, nz, x3d(:,:,:,k), 1-s)
                end do
                cova_item=reshape(x3d, (/nx*ny*nz*nvar/))
                if(do_zone_abs_average) then
                    cova=cova+abs(cova_item)
                else
                    cova=cova+cova_item
                end if
                if(do_corr) then
                    corr_item=column/(stdv*stdv(idx)*(ensPx%n_mem-1))
                    x3d=reshape(corr_item,(/nx,ny,nz,nvar/))
                    do k=1, nvar
                        call xshift_cova_corr(nx, ny, nz, x3d(:,:,:,k), 1-s)
                    end do
                    TRACEX(x3d(1:3,lat_idx(j),lev_idx(j),1))
                    corr_item=reshape(x3d, (/nx*ny*nz*nvar/))
                    if(do_zone_abs_average) then
                        corr=corr+abs(corr_item)
                    else
                        corr=corr+corr_item
                    end if
                end if
            end do
            cova=cova/nx
            corr=corr/nx
        else 
            idx=modelstate_get_idx(md, 1, lon_idx(j), lat_idx(j), lev_idx(j))
            RP_LOG("--get row")
            call Px_get_row(md, ensPx, 1, lon_idx(j), lat_idx(j), lev_idx(j), row)
            RP_LOG("--matmul")
            call Px_multiply_vector(md, ensPx, row, column)
            cova=column/(ensPx%n_mem-1)
            corr=cova/(stdv*stdv(idx))
        end if

        if(do_cova) then
            write(cova_file, *) trim(prefix)//trim(cova_prefix)//"psi_"//trim(pts_name(j))//".nc"
            cova_file=adjustl(cova_file)
            call system("cp "//trim(files(1))//" "//trim(cova_file))
            call modelstate_output(md, trim(cova_file), cova)
            RP_LOG("--cova file output --")
        end if
        if(do_corr) then
            write(corr_file, *) trim(prefix)//trim(corr_prefix)//"psi_"//trim(pts_name(j))//".nc"
            corr_file=adjustl(corr_file)
            call system("cp "//trim(files(1))//" "//trim(corr_file))
            call modelstate_output(md, trim(corr_file), corr)
            RP_LOG("--corr file output --")
        end if
    end do
end program    


!! 
subroutine xshift_cova_corr(nx, ny, nz, field, m)  !!{{{
implicit none
    character(len=*),parameter :: PROCEDURE_NAME="xshift_cova_corr"
    integer:: nx, ny, nz
    real :: field(nx,ny,nz), outf(nx,ny,nz), xidx(nx)
    integer :: m, i, j, k
    xidx=(/(i-m,i=1,nx)/)
    where(xidx <=0)
        xidx=xidx+nx
    end where
    where(xidx > nx)
        xidx=xidx-nx
    end where
    outf=field(xidx,:,:)
    field=outf
end subroutine xshift_cova_corr  !!}}}
