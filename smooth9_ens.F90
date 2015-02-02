! --------------------------------------------------------------------------
! Description:
!
!   ./smooth9_ens.exe <ens_list> <out_dir> <do_incr_smooth>
!
! History:
!
!   2014-09-26:
!
!       [shensi]: First creation.
! --------------------------------------------------------------------------
#include "inc_common.fpp"
program smooth9_ens
use utility
use ArgsManager
use module_modelstate
use module_ensemble
use module_types
implicit none
    character(len=*),parameter :: PROCEDURE_NAME="smooth9_ens"
    character(len=200) :: ens_list, out_dir, option, px_list
    character(len=200),allocatable :: ens_files(:), new_files(:), px_files(:)
    logical :: do_incr_smooth !!do smooth on Px-Px_mean
    type(type_modeldesc) :: md
    type(type_ensPx) :: ensPx, ensPx2
    real,allocatable :: xmean(:), xptbd(:)
    integer :: n_ens, i, j, im1, ip1, jm1, jp1, igrid, ivar, offset, k, &
               im1jm1, im1jp1, ip1jp1, ip1jm1
    integer :: nvar=1, nlon=64, nlat=32, nlev=3
    real:: coef_p, coef_q !!see ncl::smth9 definition
    logical :: do_output_ptbd=.true., do_output_dxa=.true.
    integer :: opt

    coef_p=0.5
    coef_q=0.25 !!heavy smooth

    call args_readargs()
    if(args_get_narg() /= 4) then
        RP_ERROR("<usage 1> ./smooth9_ens.exe <list_gx> opt1 <out_dir> <do_incr_smooth--based on ens-mean>")
        RP_ERROR("<usage 2> ./smooth9_ens.exe <list_gx> opt2 <list_px> <out_dir> : smooth base on gx-px")
    end if
    call args_getarg(1, option)
    if(lowercase(trim(option)) == "opt1") then
        call args_getarg(2, ens_list)
        call args_getarg(3, out_dir)
        call args_getarg(4, do_incr_smooth)
        opt=1
    else if(lowercase(trim(option)) == "opt2") then
        call args_getarg(2, ens_list)
        call args_getarg(3, px_list)
        call args_getarg(4, out_dir)
        opt=2
    else
        RP_ERROR("<usage 1> ./smooth9_ens.exe <list_gx> opt1 <out_dir> <do_incr_smooth--based on ens-mean>")
        RP_ERROR("<usage 2> ./smooth9_ens.exe <list_gx> opt2 <list_px> <out_dir> : smooth base on gx-px")
    end if

    call allocread_line(trim(ens_list), n_ens, ens_files)
    if(opt == 2) then
        call allocread_line(trim(px_list), i, px_files)
        if(i/=n_ens) then
            RP_ERROR("i /= n_ens")
        end if
    end if

    call modeldesc_init(md)
    call modeldesc_getcoord(md, trim(ens_files(1)))
    call Px_read(md, ens_files, ensPx)

    allocate(xmean(md%n_total))
    allocate(new_files(n_ens))
    
    if(opt == 1 .and. do_incr_smooth) then
        RP_LOG("opt1, and do_incr_smooth")
        call Px_average(ensPx, xmean)
        call Px_full2pert(ensPx, xmean)
        
        if(do_output_ptbd) then
            call system("test -d perturbed && rm -rf perturbed ; mkdir perturbed")
            do i=1,n_ens
                new_files(i)="perturbed/"//trim(kickout_filepath(trim(ens_files(i))))
                call system("cp "//trim(ens_files(i))//" "//trim(new_files(i)))
            end do
            call Px_write(md, new_files, ensPx)
        end if
    
    !! do smooth
    !! for non-pole grid, smooth by 9 points
    !! for pole, smooth by pole and its nearest lat-circle (0.5: the pole itself, 0.5: others)
        allocate(xptbd(md%n_total))
        do k=1, n_ens
            xptbd=0
            do ivar=1, nvar
                offset=md%vars(ivar)%ibgn-1
                do igrid=1, nlon*nlat
                    j=ceiling(igrid*1.0/nlon) !!j: latidx
                    i=mod(igrid, nlon)  !!i: lonidx
                    if(i == 0) then
                        i=nlon
                    end if
                    if(j==1) then  !!{{{
                        xptbd(offset+igrid)=0.5*sum(ensPx%Px(offset+1:offset+nlon,k))/nlon + &
                                            0.5*sum(ensPx%Px(offset+nlon+1:offset+nlon*2,k))/nlon
                    else if(j==nlat) then
                        xptbd(offset+igrid)=0.5*sum(ensPx%Px(offset+nlon*(nlat-1)+1:offset+nlon*nlat,k))/nlon + &
                                            0.5*sum(ensPx%Px(offset+nlon*(nlat-2)+1:offset+nlon*(nlat-1),k))/nlon
                    else 
                        !!   (i-1,j+1): im1jp1    (i,j+1): jp1   (i+1,j+1): ip1jp1
                        !!   (i-1,j): im1         (i,j): igrid   (i+1,j): ip1
                        !!   (i-1,j-1): im1jm1    (i,j-1): jm1   (i+1,j-1): ip1jm1
                        jm1=igrid-nlon
                        jp1=igrid+nlon
                        im1=i-1
                        ip1=i+1
                        if(im1 == 0 ) im1=nlon
                        if(ip1 == nlon+1) ip1=1
                        im1=im1-i+igrid
                        ip1=ip1-i+igrid
                        im1jp1=im1+nlon
                        im1jm1=im1-nlon
                        ip1jp1=ip1+nlon
                        ip1jm1=ip1-nlon
                        xptbd(offset+igrid)=ensPx%Px(offset+igrid,k) + &
                                    coef_p/4 * (ensPx%Px(offset+im1,k)+ &
                                                ensPx%Px(offset+jp1,k)+ &
                                                ensPx%Px(offset+ip1,k)+ &
                                                ensPx%Px(offset+jm1,k) - &
                                                4*ensPx%Px(offset+igrid,k) ) + &
                                    coef_q/4 * (ensPx%Px(offset+im1jm1,k) + &
                                                ensPx%Px(offset+im1jp1,k) + &
                                                ensPx%Px(offset+ip1jm1,k) + &
                                                ensPx%Px(offset+ip1jp1,k) - &
                                                4*ensPx%Px(offset+igrid,k) )
                    end if  !!}}}
                end do
            end do
            ensPx%Px(:,k)=xptbd
        end do
    
        if(do_output_ptbd) then
            call system("test -d smth_ptbd && rm -rf smth_ptbd ; mkdir smth_ptbd")
            do i=1,n_ens
                new_files(i)="smth_ptbd/"//trim(kickout_filepath(trim(ens_files(i))))
                call system("cp "//trim(ens_files(i))//" "//trim(new_files(i)))
            end do
            call Px_write(md, new_files, ensPx)
        end if
    
    !! pert2full
        call Px_pert2full(ensPx, xmean)
    
    else if(opt == 1) then !! do_full_smooth
        RP_LOG("opt1, and do_full_smooth")
    !! do smooth
    !! for non-pole grid, smooth by 9 points
    !! for pole, smooth by pole and its nearest lat-circle (0.5: the pole itself, 0.5: others)
        allocate(xptbd(md%n_total))
        do k=1, n_ens
            xptbd=0
            do ivar=1, nvar
                offset=md%vars(ivar)%ibgn-1
                do igrid=1, nlon*nlat
                    j=ceiling(igrid*1.0/nlon) !!j: latidx
                    i=mod(igrid, nlon)  !!i: lonidx
                    if(i == 0) then
                        i=nlon
                    end if
                    if(j==1) then  !!{{{
                        xptbd(offset+igrid)=0.5*sum(ensPx%Px(offset+1:offset+nlon,k))/nlon + &
                                            0.5*sum(ensPx%Px(offset+nlon+1:offset+nlon*2,k))/nlon
                    else if(j==nlat) then
                        xptbd(offset+igrid)=0.5*sum(ensPx%Px(offset+nlon*(nlat-1)+1:offset+nlon*nlat,k))/nlon + &
                                            0.5*sum(ensPx%Px(offset+nlon*(nlat-2)+1:offset+nlon*(nlat-1),k))/nlon
                    else 
                        !!   (i-1,j+1): im1jp1    (i,j+1): jp1   (i+1,j+1): ip1jp1
                        !!   (i-1,j): im1         (i,j): igrid   (i+1,j): ip1
                        !!   (i-1,j-1): im1jm1    (i,j-1): jm1   (i+1,j-1): ip1jm1
                        jm1=igrid-nlon
                        jp1=igrid+nlon
                        im1=i-1
                        ip1=i+1
                        if(im1 == 0 ) im1=nlon
                        if(ip1 == nlon+1) ip1=1
                        im1=im1-i+igrid
                        ip1=ip1-i+igrid
                        im1jp1=im1+nlon
                        im1jm1=im1-nlon
                        ip1jp1=ip1+nlon
                        ip1jm1=ip1-nlon
                        xptbd(offset+igrid)=ensPx%Px(offset+igrid,k) + &
                                    coef_p/4 * (ensPx%Px(offset+im1,k)+ &
                                                ensPx%Px(offset+jp1,k)+ &
                                                ensPx%Px(offset+ip1,k)+ &
                                                ensPx%Px(offset+jm1,k) - &
                                                4*ensPx%Px(offset+igrid,k) ) + &
                                    coef_q/4 * (ensPx%Px(offset+im1jm1,k) + &
                                                ensPx%Px(offset+im1jp1,k) + &
                                                ensPx%Px(offset+ip1jm1,k) + &
                                                ensPx%Px(offset+ip1jp1,k) - &
                                                4*ensPx%Px(offset+igrid,k) )
                    end if  !!}}}
                end do
            end do
            ensPx%Px(:,k)=xptbd
        end do
    else  !! opt=2 
        RP_LOG("opt2, and do_smooth on FGX-PX")
        call Px_read(md, px_files, ensPx2)
        do i=1,n_ens
            ensPx%Px(:,i)=ensPx%Px(:,i)-ensPx2%Px(:,i)
        end do
        if(do_output_dxa) then
            call system("test -d dxa_nosm && rm -rf dxa_nosm ; mkdir dxa_nosm")
            do i=1,n_ens
                new_files(i)="dxa_nosm/"//trim(kickout_filepath(trim(ens_files(i))))
                call system("cp "//trim(ens_files(i))//" "//trim(new_files(i)))
            end do
            call Px_write(md, new_files, ensPx)
        end if
    !! do smooth
    !! for non-pole grid, smooth by 9 points
    !! for pole, smooth by pole and its nearest lat-circle (0.5: the pole itself, 0.5: others)
        allocate(xptbd(md%n_total))
        do k=1, n_ens
            xptbd=0
            do ivar=1, nvar
                offset=md%vars(ivar)%ibgn-1
                do igrid=1, nlon*nlat
                    j=ceiling(igrid*1.0/nlon) !!j: latidx
                    i=mod(igrid, nlon)  !!i: lonidx
                    if(i == 0) then
                        i=nlon
                    end if
                    if(j==1) then  !!{{{
                        xptbd(offset+igrid)=0.5*sum(ensPx%Px(offset+1:offset+nlon,k))/nlon + &
                                            0.5*sum(ensPx%Px(offset+nlon+1:offset+nlon*2,k))/nlon
                    else if(j==nlat) then
                        xptbd(offset+igrid)=0.5*sum(ensPx%Px(offset+nlon*(nlat-1)+1:offset+nlon*nlat,k))/nlon + &
                                            0.5*sum(ensPx%Px(offset+nlon*(nlat-2)+1:offset+nlon*(nlat-1),k))/nlon
                    else 
                        !!   (i-1,j+1): im1jp1    (i,j+1): jp1   (i+1,j+1): ip1jp1
                        !!   (i-1,j): im1         (i,j): igrid   (i+1,j): ip1
                        !!   (i-1,j-1): im1jm1    (i,j-1): jm1   (i+1,j-1): ip1jm1
                        jm1=igrid-nlon
                        jp1=igrid+nlon
                        im1=i-1
                        ip1=i+1
                        if(im1 == 0 ) im1=nlon
                        if(ip1 == nlon+1) ip1=1
                        im1=im1-i+igrid
                        ip1=ip1-i+igrid
                        im1jp1=im1+nlon
                        im1jm1=im1-nlon
                        ip1jp1=ip1+nlon
                        ip1jm1=ip1-nlon
                        xptbd(offset+igrid)=ensPx%Px(offset+igrid,k) + &
                                    coef_p/4 * (ensPx%Px(offset+im1,k)+ &
                                                ensPx%Px(offset+jp1,k)+ &
                                                ensPx%Px(offset+ip1,k)+ &
                                                ensPx%Px(offset+jm1,k) - &
                                                4*ensPx%Px(offset+igrid,k) ) + &
                                    coef_q/4 * (ensPx%Px(offset+im1jm1,k) + &
                                                ensPx%Px(offset+im1jp1,k) + &
                                                ensPx%Px(offset+ip1jm1,k) + &
                                                ensPx%Px(offset+ip1jp1,k) - &
                                                4*ensPx%Px(offset+igrid,k) )
                    end if  !!}}}
                end do
            end do
            ensPx%Px(:,k)=xptbd
        end do
        if(do_output_dxa) then
            call system("test -d dxa_smth && rm -rf dxa_smth ; mkdir dxa_smth")
            do i=1,n_ens
                new_files(i)="dxa_smth/"//trim(kickout_filepath(trim(ens_files(i))))
                call system("cp "//trim(ens_files(i))//" "//trim(new_files(i)))
            end do
            call Px_write(md, new_files, ensPx)
        end if
        do k=1,n_ens
            ensPx%Px(:,k)=ensPx%Px(:,k)+ensPx2%Px(:,k)
        end do
    end if    


    do i=1, n_ens
        new_files(i)=trim(out_dir)//"/"//trim(kickout_filepath(trim(ens_files(i))))
        call system("cp "//trim(ens_files(i))//" "//trim(new_files(i)))
    end do
    call Px_write(md, new_files, ensPx)
end program smooth9_ens
