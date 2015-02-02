#include "inc_common.fpp"
module module_config
use utility
use module_types
implicit none

contains
    
    !!get dc,ld by reading namelist.drp
    subroutine config_get_namelist(dc, ld)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="config_get_namelist"
        type(type_drpconfig) :: dc
        type(type_localdesc) :: ld
        
        integer :: n_obs_slot !!slots are 0,1..,6
        logical :: with_Bterm 
        logical :: with_Cterm !!constrain : Jc
        integer :: Balpha_option !! 1/2 (for EnKF run, 2 inflated as paper 2000) 3 (for ETKF run)
        logical :: do_Py_EOF !! for ensemble assi./fcst. , false (because ens.mem will change )?
        integer :: pertPy_option !! for ETKF(run): 1: perturb Py around yb(read from ybfile); 2: decenter Py
        logical :: do_read_yb
        integer :: inflate_option !! 0: no inflate; 1: infate by inflator; 2. adaptively inflate
        real :: inflator
        real :: low_limit_inflator, up_limit_inflator
        logical :: do_ETKF !!if to generate analysis ensemble, set true
        integer :: max_nt !! max of max iterations
        real :: eps !!iterations stop when ||grad J|| reaches eps of inital ||grad J||
        integer :: minisolver !!=0, direct solve; =1, CG
        integer :: analysis_time !!equals to one obs_slot
        logical :: is_outerloop !!is part of a outerloop 
        logical :: is_outerloop_1st !! is the first iteration of a outerloop
        integer :: n_mem
        integer :: max_n_obs
        !!input files
        character(len=100), allocatable :: pxfiles(:), pyfiles(:,:), px00files(:)
        character(len=100) :: xbfile, xgfile, xb00file, xg00file
        character(len=100) :: betagfile  !! for I/O
        character(len=100), allocatable :: obfile(:), ybfile(:)
        !!output files       
        character(len=100) :: xafile, xa00file
        character(len=100), allocatable :: gxfiles(:)
        !!diagnosis output setting
        logical :: to_output_Bbalpha=.true. !! to output Balpha & balpha
        logical :: to_output_jval=.true.
        logical :: to_output_jgrad=.true.
        logical :: to_output_ob=.true.
        logical :: to_output_inflator=.true.
        logical :: to_output_pydcorr=.true.
        !!localdesc
        logical :: ishorizon, isvertical, use_localman
!        character(len=10) :: modes_method  !!only support "readinterp"
!        character(len=10) :: l2dm_method  !!only support "readinterp"
        !!local variable
        integer :: i,j
        character(len=4) :: suffix
        character(len=2) :: suffix2


        namelist /drpsettings/ n_obs_slot ,with_Bterm ,with_Cterm ,Balpha_option ,do_Py_EOF , &
                               pertPy_option , do_read_yb, inflate_option ,inflator ,do_ETKF ,max_nt ,eps ,&
                               minisolver ,analysis_time ,is_outerloop ,is_outerloop_1st ,n_mem ,max_n_obs ,&
                               low_limit_inflator, up_limit_inflator
        namelist /runsettings/ to_output_Bbalpha, to_output_jval, to_output_jgrad, to_output_ob, &
                               to_output_inflator, to_output_pydcorr
        namelist /localization/ ishorizon, isvertical, use_localman

        RP_LOG("DRP read namelist.drp")
        open(10, file="namelist.drp")
        read(10, nml=drpsettings)
        read(10, nml=runsettings)
        read(10, nml=localization)
        close(10)
        
        !!get ld
        ld%ishorizon=ishorizon
        ld%isvertical=isvertical
        ld%use_localman=use_localman
        !ld%modes_method=modes_method
        !if(ld%modes_method == "readinterp" ) then
        !    ld%rh0=0. !!degree
        !    ld%rv0=0
        !else
        !    RP_ERROR("modes_method="//trim(modes_method)//" Not recognized or supported")
        !end if
        ld%lman_xm_file="global_xm.nc"
        ld%lman_ym_file="global_ym.nc"
        ld%lman_zm_file="global_zm.nc"

        !!get dc
        dc%n_obs_slot=n_obs_slot 
        dc%with_Bterm=with_Bterm
        dc%with_Cterm=with_Cterm
        dc%Balpha_option=Balpha_option
        dc%do_Py_EOF=do_Py_EOF
        dc%pertPy_option=pertPy_option
        dc%do_read_yb=do_read_yb
        dc%inflate_option=inflate_option
        dc%inflator=inflator
        dc%low_limit_inflator=low_limit_inflator
        dc%up_limit_inflator=up_limit_inflator
        dc%do_ETKF=do_ETKF
        dc%max_nt=max_nt
        dc%eps=eps
        dc%minisolver=minisolver
        dc%analysis_time=analysis_time
        dc%n_mem=n_mem
        dc%is_outerloop=is_outerloop
        dc%is_outerloop_1st=is_outerloop_1st
        dc%max_n_obs=max_n_obs
        allocate(dc%pxfiles(dc%n_mem))
        allocate(dc%px00files(dc%n_mem))
        allocate(dc%pyfiles(dc%n_obs_slot, dc%n_mem))
        allocate(dc%obfile(dc%n_obs_slot))
        allocate(dc%ybfile(dc%n_obs_slot))
        allocate(dc%gxfiles(dc%n_mem))
        dc%xbfile="xb"
        dc%xgfile="xg"
        dc%xb00file="xb00"
        dc%xg00file="xg00"
        dc%betagfile="betag"
        dc%xafile="xa"
        dc%xa00file="xa00"
        do i=1, dc%n_mem
            write(suffix, '(A1,I3.3)') 'e',i
            dc%pxfiles(i)="ENS/px."//trim(suffix)
            dc%px00files(i)="ENS/px00."//trim(suffix)
            do j=1, dc%n_obs_slot
                write(suffix2, '(I2.2)') j-1
                dc%pyfiles(j,i)="ENS/py"//trim(suffix2)//"."//trim(suffix)
            end do
            dc%gxfiles(i)="ENS/gx."//trim(suffix)
        end do
        do j=1, dc%n_obs_slot
            write(suffix2, '(I2.2)') j-1
            dc%ybfile(j)="yb"//trim(suffix2)
            dc%obfile(j)="ob"//trim(suffix2)
        end do
        !!diagnosis output settings
        dc%to_output_Bbalpha=to_output_Bbalpha
        dc%to_output_jval=to_output_jval
        dc%to_output_jgrad=to_output_jgrad
        dc%to_output_ob=to_output_ob
        dc%to_output_inflator=to_output_inflator
        dc%to_output_pydcorr=to_output_pydcorr
    end subroutine config_get_namelist  !!}}}


    subroutine config_check_settings(dc, ld)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="config_check_settings"
        type(type_drpconfig) :: dc
        type(type_localdesc) :: ld
        logical :: lexist
        logical :: to_io_00 !! to read/write **00 file
        logical :: to_read_yb !! if use the averge of Py, no need to read yb
        logical :: to_read_beta
        logical,allocatable :: lexist_d(:), lexist_yb(:), lexist_py(:,:)
        integer :: i,j
        RP_LOG("DRP check namelist settings")
        !!check ld
        if(.not. ld%use_localman) then
            RP_ERROR("use_localman must be set .true.")
        end if
        if(ld%ishorizon) then
            inquire(file=trim(ld%lman_xm_file), exist=lexist)
            if(.not. lexist) then
                RP_ERROR("lman_xm_file NOT exist!")
            end if
            inquire(file=trim(ld%lman_ym_file), exist=lexist)
            if(.not. lexist) then
                RP_ERROR("lman_ym_file NOT exist!")
            end if
        end if
        if(ld%isvertical) then
            inquire(file=trim(ld%lman_zm_file), exist=lexist)
            if(.not. lexist) then
                RP_ERROR("lman_zm_file NOT exist!")
            end if
        end if
        !!check dc
        if(dc%Balpha_option < 1 .or. dc%Balpha_option > 3) then
            RP_ERROR("Balpha_option not recognized")
        end if
        if(dc%pertPy_option < 1 .or. dc%pertPy_option > 2) then
            RP_ERROR("pertPy_option not recognized")
        end if
        if(dc%inflate_option < 0 .or. dc%inflate_option > 3) then
            RP_ERROR("inflate_option not recognized")
        end if
        if(dc%inflate_option == 2 .and. dc%low_limit_inflator <= 0) then
            RP_WARNING("low_limit_inflator is too small")
        end if
        if(dc%inflate_option == 2 .and. dc%up_limit_inflator <=0) then
            RP_WARNING("up_limit_inflator <= 0, so no up limit for adaptive inflation")
        end if
        if(dc%inflate_option == 2 .and. dc%up_limit_inflator > 0 .and. &
           dc%low_limit_inflator>0 .and. (dc%low_limit_inflator >= dc%up_limit_inflator)) then
            RP_ERROR("up_limit_inflator <= low_limit_inflator")
        end if
        if(dc%minisolver /= 0 .and. dc%minisolver /= 1) then
            RP_ERROR("minisolver not recognized")
        end if
        if(dc%minisolver == 0 .and. ld%ishorizon) then
            RP_ERROR("you select minisoler==0(solve directly) but use localization, that's not supported")
        end if
        if(dc%analysis_time < 0 .or. dc%analysis_time > dc%n_obs_slot-1) then
            RP_ERROR("analysis_time not in [0, n_obs_slot-1]")
        end if
        if(dc%is_outerloop_1st .and. (.not. dc%is_outerloop)) then
            RP_ERROR("is_outerloop & is_outerloop_1st conflicts!")
        end if
        if(dc%n_mem < 5) then !!minimal member is 5
            RP_ERROR("n_mem < 5, too small!")
        end if
        if(dc%max_n_obs < 1) then
            RP_ERROR("max_n_obs < 1, too small!")
        end if
        !check settings consistence
        if(dc%Balpha_option /=3 .and. dc%do_ETKF) then
            RP_WARNING("Balpha_option /= 3 but do_ETKF=.true., not reasonable")
        end if
        if(dc%Balpha_option /=3 .and. dc%pertPy_option == 2) then
            RP_WARNING("Balpha_option /= 3(seems you use DRP-EnKF), pertPy_option = 2, not reasonable")
        end if
    !    if(dc%Balpha_option /=3 .and. dc%inflate_option /=0) then
    !        RP_WARNING("Balpha_option /= 3(seems you use DRP-EnKF), inflate option /=0, so inflation is IGNORED")
    !        dc%inflate_option=0
    !    end if
        if(dc%inflate_option /= 1 .and. dc%inflator /= 1.0) then
            RP_WARNING("inflate_option /= 1 and inflator /= 1.0, the inflator value is INGORED")
        end if
        if(dc%inflate_option == 1 .and. dc%inflator < 1.0) then
            RP_WARNING("the inflator < 1.0, might be a mistake!")
        end if
        if(.not. dc%do_read_yb .and. dc%pertPy_option /= 2) then
            RP_ERROR("do_read_yb = .false. but pertPy_option /= 2")
        end if
        !!check input file existence
        !xg/xb/Px
        inquire(file=trim(dc%xbfile), exist=lexist)
        if(.not. lexist) then
            RP_ERROR("xbfile="//trim(dc%xbfile)//" not exist!")
        end if
        inquire(file=trim(dc%xgfile), exist=lexist)
        if(.not. lexist) then
            RP_ERROR("xgfile="//trim(dc%xgfile)//" not exist!")
        end if
        do i=1, dc%n_mem
            inquire(file=trim(dc%pxfiles(i)), exist=lexist)
            if(.not. lexist) then
               RP_ERROR("pxfiles(i)="//trim(dc%pxfiles(i))//" not exist!")
            end if
        end do
        !xg00/xb00/px00
        to_io_00=(dc%analysis_time /= 0 .and. dc%is_outerloop)
        if(to_io_00) then !!xb00, xg00, px00
            inquire(file=trim(dc%xb00file), exist=lexist)
            if(.not. lexist) then
                RP_ERROR("xb00file="//trim(dc%xb00file)//" not exist!")
            end if
            inquire(file=trim(dc%xg00file), exist=lexist)
            if(.not. lexist) then
                RP_ERROR("xg00file="//trim(dc%xg00file)//" not exist!")
            end if
            do i=1, dc%n_mem
                inquire(file=trim(dc%px00files(i)), exist=lexist)
                if(.not. lexist) then
                   RP_ERROR("px00files(i)="//trim(dc%px00files(i))//" not exist!")
                end if
            end do
        end if
        !beta
        to_read_beta=(dc%is_outerloop .and. (.not. dc%is_outerloop_1st))
        if(to_read_beta) then
            inquire(file=trim(dc%betagfile), exist=lexist)
            if(.not. lexist) then
                RP_ERROR("betagfile="//trim(dc%betagfile)//" not exist!")
            end if
        end if
        !on ys (d, yb(on some cases), py), their every time slots file existences must be consistent
        allocate(lexist_d(dc%n_obs_slot), lexist_yb(dc%n_obs_slot), lexist_py(dc%n_obs_slot, dc%n_mem))
        to_read_yb=dc%do_read_yb
        do i=1, dc%n_obs_slot
            inquire(file=trim(dc%obfile(i)), exist=lexist_d(i))
            inquire(file=trim(dc%ybfile(i)), exist=lexist_yb(i))
            do j=1, dc%n_mem
                inquire(file=trim(dc%pyfiles(i,j)), exist=lexist_py(i,j))
            end do
            if(.not. all(lexist_py(i,:) == lexist_d(i))) then
                RP_ERROR("the existences of pyfiles & obfile are not consistent")
            end if
            if(to_read_yb .and. lexist_d(i) /= lexist_yb(i)) then
                RP_ERROR("the existences of obfile & ybfile are not consistent")
            end if
        end do
        !input_inflators.csv
        if(dc%inflate_option == 3) then
            inquire(file=trim("input_inflators.csv"), exist=lexist)
            if(.not. lexist) then
                RP_ERROR("use inflate_option==3 but input_inflators.csv NOT exist!")
            end if
        end if
    end subroutine config_check_settings  !!}}}

    !!depressed: ld is set in code
    subroutine localdesc_setup(ld)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="localdesc_setup"
        type(type_localdesc) :: ld
        ld%ishorizon=.true.
        ld%isvertical=.false.
!       ld%modes_method="readinterp"
        ld%rh0=10. !!degree
        ld%rv0=1
 !      ld%xmode_file="global_xmodes.csv"
 !      ld%ymode_file="global_ymodes.csv"
 !      ld%zmode_file="UNTITLED"
        ld%lman_xm_file="global_xm.nc"
        ld%lman_ym_file="global_ym.nc"
        ld%lman_zm_file="global_zm.nc"
    end subroutine localdesc_setup  !!}}}


end module


