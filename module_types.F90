module module_types
implicit none

    ! variables
    !! on model-descriptor
    type :: type_modelvardesc
        integer :: n_grids, ibgn, iend
        character(len=8) :: varname
        integer :: xcoord_id, ycoord_id, zcoord_id
    end type

    type :: type_coord
        integer :: n
        real,allocatable :: val(:)
        logical :: iscycle
    end type

    type :: type_modeldesc
        integer :: n_total
        integer :: n_xcoord, n_ycoord, n_zcoord
        type(type_coord),allocatable :: xcoords(:), ycoords(:), zcoords(:)
        integer :: nvar
        type(type_modelvardesc),allocatable :: vars(:)
        logical :: has_init
    end type

    !! on obs-descriptor
    type :: type_obsgroup
        character(len=8) :: group_name
        integer :: n
        integer,allocatable :: idx(:)
    end type

    type :: type_obsdesc
        logical :: has_read
        integer :: n_obs !!total
        real,allocatable :: xpos(:), ypos(:), zpos(:), error(:)
        integer,allocatable :: qcflag(:)
        integer,allocatable :: tidx(:)
        character(len=8),allocatable :: varname(:)
        !! groups: if want no group, has_group=.false.
        logical :: has_group=.true.
        logical :: has_group_set=.false.
        integer :: n_group=3
        type(type_obsgroup) :: group(3)
    end type  

    integer,parameter :: OG_L2=1, OG_L5=2, OG_L8=3 !! od%og(:)

    !!on local
    type :: type_localmodes
        integer :: lenn=0, m=0
        real,allocatable :: modes(:,:) !!lenn, m
    end type

    type :: type_localdesc
        logical :: ishorizon, isvertical
        character(len=10) :: modes_method !!"decomp1d","decomp2d","read","sin"..
        real :: rh0, rv0 !!when modes is generated(decomp1d, decomp2d , sin..
        !! localman: if use localman, xm/ym/zm files are generated before
        logical :: use_localman !!
        character(len=100) :: lman_xm_file, lman_ym_file, lman_zm_file
    end type


    type :: type_localman
        !!localization domain info
        integer :: nx, ny, nz
        real,allocatable :: xcoord(:), ycoord(:), zcoord(:)
        logical :: is_x_cycle=.true.
        !!a even coarser grid
        integer :: my, mz 
        integer,allocatable :: ybgn(:), yend(:) !!my (1..my)->(1..ny)
        integer,allocatable :: zbgn(:), zend(:) !!mz
        !!localization localmodes info 
        type(type_localmodes),allocatable :: x3dm(:,:) !! my*mz
        type(type_localmodes),allocatable :: y2dm(:) !! mz
        type(type_localmodes) :: zm !!
        logical :: has_setdomain = .false., &
                   has_set_xm = .false., &
                   has_set_ym = .false., &
                   has_set_zm = .false., &
                   is_x_loc=.false., &
                   is_y_loc=.false., &
                   is_z_loc=.false. !! xm can be set nolocal (has_set_xm=T, is_x_loc=F)
    end type

    !! on ensemble
    type :: type_ensPx
        integer :: n_mem, n_state, n_total_mem
        real,allocatable :: Px(:,:) !!n_state * n_mem
        !!on localization
!        type(type_localmodes),allocatable :: xlm(:), ylm(:), zlm(:)
!        type(type_local2dman) :: l2dm
        type(type_localman) :: lman       
        logical :: has_localize=.false. !!including
!        logical :: use_l2dm=.false.
    end type

    type :: type_ensPy
        integer :: n_mem, n_state, n_total_mem
        real,allocatable :: Py(:,:) !!n_state * n_mem
        !!on localization
        type(type_localmodes) :: xlm, ylm, zlm
        logical :: has_localize=.false.
    end type

    type :: type_beta
        integer :: m1, m2, m3, m4
        real, allocatable :: array(:)
    end type

    !! some can be read from namelist while some is set in program
    type :: type_drpconfig
       ! real :: len_window  
        integer :: n_obs_slot !!slots are 0,1..,6
        logical :: with_Bterm 
        logical :: with_Cterm !!constrain : Jc
        integer :: Balpha_option !! 1/2 (for EnKF run, 2 inflated as paper 2000) 3 (for ETKF run)
        logical :: do_Py_EOF !! for ensemble assi./fcst. , false (because ens.mem will change )?
        integer :: pertPy_option !! for ETKF(run): 1: perturb Py around yb(read from ybfile), ybase=yb; 2: decenter Py: ybase=mean(Py)
        logical :: do_read_yb !! true: read yb (from ybfile), used by calculate d(yo-yb); false, use d=yo-mean(Py), i.e. yb=mean(Py)
        integer :: inflate_option !! 0: no inflate; 1: infate by inflator; 2. adaptively inflate; 
        !! && 3. adpatively inflate (use a smooth: 0.2*a(-2) + 0.3*a(-1) + 0.5*a(0)) : all these only in ETKF-run
        real :: inflator !! for inflate_option=1
        real :: low_limit_inflator !!for inflate_option=2
        real :: up_limit_inflator !! for inflate_option=2, set < 0 means no up limit
        logical :: do_ETKF !!if to generate analysis ensemble, set true
        integer :: max_nt !! max of max iterations
        real :: eps !!iterations stop when ||grad J|| reaches eps of inital ||grad J||
        integer :: minisolver !!=1, CG ; =0, direct solve
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
        logical :: to_output_Bbalpha !! to output Balpha * balpha
        logical :: to_output_jval  
        logical :: to_output_jgrad
        logical :: to_output_ob
        logical :: to_output_inflator
        logical :: to_output_pydcorr
    end type

    !!global vars
    real :: timer_bgn, timer_last=-1, timer

end module module_types
