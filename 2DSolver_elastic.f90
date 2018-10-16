!
!---------------------------------!
! Wavefield extrapolation example !
!---------------------------------!
!
! Written by Yder MASSON, september 20, 2018
! email: yder_masson@berkeley.edu
!
!-----------
program main
!-----------

  implicit none

  ! Grid parameters
  integer :: nx            ! grid dimension x
  integer :: ny            ! grid dimension y
  double precision :: dlx   ! grid spacing x
  double precision :: dly   ! grid spacing y

  ! physical properties of the propagating medium
  double precision, allocatable :: vp(:,:)     ! p-wave velocity
  double precision, allocatable :: vs(:,:)     ! s-wave velocity
  double precision, allocatable :: lambda(:,:) ! lame parameter
  double precision, allocatable :: mu(:,:)     ! shear modulus
  double precision, allocatable :: rho(:,:)    ! density

  ! source parameters
  integer :: n_src ! number of sources

  integer, allocatable :: i_src(:) ! source position x
  integer, allocatable :: j_src(:) ! source position y

  integer :: point_per_wavelength  ! minimum number of point per wavelength 
  double precision :: frequency    ! central frequency of source time function

  double precision, allocatable :: stf_src_x(:,:) ! Source time functions x (for all sources)
  double precision, allocatable :: stf_src_y(:,:) ! Source time functions y (for all sources)

  ! receiver parameters
  integer :: n_rec ! number of receivers

  integer, allocatable :: i_rec(:) ! receiver position x
  integer, allocatable :: j_rec(:) ! receiver position y

  double precision, allocatable :: ux_rec(:,:) ! for storing displacement seismograms x
  double precision, allocatable :: uy_rec(:,:) ! for storing displacement seismograms y

  ! time parameters
  integer :: nt          ! number of time steps
  double precision :: dt ! time sampling interval

  ! other parameters
  integer :: is,ir,it ! looping indices

  ! Mirror parameters
  double precision, allocatable :: win_mir(:,:) ! mirror window function

  integer :: i1m ! bottom left corner of mirror window function index i
  integer :: i2m ! top right   corner of mirror window function index i
  integer :: j1m ! bottom left corner of mirror window function index j
  integer :: j2m ! top right   corner of mirror window function index j

  integer, allocatable :: i_mir(:) ! position indices i of mirror points (where Fm/=0)
  integer, allocatable :: j_mir(:) ! position indices j of mirror points (where Fm/=0)
  double precision, allocatable :: fx_mir(:,:) ! Mirror force Fxm, for all mirror points for all time steps
  double precision, allocatable :: fy_mir(:,:) ! Mirror force Fym, for all mirror points for all time steps
  
  ! For single point test
  double precision, allocatable :: tmp_x(:,:), tmp_y(:,:)
  double precision, allocatable :: ux_single_point(:)
  double precision, allocatable :: uy_single_point(:)
  
  double precision, allocatable :: green_xx(:)
  double precision, allocatable :: green_xy(:)
  double precision, allocatable :: green_yx(:)
  double precision, allocatable :: green_yy(:)
  
  ! Green functions parameters
  integer, allocatable :: i_green(:) ! position indices i of mirror points (where Fm/=0)
  integer, allocatable :: j_green(:) ! position indices j of mirror points (where Fm/=0)
  double precision, allocatable :: fxx_green(:,:) ! To store green's functions (response ux for an excitation fx)
  double precision, allocatable :: fxy_green(:,:) ! To store green's functions (response uy for an excitation fx)
  double precision, allocatable :: fyx_green(:,:) ! To store green's functions (response ux for an excitation fy)
  double precision, allocatable :: fyy_green(:,:) ! To store green's functions (response uy for an excitation fy)
  double precision, allocatable :: ux_rec_extrapolated(:) ! extrapolated displacement seismograms x
  double precision, allocatable :: uy_rec_extrapolated(:) ! extrapolated diisplacment seismograms y

  ! array for storing the delta function (for computing the Green's functions)
  double precision, allocatable :: delta_x(:,:) ! it is 2D for consistency with the solverbut really a 1D array
  double precision, allocatable :: delta_y(:,:) ! it is 2D for consistency with the solverbut really a 1D array
  double precision :: f0
  integer :: n_mir
  integer :: ic
  ! define PI
  double precision :: PI = 4.D0*datan(1.D0) !< this is PI

  ! setup grid 
  nx = 500 ! grid dimension x
  ny = 500 ! grid dimension y

  dlx = 1000.D0 ! grid spacing x
  dly = 1000.D0 ! grid spacing y

  ! setup propagating medium
  allocate(vp(nx,ny)) 
  allocate(vs(nx,ny)) 

  allocate(lambda(nx,ny)) 
  allocate(mu(nx,ny)) 
  allocate(rho(nx,ny)) 

  vp(:,:)  = 1000.D0 ! set constant velocity
  vs(:,:)  =  750.D0 ! set constant velocity
  rho(:,:) = 1000.d0 ! set constant density

  mu = vs*vs*rho
  lambda = vp*vp*rho-2*mu

  ! setup source attributes
  point_per_wavelength = 20 ! we want at least 20 points per wavelength to avoid numerical dispersion
  frequency = minval((/vs,vp/))/(point_per_wavelength*max(dlx,dly)) ! setup the central frequency to have enough point per wavelength

  ! setup time parameters
  nt = 600 ! number of time steps 
  dt = min(dlx,dly)/maxval((/vp,vs/))/2. ! dt must satisfy the Courant criterion

  print*, dt,frequency
  ic = 1600
  ! setup source 
  n_src = 1 ! use a single source

  allocate(stf_src_x(nt,n_src))
  allocate(stf_src_y(nt,n_src))

  allocate(i_src(n_src))
  allocate(j_src(n_src))
  do is = 1,n_src
     i_src(is) = nx/2 ! source position x
     j_src(is) = ny/2 ! source position y
  enddo

  ! setup receiver
  n_rec = 1 ! use a single station

  allocate(ux_rec(nt,n_rec))
  allocate(uy_rec(nt,n_rec))

  allocate(i_rec(n_rec))
  allocate(j_rec(n_rec))
  do ir = 1,n_rec
     i_rec(ir) = nx/2 ! station position x
     j_rec(ir) = ny/4   ! station position y
  enddo

  ! setup source time function
  open(1,file='source_x.txt',status='unknown')
  open(2,file='source_y.txt',status='unknown')
  do is = 1,n_src
     call gen_source_time_function(frequency,dt,stf_src_x(:,is)) ! compute source time function
     stf_src_y(:,is) = 0.D0 ! fy = 0 in this example
     do it = 1,nt
        write(1,*)(it-1)*dt,stf_src_x(it,is) ! write source time functions to file
        write(2,*)(it-1)*dt,stf_src_y(it,is) ! write source time functions to file
     enddo
     write(1,*)
     write(2,*)
  enddo
  close(1)
  close(2)

  ! setup mirror window
  i1m = nx/2-nx/10 ! 
  i2m = nx/2+nx/10 ! mirror window bounds
  j1m = ny/2-ny/10 !
  j2m = ny/2+ny/10 !

  allocate(win_mir(nx,ny))
  win_mir(:,:) = 0.D0              ! window is zero outside the mirror
  win_mir(i1m:i2m,j1m:j2m) = 1.D0  ! window is one inside the mirror

  !------------------!
  ! first simulation !
  !------------------!

  !
  ! source is inside the mirror
  ! we recorde the mirror force
  ! This would be the local simulation done with Regsem
  ! however I use a global simulation here
  ! so we get the seismogram p_rec the we want to reproduce at the station

  print*,'test 1'
  call solver_FD2D_elastic(      nt        & 
       ,                         lambda    &
       ,                         mu        &
       ,                         i_src     &
       ,                         j_src     &
       ,                         stf_src_x &
       ,                         stf_src_y &
       ,                         i_rec     &
       ,                         j_rec     &
       ,                         ux_rec    &
       ,                         uy_rec    &
       ,                         n_mir     &
       ,                         1         &
       ,                         win_mir   &
       ,                         i_mir     &
       ,                         j_mir     &
       ,                         fx_mir    &   
       ,                         fy_mir    )

  ! setup the delta function
  allocate(delta_x(nt,1))
  allocate(delta_y(nt,1))
  
  delta_x(:,:) = 0.D0
  delta_y(:,:) = 0.D0
  
  ! Ricker
  f0 = 25.D-2/dt
  open(24,file='green_stf_ricker_y.dat',form='formatted',status='unknown')
  do it = 1,nt
     delta_y(it,1) = rickr((it-10)*dt,f0)
     write(24,*)(it-1)*dt,delta_y(it,1) 
  enddo
  close(24)
  
  ! Convolve with the source time function of the green's function
  ux_rec(:,1) = convolve(ux_rec(:,1), delta_y(:,1), dt)
  uy_rec(:,1) = convolve(uy_rec(:,1), delta_y(:,1), dt)
    
  ! Collect the seismograms 
  open(1,file='ux_rec_ref.dat',status='unknown')
  open(2,file='uy_rec_ref.dat',status='unknown')
  do it = 1,nt
     write(1,*)ux_rec(it,1)
     write(2,*)uy_rec(it,1)
  enddo
  close(1)
  close(2)
  
  ! Here is where the single point test should take place
  allocate(tmp_x(size(fx_mir,1),size(fx_mir,2)))
  allocate(tmp_y(size(fx_mir,1),size(fx_mir,2)))
  tmp_x(:,:) = 0.D0
  tmp_y(:,:) = 0.D0

  ! Now just the one point that we actually care about.
  tmp_x(:,ic) = fx_mir(:,ic)
  tmp_y(:,ic) = fy_mir(:,ic)
 
 ! Uncomment the call below to view the extrapolated wavefield
 ! notice the use of -fx_mir, -fy_mir
 call solver_FD2D_elastic(      nt        & 
      ,                         lambda    &
      ,                         mu        &
      ,                         i_mir     &
      ,                         j_mir     &
      ,                         -tmp_x    &
      ,                         -tmp_y    &
      ,                         i_rec     &
      ,                         j_rec     &
      ,                         ux_rec    &
      ,                         uy_rec    &
      ,                         n_mir     &
      ,                         4         )

  ! store the single point seismogram and fem 
  open(7 ,file='ux_rec_single_point.dat',status='unknown')
  open(8 ,file='uy_rec_single_point.dat',status='unknown')
  open(9 ,file='fem_x_single_point.dat',status='unknown')
  open(10,file='fem_y_single_point.dat',status='unknown')
  do it = 1,nt
     write(7 ,*)(it-1)*dt,ux_rec(it,1)
     write(8 ,*)(it-1)*dt,uy_rec(it,1)
     write(9 ,*)tmp_x(it,ic)
     write(10,*)tmp_y(it,ic)
  enddo
  close(7)
  close(8)
  close(9)
  close(10)
  
  !-------------------!
  ! second simulation !
  !-------------------!

  ! here we compute the green functions
  ! i.e. we record the wavefield at the mirror for a point source excitation

  ! allocate space for Green's functions
  allocate(i_green(size(i_mir)))
  allocate(j_green(size(i_mir)))
  allocate(fxx_green(nt,size(i_mir)))
  allocate(fxy_green(nt,size(i_mir)))
  allocate(fyx_green(nt,size(i_mir)))
  allocate(fyy_green(nt,size(i_mir)))

  ! we want to record at the mirror
  i_green = i_mir
  j_green = j_mir

  ! run the simulation

  ! Green's function x
  
  delta_x(:,:) = 0.D0
  delta_y(:,:) = 0.D0
  
  ! Ricker
  open(23,file='green_stf_ricker_x.dat',form='formatted',status='unknown')
  do it = 1,nt
     delta_x(it,1) = rickr((it-10)*dt,f0)
     write(23,*)(it-1)*dt,delta_x(it,1) 
  enddo
  close(23)

  ! Simple dirac (as simple as it gets)
  ! delta_x(1,1) = 1.D0/dt
  
  print*,'test 2'
  call solver_FD2D_elastic(      nt        & 
       ,                         lambda    &
       ,                         mu        &
       ,                         i_rec     &
       ,                         j_rec     &
       ,                         delta_x   &
       ,                         delta_y   &
       ,                         i_green   &
       ,                         j_green   &
       ,                         fxx_green &
       ,                         fxy_green &
       ,                         n_mir     &
       ,                         2         )
 
  allocate(green_xx(size(ux_rec,1)))
  allocate(green_xy(size(ux_rec,1)))
  allocate(green_yx(size(ux_rec,1)))
  allocate(green_yy(size(ux_rec,1)))

  green_xx(:) = fxx_green(:,ic)
  green_xy(:) = fxy_green(:,ic)

  ! Green's function y
  open(3,file='green_1_ux.dat',status='unknown')
  open(4,file='green_1_uy.dat',status='unknown')
  do it = 1,nt
     write(3,*)green_xx(it)
     write(4,*)green_xy(it)
  enddo
  close(3)
  close(4)
    
  delta_x(:,:) = 0.D0
  delta_y(:,:) = 0.D0

  ! Ricker
  open(24,file='green_stf_ricker_y.dat',form='formatted',status='unknown')
  do it = 1,nt
     delta_y(it,1) = rickr((it-10)*dt,f0)
     write(24,*)(it-1)*dt,delta_y(it,1) 
  enddo
  close(24)

  ! Simple dirac
  ! delta_y(1,1) = 1.D0/dt
  
  call solver_FD2D_elastic(      nt        & 
       ,                         lambda    &
       ,                         mu        &
       ,                         i_rec     &
       ,                         j_rec     &
       ,                         delta_x   &
       ,                         delta_y   &
       ,                         i_green   &
       ,                         j_green   &
       ,                         fyx_green &
       ,                         fyy_green &
       ,                         n_mir     &
       ,                         3         )
       

  green_yx(:) = fyx_green(:,ic)
  green_yy(:) = fyy_green(:,ic)
  
  open(5,file='green_2_ux.dat',status='unknown')
  open(6,file='green_2_uy.dat',status='unknown')
  do it = 1,nt
     write(5,*)green_yx(it)
     write(6,*)green_yy(it)
  enddo
  close(5)
  close(6)
  
  ! now perform wavefield extrappolation
  allocate(ux_rec_extrapolated(nt))
  allocate(uy_rec_extrapolated(nt))
  allocate(ux_single_point(nt))
  allocate(uy_single_point(nt))
  
  ux_rec_extrapolated = 0.D0
  uy_rec_extrapolated = 0.D0
  
  ! open(35,file='mirror.001')
  ! open(36,file='mirror.002')
  ! open(37,file='mirror.003')
  
  ux_single_point = 0.D0
  ux_single_point = 0.D0

  ux_single_point = - convolve(tmp_x(:,ic), green_xx(:), dt)
  ux_single_point = ux_single_point - convolve(tmp_y(:,ic), green_yx(:), dt)
  uy_single_point = - convolve(tmp_x(:,ic), green_xy(:), dt)
  uy_single_point = uy_single_point - convolve(tmp_y(:,ic), green_yy(:), dt)

  do ir = 1,size(i_mir) ! sum convolution for all mirror points
     ux_rec_extrapolated = ux_rec_extrapolated - convolve(fx_mir(:,ir),fxx_green(:,ir),dt)
     ux_rec_extrapolated = ux_rec_extrapolated - convolve(fy_mir(:,ir),fxy_green(:,ir),dt)
  enddo
  
  do ir = 1,size(i_mir) ! sum convolution for all mirror points   
     uy_rec_extrapolated = uy_rec_extrapolated - convolve(fx_mir(:,ir),fyx_green(:,ir),dt)
     uy_rec_extrapolated = uy_rec_extrapolated - convolve(fy_mir(:,ir),fyy_green(:,ir),dt)
  enddo
  close(35)
  close(36)
  close(37)
  
  ! store the extrapolated seismogram
  open(1,file='ux_rec_extrapolated.dat',status='unknown')
  open(2,file='uy_rec_extrapolated.dat',status='unknown')
  open(11,file='ux_rec_sp_extrapol.dat',status='unknown')
  open(12,file='uy_rec_sp_extrapol.dat',status='unknown')
  do it = 1,size(ux_rec_extrapolated)
     write(1,*)ux_rec_extrapolated(it)
     write(2,*)uy_rec_extrapolated(it)
     write(11,*)(it-1)*dt,ux_single_point(it)
     write(12,*)(it-1)*dt,uy_single_point(it)
  enddo
  close(1)
  close(2)
  close(11)
  close(12)

  ! Let's clean up after ourselves, we're not savages.
  deallocate(uy_rec_extrapolated)
  deallocate(ux_rec_extrapolated)
  deallocate(vp)
  deallocate(vs)
  deallocate(lambda)
  deallocate(mu)
  deallocate(rho)
  deallocate(stf_src_x)
  deallocate(stf_src_y)
  deallocate(i_src)
  deallocate(j_src)
  deallocate(ux_rec)
  deallocate(uy_rec)
  deallocate(i_rec)
  deallocate(j_rec)
  deallocate(win_mir)
  deallocate(tmp_x)
  deallocate(tmp_y)

  deallocate(delta_x)
  deallocate(delta_y)
  
  deallocate(i_green)
  deallocate(j_green)
  deallocate(fxx_green)
  deallocate(fxy_green)
  deallocate(fyx_green)
  deallocate(fyy_green)

  deallocate(green_xx)
  deallocate(green_xy)
  deallocate(green_yx)
  deallocate(green_yy)

  !----------!
  ! finished !
  !----------! 

  stop

contains

  !----------------------------------------------------
  subroutine gen_source_time_function(frequency,dt,stf)
  !----------------------------------------------------
    !
    ! compute the source time function
    !
    implicit none
    !
    double precision :: frequency, dt, stf(:)
    double precision :: t_used, aval
    !
    integer :: it
    !
    aval = PI*PI*frequency*frequency
    !
    do it = 1,size(stf)
       t_used = (it-1)*dt-1.2D0/frequency
       stf(it) = -(1.D0-2.D0*aval*t_used**2) * exp(-aval*t_used**2)
    enddo
    !
    return
  !--------------------------------------
  end subroutine gen_source_time_function
  !--------------------------------------

  !-----------------------------------------
  subroutine solver_FD2D_elastic(nt        & 
       ,                         lambda    &
       ,                         mu        &
       ,                         i_src     &
       ,                         j_src     &
       ,                         stf_src_x &
       ,                         stf_src_y &
       ,                         i_rec     &
       ,                         j_rec     &
       ,                         ux_rec    &
       ,                         uy_rec    &
       ,                         n_mir     &
       ,                         optn      &
       ,                         win_mir   &
       ,                         i_mir     &
       ,                         j_mir     &
       ,                         fx_mir    &   
       ,                         fy_mir    )

  !----------------------------------------- 
    !
    implicit none
    !
    ! win_mir, i_mir, j_mir and f_mir are optional arguments
    ! if they are present, the solver will record the mirror force f_mir(i_mir(i),j_mir(j))
    ! for all points i on the mirror
    !
    double precision, optional              :: win_mir(:,:) ! mirror window function
    integer, allocatable, optional          :: i_mir(:)     ! mirror points position indices i
    integer, allocatable, optional          :: j_mir(:)     ! mirror points position indices j
    double precision, allocatable, optional :: fx_mir(:,:)   ! mirror excitation force x
    double precision, allocatable, optional :: fy_mir(:,:)   ! mirror excitation force y
    !
    integer :: nt        ! number of time stemps in the simulation
    integer :: i_src(:)  ! sources position indices i
    integer :: j_src(:)  ! sources position indices j
    integer :: i_rec(:)  ! receivers position indices i
    integer :: j_rec(:)  ! receivers position indices j
    !
    double precision :: lambda(:,:)        ! lame parameter in the propagating medium
    double precision :: mu(:,:)        ! shear modulus in the propagating medium
    double precision :: stf_src_x(:,:)  ! source time functions
    double precision :: stf_src_y(:,:)  ! source time functions
    double precision :: ux_rec(:,:)    ! recordings at the receivers
    double precision :: uy_rec(:,:)    ! recordings at the receivers
    !
    ! wavefield variables
    double precision, dimension(size(mu,1),size(mu,2)) :: fx      ! internal forces x
    double precision, dimension(size(mu,1),size(mu,2)) :: fy      ! internal forces y
    double precision, dimension(size(mu,1),size(mu,2)) :: ux      ! displacement x at the current time step
    double precision, dimension(size(mu,1),size(mu,2)) :: ux_old  ! displacement x at the previous time step
    double precision, dimension(size(mu,1),size(mu,2)) :: ux_new  ! displacement x at the next time step
    double precision, dimension(size(mu,1),size(mu,2)) :: uy      ! displacement y at the current time step
    double precision, dimension(size(mu,1),size(mu,2)) :: uy_old  ! displacement y at the previous time step
    double precision, dimension(size(mu,1),size(mu,2)) :: uy_new  ! displacement y at the next time step
    double precision, dimension(size(mu,1),size(mu,2)) :: fxm     ! mirror excitation force x (computed at all grid points)
    double precision, dimension(size(mu,1),size(mu,2)) :: fym     ! mirror excitation force y (computed at all grid points)
    !
    ! Array for storing the displacement in 2D at all points of the mirror
    double precision, allocatable :: total_disp(:,:),total_disp_spln(:,:,:)
    double precision, allocatable :: diag(:,:),ddiag(:)      
    integer :: count_rec_mir
   
    ! Parameters for the B-spline compression
    integer :: spln_order_mir = 4
    double precision, allocatable :: bmn_coeff(:)
    integer :: decim_fact_mir
    integer :: rec_len_mir
    character :: fname_mirror*100
    integer :: optn

    double precision, dimension(size(mu,1),size(mu,2)) :: map     ! for plotting the wavefield
    !
    integer :: i,j,ii,it,is,ir,n_mir ! looping indices and counters
    !
    ! find mirror points if win_mir is present
    !
    
    if(present(win_mir))then
       call force_int_elastic(lambda,mu,dlx,dly,win_mir,win_mir,fx,fy) ! when (fx,fy)/=0 the we want to record the mirror
       n_mir = 0
       do j = 1,size(mu,2)
          do i = 1,size(mu,1)
             if(abs(fx(i,j))>=1.d-12.or.abs(fy(i,j))>=1.d-12)then ! I have had precision issues
                n_mir = n_mir+1 ! count points on the mirror
             endif
          enddo
       enddo
       allocate(i_mir(n_mir)) ! allocate space for mirror
       allocate(j_mir(n_mir))
       allocate(fx_mir(nt,n_mir))
       allocate(fy_mir(nt,n_mir))
       n_mir = 0
       do j = 1,size(mu,2)
          do i = 1,size(mu,1)
             if(abs(fx(i,j))>=1.d-12.or.abs(fy(i,j))>=1.d-12)then ! I have had precision issues
                n_mir = n_mir+1 ! store the positions of the mirror points
                i_mir(n_mir) = i
                j_mir(n_mir) = j
             endif
          enddo
       enddo
    endif
    !
    ! open graphic device
    call pgopen('/xwin')
    call pgsvp(0.,1.,0.,1.)
    call pgwnad(0.5,size(mu,1)+0.5,0.5,size(mu,2)+0.5)
    !

    ux_old(:,:) = 0.D0 ! init field variables
    ux(:,:)     = 0.D0
    ux_new(:,:) = 0.D0
    fx(:,:)     = 0.D0
    uy_old(:,:) = 0.D0 
    uy(:,:)     = 0.D0
    uy_new(:,:) = 0.D0
    fy(:,:)     = 0.D0
    !
  
    do it = 1,nt ! time loop
       !
       ! write(*,*)it,minval(ux),maxval(ux),minval(uy),maxval(uy)
       ! 
       ! Constructing the array for the global wavefield
       if (it==1)then
          decim_fact_mir = 10 
          count_rec_mir = 0
          print*,n_mir,decim_fact_mir
          allocate(bmn_coeff((spln_order_mir+1)*decim_fact_mir+1))
          call bmn(bmn_coeff ,decim_fact_mir,spln_order_mir)
          
          allocate(total_disp(2,n_mir))
          allocate(total_disp_spln(2,n_mir,spln_order_mir+1))
          
          inquire(iolength=rec_len_mir)total_disp_spln(:,:,1)
          print*,n_mir
          write(fname_mirror,"(a,I3.3)")"mirror.",optn
          open(33,file=fname_mirror,access='direct',form="unformatted",&
               recl=rec_len_mir,status='replace')
       endif

       ! Record wavefield
       do ir = 1,size(i_rec)
          ux_rec(it,ir) = ux(i_rec(ir),j_rec(ir))
          uy_rec(it,ir) = uy(i_rec(ir),j_rec(ir))
       enddo       
          
       ! compute/record mirror force if needed (eq.55 in Masson et al., 2018)
       if(present(win_mir))then
          !
          call force_int_elastic(lambda,mu,dlx,dly,ux,uy,fx,fy)
          !
          fxm =win_mir*fx
          fym =win_mir*fy
          !
          call force_int_elastic(lambda,mu,dlx,dly,ux*win_mir,uy*win_mir,fx,fy)
          !
          fxm = fxm-fx
          fym = fym-fy
          !
          do ir = 1,n_mir
             fx_mir(it,ir) = fxm(i_mir(ir),j_mir(ir))
             fy_mir(it,ir) = fym(i_mir(ir),j_mir(ir))
          enddo
       endif
       !
       ! compute "internal forces" (i.e. right hand side in eq.55; Masson et al., 2018)
       call force_int_elastic(lambda,mu,dlx,dly,ux,uy,fx,fy)
       !
       ! add sources
       do is = 1,size(i_src)
          fx(i_src(is),j_src(is)) = fx(i_src(is),j_src(is))+stf_src_x(it,is)
          fy(i_src(is),j_src(is)) = fy(i_src(is),j_src(is))+stf_src_y(it,is)
       enddo
       !
       ux_new = 2*ux-ux_old+dt*dt/rho*fx ! time integration
       uy_new = 2*uy-uy_old+dt*dt/rho*fy ! time integration
       !
       ux_old = ux ! update the time variables
       ux = ux_new
       uy_old = uy
       uy = uy_new
       !
       if(present(win_mir))then
          do ii = 1,n_mir
             total_disp(1,ii) = ux(i_mir(ii), j_mir(ii))
             total_disp(2,ii) = uy(i_mir(ii), j_mir(ii))
          enddo
       else
          do ii = 1,size(i_rec)
             total_disp(1,ii) = ux(i_rec(ii), j_rec(ii))
             total_disp(1,ii) = uy(i_rec(ii), j_rec(ii))
          enddo
       endif

       ! Now the spline compression
       do j = 1,spln_order_mir+1
          i = mod(it-1,decim_fact_mir)+(spln_order_mir+1-j)*decim_fact_mir+1
          total_disp_spln(:,:,j) = total_disp_spln(:,:,j) + bmn_coeff(i)*total_disp(:,:)
       enddo
       
       ! Decimate 
       if(mod(it-1,decim_fact_mir) == decim_fact_mir-1)then
          count_rec_mir = count_rec_mir+1
          print*,count_rec_mir
          write(33,rec=count_rec_mir)total_disp_spln(:,:,1)
          do j = 1,spln_order_mir
             total_disp_spln(:,:,j) = total_disp_spln(:,:,j+1)
          enddo
          total_disp_spln(:,:,spln_order_mir+1) = 0
       endif
        
       if (it==nt)then
          do j = 1,spln_order_mir
             count_rec_mir = count_rec_mir+1
             write(33,rec=count_rec_mir)total_disp_spln(:,:,j)
          enddo
       
          ! obtain B-spline coefficients
          deallocate(total_disp_spln)
          allocate(diag(spln_order_mir+1,count_rec_mir))
          allocate(ddiag(count_rec_mir))
          allocate(total_disp_spln(2,n_mir,0:spln_order_mir))
          call fill_spmat_diag(decim_fact_mir,spln_order_mir,count_rec_mir,nt,diag)
          call bchfac(diag, spln_order_mir+1, count_rec_mir, ddiag)
          call bchslv(total_disp_spln,2,n_mir,spln_order_mir,diag,spln_order_mir+1,count_rec_mir,33)
       
          deallocate(diag)
          deallocate(ddiag)
          deallocate(total_disp)
          deallocate(total_disp_spln)
       endif
       ! below is purely graphics
       !
       ! plot wavefield
       map = ux**2+uy**2
       call pgimag(sngl(map),size(map,1),size(map,2),1,size(map,1),1,size(map,2), &
            -maxval(sngl(abs(map))),maxval(sngl(abs(map))),(/0.,1.,0.,0.,0.,1./))
       !
       ! plot mirror window if present
       ! 
       if(present(win_mir))then
          do i = 1,size(mu,1)-1
             do j = 1,size(mu,2)
                if(win_mir(i,j)/=win_mir(i+1,j))then
                   call pgmove(real(i)+0.5,real(j)-0.5)
                   call pgdraw(real(i)+0.5,real(j+1)-0.5)
                endif
             enddo 
          enddo
          do i = 1,size(mu,1)
             do j = 1,size(mu,2)-1
                if(win_mir(i,j)/=win_mir(i,j+1))then
                   call pgmove(real(i)-0.5,real(j)+0.5)
                   call pgdraw(real(i+1)-0.5,real(j)+0.5)
                endif
             enddo
          enddo
       endif
       !
       ! plot sources
       do i = 1,size(i_src)
          call pgpt1(real(i_src(i))+0.5,real(j_src(i))-0.5,6)
       enddo
       !
       ! plot receivers
       if(present(win_mir))then
          do i = 1,size(i_rec)
             call pgpt1(real(i_rec(i))+0.5,real(j_rec(i))-0.5,5)
          enddo
       else
          do i = 1,size(i_rec)
             call pgpt1(real(i_rec(i)),real(j_rec(i)),-2)
             !if(i_rec(i+i)==i_rec(i)) then
              !  call pgdraw(real(i_rec(i+1)),real(j_rec(i+1)))
             !endif
          enddo
       endif
       !
    enddo
    close(33)
    deallocate(bmn_coeff)

    !
    call pgclos
    !
    return
    !
  !---------------------------------
  end subroutine solver_FD2D_elastic
  !---------------------------------
  
  !------------------------------------------------------
  subroutine force_int_elastic(lambda,mu,dlx,dly,ux,uy,fx,fy)
  !------------------------------------------------------
    !
    implicit none
    !
    double precision :: lambda(:,:),mu(:,:),dlx,dly
    double precision :: ux(:,:),uy(:,:),fx(:,:), fy(:,:)
    double precision, dimension(size(ux,1),size(ux,2)) :: txx, tyy, txy
    !
    fx(:,:) = 0.D0
    fy(:,:) = 0.D0
    !
    txx = (lambda+2*mu)*dx(ux) + lambda*dy(uy)
    tyy = (lambda+2*mu)*dy(uy) + lambda*dx(ux) 
    txy = mu*(dy(ux) + dx(uy))
    !
    fx = dxt(txx) + dyt(txy)
    fy = dyt(tyy) + dxt(txy)
    !
    return
    !
  !-------------------------------
  end subroutine force_int_elastic
  !-------------------------------

  !-------------
  function dx(x)
  !-------------
    !
    implicit none
    !
    double precision :: x(:,:),dx(size(x,1),size(x,2))
    integer :: j,nx,ny
    !
    nx = size(x,1)
    ny = size(x,2)
    !
    do j = 1,ny
       dx(:,j) = dev(x(:,j))/dlx
    enddo
    !
    return
    !
  !--------------
  end function dx
  !--------------

  !--------------
  function dy(x)
  !--------------
    !
    implicit none
    !
    double precision :: x(:,:),dy(size(x,1),size(x,2))
    integer :: i,nx,ny
    !
    nx = size(x,1)
    ny = size(x,2)
    !
    do i = 1,nx
       dy(i,:) = dev(x(i,:))/dly
    enddo
    !
    return
    !
  !--------------
  end function dy
  !--------------
  
  !--------------
  function dxt(x)
  !--------------
    !
    implicit none
    !
    double precision :: x(:,:),dxt(size(x,1),size(x,2))
    integer :: j,nx,ny
    !
    nx = size(x,1)
    ny = size(x,2)
    !
    do j = 1,ny
       dxt(:,j) = devt(x(:,j))/dlx
    enddo
    !
    return
    !
  !---------------
  end function dxt
  !--------------
  
  !--------------
  function dyt(x)
  !--------------
    !
    implicit none
    !
    double precision :: x(:,:),dyt(size(x,1),size(x,2))
    integer :: i,nx,ny
    !
    nx = size(x,1)
    ny = size(x,2)
    !
    do i = 1,nx
       dyt(i,:) = devt(x(i,:))/dly
    enddo
    !
    return
    !
  !---------------
  end function dyt
  !--------------

  !--------------
  function dev(x)
  !--------------
    !
    implicit none
    !
    double precision :: x(:),dev(size(x))
    integer :: i,n
    double precision, parameter :: c1 = -1.d0/4.D0
    double precision, parameter :: c2 = -5.d0/6.D0
    double precision, parameter :: c3 = +3.d0/2.D0
    double precision, parameter :: c4 = -1.d0/2.D0
    double precision, parameter :: c5 = +1.d0/12.D0
    !
    n = size(x)
    dev(:) = 0.D0
    !
    do i = 4,size(x)-3
       dev(i) = c1*x(i-1)+c2*x(i)+c3*x(i+1)+c4*x(i+2)+c5*x(i+3)
    enddo
    !
    return
    !
  !---------------
  end function dev
  !---------------

  
  !---------------
  function devt(x)
  !---------------
    !
    implicit none
    !
    double precision :: x(:),devt(size(x))
    integer :: i,n
    double precision, parameter :: c1 = -1.d0/12.D0
    double precision, parameter :: c2 = +1.d0/2.D0
    double precision, parameter :: c3 = -3.d0/2.D0
    double precision, parameter :: c4 = +5.d0/6.D0
    double precision, parameter :: c5 = +1.d0/4.D0
    !
    n = size(x)
    devt(:) = 0.D0
    !
    do i = 4,size(x)-3
       devt(i) = c1*x(i-3)+c2*x(i-2)+c3*x(i-1)+c4*x(i)+c5*x(i+1)
    enddo
    !
    return
    !
  !----------------
  end function devt
  !----------------
  
  !------------------------
  function convolve(a,b,dt)
  !------------------------
    !
    ! convolve a with b 
    ! keep only the causal part
    ! convolution is cut so it has the same duration as a
    !
    implicit none
    !
    double precision :: a(:),b(:)
    double precision :: convolve(size(a))
    double precision :: buffer(size(a)+size(b))
    double precision :: dt
    integer :: i,na,nb
    !
    na = size(a)
    nb = size(b)
    !
    buffer(:) = 0.D0
    !
    do i = 1,size(a)
       ! for all time steps t in a we add the delayed response b scaled by the value of a(t)
       buffer(i:i+nb-1) = buffer(i:i+nb-1)+a(i)*b
    enddo
    !
    convolve(:) = buffer(1:na)*dt
    !
    return
    !
    !--------------------
  end function convolve
  !--------------------

  
  !====================
  subroutine bmn(b,m,n)
    !====================
    !
    implicit none
    !
    integer :: m,n,i
    double precision :: b((n+1)*m+1),bspln
    !
    do i = 1,(n+1)*m+1
       b(i) = bspln(0,n,dble(i-1)/dble((n+1)*m)*dble(n+1))
    enddo
    !
    return
    !=================
  end subroutine bmn
  !=================
!===========================================
subroutine fill_spmat_diag(m,n,nsp,nt,spmat)
!===========================================
!
implicit none
!
integer :: m,n,nsp,nt
double precision :: spmat(n+1,nsp)
double precision :: b((n+1)*m+1)
!
integer :: i,j,i1,i2,i3,i4
!
call bmn(b,m,n)
!
do j = 1,nsp
   do i = 1,n+1
!
      i1 = max(1+(i-1)*m,(n+1-j)*m+1)
      i2 = min(1+(n+1)*m,(n-j+1)*m+nt)
      i3 = max(1,(n+2-j-i)*m+1)
      i4 = min(1+(n+2-i)*m,(n+2-j-i)*m+nt)
!
      if(j+i-1<=nsp)then
         spmat(i,j) = dot_product(b(i1:i2),b(i3:i4))
      else
         spmat(i,j) = 0.d0
      endif
!
   enddo
enddo
!
return
!=============================
end subroutine fill_spmat_diag
!=============================  

!==========================================
subroutine bchfac ( w, nbands, nrow, diag )
!==========================================

!*****************************************************************************80
!
!! BCHFAC constructs a Cholesky factorization of a matrix.
!
!  Discussion:
!
!    The factorization has the form
!
!      C = L * D * L'
!  
!    with L unit lower triangular and D diagonal, for a given matrix C of 
!    order NROW, where C is symmetric positive semidefinite and banded, 
!    having NBANDS diagonals at and below the main diagonal.
! 
!    Gauss elimination is used, adapted to the symmetry and bandedness of C.
! 
!    Near-zero pivots are handled in a special way.  The diagonal 
!    element C(N,N) = W(1,N) is saved initially in DIAG(N), all N. 
! 
!    At the N-th elimination step, the current pivot element, W(1,N), 
!    is compared with its original value, DIAG(N).  If, as the result 
!    of prior elimination steps, this element has been reduced by about 
!    a word length, that is, if W(1,N) + DIAG(N) <= DIAG(N), then the pivot 
!    is declared to be zero, and the entire N-th row is declared to
!    be linearly dependent on the preceding rows.  This has the effect 
!    of producing X(N) = 0 when solving C * X = B for X, regardless of B.
! 
!    Justification for this is as follows.  In contemplated applications 
!    of this program, the given equations are the normal equations for 
!    some least-squares approximation problem, DIAG(N) = C(N,N) gives 
!    the norm-square of the N-th basis function, and, at this point, 
!    W(1,N) contains the norm-square of the error in the least-squares 
!    approximation to the N-th basis function by linear combinations 
!    of the first N-1.  
!
!    Having W(1,N)+DIAG(N) <= DIAG(N) signifies that the N-th function 
!    is linearly dependent to machine accuracy on the first N-1 
!    functions, therefore can safely be left out from the basis of 
!    approximating functions.
!
!    The solution of a linear system C * X = B is effected by the 
!    succession of the following two calls:
! 
!      call bchfac ( w, nbands, nrow, diag )
!
!      call bchslv ( w, nbands, nrow, b, x )
!
!  Modified:
!
!    14 February 2007
!
!  Author:
!
!    Carl DeBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) W(NBANDS,NROW).
!    On input, W contains the NBANDS diagonals in its rows, 
!    with the main diagonal in row 1.  Precisely, W(I,J) 
!    contains C(I+J-1,J), I=1,...,NBANDS, J=1,...,NROW.
!    For example, the interesting entries of a seven diagonal
!    symmetric matrix C of order 9 would be stored in W as
! 
!      11 22 33 44 55 66 77 88 99
!      21 32 43 54 65 76 87 98  *
!      31 42 53 64 75 86 97  *  *
!      41 52 63 74 85 96  *  *  *
!
!    Entries of the array not associated with an
!    entry of C are never referenced.
!    On output, W contains the Cholesky factorization 
!    C = L*D*L', with W(1,I) containing 1/D(I,I) and W(I,J) 
!    containing L(I-1+J,J), I=2,...,NBANDS.
!
!    Input, integer ( kind = 4 ) NBANDS, indicates the bandwidth of the
!    matrix C, that is, C(I,J) = 0 for NBANDS < abs(I-J).
! 
!    Input, integer ( kind = 4 ) NROW, is the order of the matrix C.
! 
!    Work array, real ( kind = 8 ) DIAG(NROW).
!
  implicit none

  integer nbands
  integer nrow

  double precision diag(nrow)
  integer i
  integer imax
  integer j
  integer jmax
  integer n
  double precision ratio
  double precision w(nbands,nrow)

  if ( nrow <= 1 ) then
    if ( 0.0D+00 < w(1,1) ) then
      w(1,1) = 1.0D+00 / w(1,1)
    end if
    return
  end if
!
!  Store the diagonal.
!
  diag(1:nrow) = w(1,1:nrow)
!
!  Factorization.
!
  do n = 1, nrow
 
    if ( w(1,n) + diag(n) <= diag(n) ) then
      w(1:nbands,n) = 0.0D+00
    else
 
      w(1,n) = 1.0D+00 / w(1,n)
 
      imax = min ( nbands - 1, nrow - n )
 
      jmax = imax
 
      do i = 1, imax
 
        ratio = w(i+1,n) * w(1,n)
 
        do j = 1, jmax
          w(j,n+i) = w(j,n+i) - w(j+i,n) * ratio
        end do
 
        jmax = jmax-1
        w(i+1,n) = ratio
 
      end do
 
    end if
 
  end do
 
  return

!==
end
!==


!==========================================================
subroutine bchslv (tmp, n1, n2, n3, w, nbands, nrow, lunit)
!==========================================================

!*****************************************************************************
!
!! BCHSLV solves a banded symmetric positive definite system.
!
!  Discussion:
!
!    The system is of the form:
!
!      C * X = B 
!  
!    and the Cholesky factorization of C has been constructed 
!    by BCHFAC.
! 
!    With the factorization 
!
!      C = L * D * L'
!
!    available, where L is unit lower triangular and D is diagonal, 
!    the triangular system 
!
!      L * Y = B 
!
!    is solved for Y (forward substitution), Y is stored in B, the 
!    vector D**(-1)*Y is computed and stored in B, then the 
!    triangular system L'*X = D**(-1)*Y is solved for X 
!    (back substitution).
!
!  Modified:
!
!    14 February 2007
!
!  Author:
!
!    Carl DeBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) W(NBANDS,NROW), the Cholesky factorization for C, 
!    as computed by BCHFAC.
! 
!    Input, integer ( kind = 4 ) NBANDS, the bandwidth of C.
!
!    Input, integer ( kind = 4 ) NROW, the order of the matrix C.
! 
!    Input/output, real ( kind = 8 ) B(NROW).
!    On input, the right hand side.
!    On output, the solution.
!
  
  implicit none

  integer :: j
  integer :: n,n1,n2,n3
  integer :: nrow,nbands,lunit
  double precision :: w(nbands,nrow)
  ! character :: opt*5      
  double precision :: tmp(n1,n2,0:n3)   
! 
!nbands = mir%spln_order+1, nrow =  mir%count_rec
!
!  allocate(mir%tmp(0:mir%recl_mirror-1,0:mir%spln_order))
!
  if ( nrow <= 1 ) then
    write(*,*)'warning few time steps after decimation'
    return
  end if
!
!  Forward substitution. 
!  Solve L*Y = B.
!
  do n = 1, nrow-nbands+1
     if(n==1)then
        do j = 0, nbands - 1
           read(lunit,rec=j+n)tmp(:,:,j)
        enddo
     else
        do j = 0,nbands-2
           tmp(:,:,j) = tmp(:,:,j+1)
        enddo
         read(lunit,rec=n+nbands-1)tmp(:,:,nbands-1)
     endif
     do j = 1, nbands - 1
        tmp(:,:,j) = tmp(:,:,j) - w(j+1,n) * tmp(:,:,0)
     end do
     
     write(lunit,rec=n)tmp(:,:,0)
    
 end do

  do n = nrow-nbands+2, nrow

     do j = 0,nbands-2
        tmp(:,:,j) = tmp(:,:,j+1)
     enddo

     do j = 1, nrow - n
        tmp(:,:,j) = tmp(:,:,j) - w(j+1,n) * tmp(:,:,0)
     end do
     
     write(lunit,rec=n)tmp(:,:,0)
     
  end do
!
!  Back substitution. 
!  Solve L'*X = D**(-1)*Y.
!
  do n = nrow, nrow-nbands+2, -1

     do j = nrow-n,1,-1
        tmp(:,:,j) = tmp(:,:,j-1)
     enddo
     read(lunit,rec=n)tmp(:,:,0)
     
    tmp(:,:,0) = tmp(:,:,0) * w(1,n)

    do j = 1, nrow - n
      tmp(:,:,0) = tmp(:,:,0) - w(j+1,n) * tmp(:,:,j)
    end do

     write(lunit,rec=n)tmp(:,:,0)

  end do

  do n = nrow-nbands+1, 1, -1

     do j = nbands-1,1,-1
        tmp(:,:,j) = tmp(:,:,j-1)
     enddo
     read(lunit,rec=n)tmp(:,:,0)
     
    tmp(:,:,0) = tmp(:,:,0) * w(1,n)

    do j = 1, nbands - 1
      tmp(:,:,0) = tmp(:,:,0) - w(j+1,n) * tmp(:,:,j)
    end do

     write(lunit,rec=n)tmp(:,:,0)

  end do

!  deallocate(mir%tmp)
  return

!====================
end subroutine bchslv
!====================

  double precision function rickr(t,f0)
  
  implicit none
  
  ! double precision :: rickr
  double precision t,f0

  ! Ricker
  rickr = (1.d0 - 2.d0*PI*PI*f0*f0*t*t ) &
       * exp( -PI*PI*f0*f0*t*t )
  
  end function rickr


  !
  !---------------
end program main
!---------------


!========================================
  recursive function bspln(i,k,x) result(b)
    !========================================
    
    implicit none
    
    double precision :: x,b
    integer :: k,i
    
    b=0.d0
    if(k+1>0)then
       if(k+1==1)then
          if(x>=i.and.x<i+1)b=1.d0   
       else
          b=(x-i)*bspln(i,k-1,x)/(k+1-1)+(i+k+1-x)*bspln(i+1,k-1,x)/(k+1-1)
          if(k==0)b=0
       endif
    endif
    
    !=================
  end function bspln
  !=================
  
