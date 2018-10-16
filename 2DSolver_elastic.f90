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

  ! define PI

  double precision :: PI = 4.D0*datan(1.D0) !< this is PI

  ! setup grid 

  nx = 500 ! grid dimension x
  ny = 500 ! grid dimension y

  dlx = 0.25D0 ! grid spacing x
  dly = 0.25D0 ! grid spacing y

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

  ! setup source 

  n_src = 1 ! use a single source

  allocate(stf_src_x(nt,n_src))
  allocate(stf_src_y(nt,n_src))

  allocate(i_src(n_src))
  allocate(j_src(n_src))
  do is = 1,n_src
     i_src(is) = nx/3 ! source position x
     j_src(is) = ny/3 ! source position y
  enddo

  ! setup receiver

  n_rec = 1 ! use a single station

  allocate(ux_rec(nt,n_rec))
  allocate(uy_rec(nt,n_rec))

  allocate(i_rec(n_rec))
  allocate(j_rec(n_rec))
  do ir = 1,n_rec
     i_rec(ir) = 1*nx/2 ! station position x
     j_rec(ir) = 1*ny/2   ! station position y
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
       ,                         win_mir   &
       ,                         i_mir     &
       ,                         j_mir     &
       ,                         fx_mir    &   
       ,                         fy_mir    )

  ! Uncomment the call below to view the extrapolated wavefield
  ! notice the use of -fx_mir, -fy_mir
 ! call solver_FD2D_elastic(      nt        & 
 !      ,                         lambda    &
 !      ,                         mu        &
 !      ,                         i_mir     &
 !      ,                         j_mir     &
 !      ,                        -fx_mir    &
 !      ,                        -fy_mir    &
 !      ,                         i_rec     &
 !      ,                         j_rec     &
 !      ,                         ux_rec    &
 !      ,                         uy_rec    )

  ! store the reference seismogram

  open(1,file='ux_rec_ref.dat',status='unknown')
  open(2,file='uy_rec_ref.dat',status='unknown')
  do it = 1,nt
     write(1,*)ux_rec(it,1)
     write(2,*)uy_rec(it,1)
  enddo
  close(1)
  close(2)

  !-------------------!
  ! second simulation !
  !-------------------!

  ! here we compute the green functions
  ! i.e. we record the wavefield at the mirror for a point source excitation

  ! setup the delta function
  allocate(delta_x(nt,1))
  allocate(delta_y(nt,1))

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
  delta_x(1,1) = 1.D0/dt
  
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
       ,                         fxy_green )
  
  ! Green's function y
  
  delta_x(:,:) = 0.D0
  delta_y(:,:) = 0.D0
  delta_y(1,1) = 1.D0/dt
  
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
       ,                         fyy_green )


  ! now perform wavefield extrappolation

  allocate(ux_rec_extrapolated(nt))
  allocate(uy_rec_extrapolated(nt))
  
  ux_rec_extrapolated = 0.D0
  uy_rec_extrapolated = 0.D0
  
  do ir = 1,size(i_mir) ! sum convolution for all mirror points
     ux_rec_extrapolated = ux_rec_extrapolated - convolve(fx_mir(:,ir),fxx_green(:,ir),dt)
     ux_rec_extrapolated = ux_rec_extrapolated - convolve(fy_mir(:,ir),fyx_green(:,ir),dt)
  enddo
  
  do ir = 1,size(i_mir) ! sum convolution for all mirror points
     uy_rec_extrapolated = uy_rec_extrapolated - convolve(fx_mir(:,ir),fxy_green(:,ir),dt)
     uy_rec_extrapolated = uy_rec_extrapolated - convolve(fy_mir(:,ir),fyy_green(:,ir),dt)
  enddo

  ! store the extrapolated seismogram
  open(1,file='ux_rec_extrapolated.dat',status='unknown')
  open(2,file='uy_rec_extrapolated.dat',status='unknown')
  do it = 1,size(ux_rec_extrapolated)
     write(1,*)ux_rec_extrapolated(it)
     write(2,*)uy_rec_extrapolated(it)
  enddo
  close(1)
  close(2)

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
    double precision, dimension(size(mu,1),size(mu,2)) :: map     ! for plotting the wavefield
    !
    integer :: i,j,it,is,ir,n_mir ! looping indices and counters
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
       ! below is purely graphics
       !
       ! plot wavefield
       map = ux**2+uy**2
       call pgimag(sngl(map),size(map,1),size(map,2),1,size(map,1),1,size(map,2), &
            -maxval(sngl(abs(map))),maxval(sngl(abs(map))),(/0.,1.,0.,0.,0.,1./))
       !
       ! plot mirror window if present
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
       do i = 1,size(i_rec)
          call pgpt1(real(i_rec(i))+0.5,real(j_rec(i))-0.5,5)
       enddo
       !
    enddo
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
    integer :: i,j
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
    integer :: i,j,nx,ny
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
    integer :: i,j,nx,ny
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
    integer :: i,j,nx,ny
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
    integer :: i,j,nx,ny
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
    integer :: i,j,na,nb
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
  !
  !---------------
end program main
!---------------


