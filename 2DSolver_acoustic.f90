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
  double precision :: dx   ! grid spacing x
  double precision :: dy   ! grid spacing y

  ! physical properties of the propagating medium

  double precision, allocatable :: c(:,:) ! sound speed

  ! source parameters

  integer :: n_src ! number of sources

  integer, allocatable :: i_src(:) ! source position x
  integer, allocatable :: j_src(:) ! source position y

  integer :: point_per_wavelength  ! minimum number of point per wavelength 
  double precision :: frequency    ! central frequency of source time function

  double precision, allocatable :: stf_src(:,:) ! Source time functions (for all sources)

  ! receiver parameters

  integer :: n_rec ! number of receivers

  integer, allocatable :: i_rec(:) ! receiver position x
  integer, allocatable :: j_rec(:) ! receiver position y

  double precision, allocatable :: p_rec(:,:) ! for storing pressure seismograms

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
  double precision, allocatable :: f_mir(:,:) ! Mirror force Fm, for all mirror points for all time steps

  ! Green functions parameters

  integer, allocatable :: i_green(:) ! position indices i of mirror points (where Fm/=0)
  integer, allocatable :: j_green(:) ! position indices j of mirror points (where Fm/=0)
  double precision, allocatable :: f_green(:,:) ! To store green's functions
  double precision, allocatable :: p_rec_extrapolated(:) ! extrapolated pressure seismograms

  ! array for storing the delta function (for computing the Green's functions)

  double precision, allocatable :: delta(:,:) ! it is 2D for consistency with the solverbut really a 1D array

  ! define PI

  double precision :: PI = 4.D0*datan(1.D0) !< this is PI

  ! setup grid 

  nx = 500 ! grid dimension x
  ny = 500 ! grid dimension y

  dx = 0.25D0 ! grid spacing x
  dy = 0.25D0 ! grid spacing y

  ! setup propagating medium

  allocate(c(nx,ny)) 
  c(:,:) = 1000.D0 ! set constant velocity

  ! setup source attributes

  point_per_wavelength = 20 ! we want at least 20 points per wavelength to avoid numerical dispersion
  frequency = minval(c)/(point_per_wavelength*max(dx,dy)) ! setup the central frequency to have enough point per wavelength

  ! setup time parameters

  nt = 600 ! number of time steps 
  dt = min(dx,dy)/maxval(c)/2. ! dt must satisfy the Courant criterion

  ! setup source 

  n_src = 1 ! use a single source

  allocate(stf_src(nt,n_src))

  allocate(i_src(n_src))
  allocate(j_src(n_src))
  do is = 1,n_src
     i_src(is) = nx/2 ! source position x
     j_src(is) = ny/2 ! source position y
  enddo

  ! setup receiver

  n_rec = 1 ! use a single station

  allocate(p_rec(nt,n_rec))

  allocate(i_rec(n_rec))
  allocate(j_rec(n_rec))
  do ir = 1,n_rec
     i_rec(ir) = 3*nx/4 ! station position x
     j_rec(ir) = ny/2   ! station position y
  enddo

  ! setup source time function

  open(1,file='source.txt',status='unknown')
  do is = 1,n_src
     call gen_source_time_function(frequency,dt,stf_src(:,is)) ! compute source time function
     do it = 1,nt
        write(1,*)(it-1)*dt,stf_src(it,is) ! write source time functions to file
     enddo
     write(1,*)
  enddo
  close(1)

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

  call solver_FD2D_acoustic(nt,c,i_src,j_src,stf_src,&
       i_rec,j_rec,p_rec,win_mir,i_mir,j_mir,f_mir)

  ! store the reference seismogram

  open(1,file='p_rec_ref.dat',status='unknown')
  do it = 1,nt
     write(1,*)p_rec(it,1)
  enddo
  close(1)

  !-------------------!
  ! second simulation !
  !-------------------!

  ! here we compute the green functions
  ! i.e. we record the wavefield at the mirror for a point source excitation

  ! setup the delta function
  allocate(delta(nt,1))
  delta(:,:)  = 0
  delta(1,1) = 1./dt ! the integral of delta versus time must be one

  ! allocate space for Green's functions
  allocate(i_green(size(i_mir)))
  allocate(j_green(size(i_mir)))
  allocate(f_green(nt,size(i_mir)))

  ! we want to record at the mirror
  i_green = i_mir
  j_green = j_mir

  ! run the simulation
  call solver_FD2D_acoustic(nt,c,i_rec,j_rec,delta,&
       i_green,j_green,f_green)

  ! now perform wavefield extrappolation

  allocate(p_rec_extrapolated(nt))
  p_rec_extrapolated = 0.D0
  do ir = 1,size(i_mir) ! sum convolution for all mirror points
     p_rec_extrapolated = p_rec_extrapolated - convolve(f_mir(:,ir),f_green(:,ir),dt)
  enddo

  ! store the extrapolated seismogram
  open(1,file='p_rec_extrapolated.dat',status='unknown')
  do it = 1,size(p_rec_extrapolated)
     write(1,*)p_rec_extrapolated(it)
  enddo
  close(1)

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
       stf_src(it,n_src) = -(1.D0-2.D0*aval*t_used**2) * exp(-aval*t_used**2)
    enddo
    !
    return
  !--------------------------------------
  end subroutine gen_source_time_function
  !--------------------------------------

  !----------------------------------------
  subroutine solver_FD2D_acoustic(nt      & 
       ,                          c       &
       ,                          i_src   &
       ,                          j_src   &
       ,                          stf_src &
       ,                          i_rec   &
       ,                          j_rec   &
       ,                          p_rec   &
       ,                          win_mir &
       ,                          i_mir   &
       ,                          j_mir   &
       ,                          f_mir   )     
  !---------------------------------------- 
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
    double precision, allocatable, optional :: f_mir(:,:)   ! mirror excitation force
    !
    integer :: nt        ! number of time stemps in the simulation
    integer :: i_src(:)  ! sources position indices i
    integer :: j_src(:)  ! sources position indices j
    integer :: i_rec(:)  ! receivers position indices i
    integer :: j_rec(:)  ! receivers position indices j
    !
    double precision :: c(:,:)        ! sound speed in the propagating medium
    double precision :: stf_src(:,:)  ! source time functions
    double precision :: p_rec(:,:)    ! recordings at the receivers
    !
    ! wavefield variables
    double precision, dimension(size(c,1),size(c,2)) :: f      ! internal forces
    double precision, dimension(size(c,1),size(c,2)) :: p      ! pressure field at the current time step
    double precision, dimension(size(c,1),size(c,2)) :: p_old  ! pressure field at the previous time step
    double precision, dimension(size(c,1),size(c,2)) :: p_new  ! pressure field at the next time step
    double precision, dimension(size(c,1),size(c,2)) :: fm     ! mirror excitation force computed at all locations (not only at the mirror points)      
    !
    integer :: i,j,it,is,ir,n_mir ! looping indices and counters
    !
    ! find mirror points if win_mir is present
    !
    if(present(win_mir))then
       f = force_int(dble(win_mir),c,dx,dy) ! when f/=0 the we want to record the mirror
       n_mir = 0
       do j = 1,size(c,2)
          do i = 1,size(c,1)
             if(f(i,j)/=0)then
                n_mir = n_mir+1 ! count points on the mirror
             endif
          enddo
       enddo
       allocate(i_mir(n_mir)) ! allocate space for mirror
       allocate(j_mir(n_mir))
       allocate(f_mir(nt,n_mir))
       n_mir = 0
       do j = 1,size(c,2)
          do i = 1,size(c,1)
             if(f(i,j)/=0)then
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
    call pgwnad(0.5,size(c,1)+0.5,0.5,size(c,2)+0.5)
    !
    p_old(:,:) = 0.D0 ! init field variables
    p(:,:)     = 0.D0
    p_new(:,:) = 0.D0
    f(:,:)     = 0.D0
    !
    do it = 1,nt ! time loop
       !
       ! Record wavefield
       do ir = 1,size(i_rec)
          p_rec(it,ir) = p(i_rec(ir),j_rec(ir))
       enddo
       !
       ! compute "internal forces" (i.e. right hand side in eq.55; Masson et al., 2018)
       f = force_int(p,c,dx,dy)

       ! compute/record mirror force if needed (eq.55 in Masson et al., 2018)
       if(present(win_mir))then
          fm = -force_int(p*win_mir,c,dx,dy)+win_mir*force_int(p,c,dx,dy)
          do ir = 1,n_mir
             f_mir(it,ir) = fm(i_mir(ir),j_mir(ir))
          enddo
       endif
       !
       ! add sources
       do is = 1,size(i_src)
          f(i_src(is),j_src(is)) = f(i_src(is),j_src(is))+stf_src(it,is)
       enddo
       !
       p_new = 2*p-p_old+dt*dt*f ! time integration
       !
       p_old = p ! update the time variables
       p = p_new
       !
       ! below is purely graphics
       !
       ! plot wavefield
       call pgimag(sngl(p),size(p,1),size(p,2),1,size(p,1),1,size(p,2), &
            -maxval(sngl(abs(p))),maxval(sngl(abs(p))),(/0.,1.,0.,0.,0.,1./))
       !
       ! plot mirror window if present
       if(present(win_mir))then
          do i = 1,size(c,1)-1
             do j = 1,size(c,2)
                if(win_mir(i,j)/=win_mir(i+1,j))then
                   call pgmove(real(i)+0.5,real(j)-0.5)
                   call pgdraw(real(i)+0.5,real(j+1)-0.5)
                endif
             enddo
          enddo
          do i = 1,size(c,1)
             do j = 1,size(c,2)-1
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
  !----------------------------------
  end subroutine solver_FD2D_acoustic
  !----------------------------------

  !----------------------------
  function force_int(p,c,dx,dy)
    !----------------------------
    !
    ! compute the right hand side of eq.54 in Masson et al., 2014 (without the source term) 
    !
    implicit none
    !
    double precision :: p(:,:),c(:,:),dx,dy
    double precision :: force_int(size(p,1),size(p,2))
    !
    integer :: i,j
    !
    force_int(:,:) = 0.D0
    !
    do j = 2,size(p,2)-1
       do i = 2,size(p,1)-1
          force_int(i,j) = c(i,j)**2*( (p(i+1,j)-2*p(i,j)+p(i-1,j))/dx/dx &
               +(p(i,j+1)-2*p(i,j)+p(i,j-1))/dy/dy)
       enddo
    enddo
    !
    return
    !
    !---------------------
  end function force_int
  !---------------------

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


