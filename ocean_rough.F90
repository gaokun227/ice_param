
module ocean_rough_mod

!-----------------------------------------------------------------------

use utilities_mod, only: error_mesg, FATAL, file_exist, open_file,  &
                         check_nml_error, get_my_pe, close_file
use constants_mod, only: grav, vonkarm

implicit none
private

public :: compute_ocean_roughness, fixed_ocean_roughness

!-----------------------------------------------------------------------
character(len=256) :: version = '$Id: ocean_rough.F90,v 1.2 2000/07/28 20:37:24 fms Exp $'
character(len=256) :: tag = '$Name: bombay $'
!-----------------------------------------------------------------------
!----- namelist -----

  real    :: charnock       = 0.32
  real    :: roughness_init = 0.00044   ! not used in this version
  real    :: roughness_min  = 1.e-6
  real    :: rho_atm        = 1.13

  logical :: use_fixed_rough = .true.
  real    :: roughness_mom   = 5.8e-5
  real    :: roughness_heat  = 5.8e-5   ! was 4.00e-4

namelist /ocean_rough_nml/ roughness_init, roughness_heat,  &
                           roughness_mom,  roughness_min,   &
                           use_fixed_rough, charnock, rho_atm

!-----------------------------------------------------------------------

  real    :: zconst
  logical :: do_init = .true.

  real, parameter :: z_ref  = 10.
  real, parameter :: cd_ref = 1.1e-3
  real, parameter :: ksq    = vonkarm*vonkarm

contains

!#######################################################################

 subroutine compute_ocean_roughness ( ocean, tau_x, tau_y,  &
                                      rough_mom, rough_heat )

 logical, intent(in)  :: ocean(:,:)
 real,    intent(in)  :: tau_x(:,:), tau_y(:,:)
 real,    intent(out) :: rough_mom(:,:), rough_heat(:,:)

!-----------------------------------------------------------------------
!  computes ocean roughness for momentum using wind stress
!  and sets roughness for heat/moisture using namelist value
!-----------------------------------------------------------------------

   if (do_init) call ocean_rough_init


   if (use_fixed_rough) then

!  --- set roughness for momentum and heat/moisture ---

      call fixed_ocean_roughness ( ocean, rough_mom, rough_heat )

   else

!  --- compute roughness for momentum ---
!  --- back-compute roughness for heat/moisture -----
!          assuming neutral conditions

      where (ocean)
          rough_mom  = zconst * sqrt(tau_x*tau_x + tau_y*tau_y)
          rough_mom  = max (rough_mom , roughness_min)
!         --- use momentum roughness for heat ---
          rough_heat = rough_mom
!!!!      rough_heat = z_ref * exp(-ksq / (cd_ref*log(z_ref/rough_mom)))
!!!!      rough_heat = max (rough_heat, roughness_min)
      elsewhere
          rough_mom  = 0.0
          rough_heat = 0.0
      endwhere

   endif

!-----------------------------------------------------------------------

 end subroutine compute_ocean_roughness

!#######################################################################

 subroutine fixed_ocean_roughness ( ocean, rough_mom, rough_heat )

 logical, intent(in)  :: ocean(:,:)
 real,    intent(out) :: rough_mom(:,:), rough_heat(:,:)

   if (do_init) call ocean_rough_init

    where (ocean)
       rough_mom  = roughness_mom
       rough_heat = roughness_heat
    elsewhere
       rough_mom  = 0.0
       rough_heat = 0.0
    endwhere

 end subroutine fixed_ocean_roughness

!#######################################################################

 subroutine ocean_rough_init

   integer :: unit, ierr, io

!   ----- read and write namelist -----

    if ( file_exist('input.nml')) then
        unit = open_file (file='input.nml', action='read')
        ierr=1; do while (ierr /= 0)
           read  (unit, nml=ocean_rough_nml, iostat=io, end=10)
           ierr = check_nml_error(io,'ocean_rough_nml')
        enddo
 10     call close_file (unit)
    endif

    unit = open_file ('logfile.out', action='append')
    if ( get_my_pe() == 0 ) then
         write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
         write (unit, nml=ocean_rough_nml)
    endif
    call close_file (unit)

!------ constants -----

    zconst         = charnock / (grav * rho_atm)
    roughness_heat = max (roughness_heat, roughness_min)
    roughness_mom  = max (roughness_mom , roughness_min)

    do_init = .false.

 end subroutine ocean_rough_init

!#######################################################################

end module ocean_rough_mod

