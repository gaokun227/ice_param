
module ocean_rough_mod

!-----------------------------------------------------------------------

use utilities_mod, only: error_mesg, FATAL, file_exist, open_file,  &
                         check_nml_error, get_my_pe, close_file
use constants_mod, only: grav

implicit none
private

public :: compute_ocean_roughness, fixed_ocean_roughness

!-----------------------------------------------------------------------
character(len=256) :: version = '$Id: ocean_rough.F90,v 1.4 2002/07/16 22:48:06 fms Exp $'
character(len=256) :: tag = '$Name: inchon $'
!-----------------------------------------------------------------------
!----- namelist -----

  real    :: roughness_init = 0.00044   ! not used in this version
  real    :: roughness_min  = 1.e-6
  real    :: charnock       = 0.032
  real    :: rho_atm        = 1.13

  real    :: roughness_mom   = 5.8e-5
  real    :: roughness_heat  = 5.8e-5   ! was 4.00e-4
  real    :: roughness_moist = 5.8e-5

  character(len=32) :: rough_scheme = 'fixed'   ! possible values:
                                                !   'fixed'
                                                !   'charnock'
                                                !   'beljaars'

namelist /ocean_rough_nml/ roughness_init, roughness_heat,  &
                           roughness_mom,  roughness_moist, &
                           roughness_min,  rho_atm,         &
                           charnock,                        &
                           rough_scheme

!-----------------------------------------------------------------------

  logical :: do_init = .true.

!-----------------------------------------------------------------------
! ---- constants ----

! ..... high wind speed - rough sea
  real, parameter :: zcom1 = 1.8e-2    ! Charnock's constant
  real, parameter :: zcoh1 = 0.0       ! Beljaars 1994 values
  real, parameter :: zcoq1 = 0.0
! real, parameter :: zcoh1 = 1.4e-5
! real, parameter :: zcoq1 = 1.3e-4
! ..... low wind speed - smooth sea
  real, parameter :: gnu   = 1.5e-5
  real, parameter :: zcom2 = 0.11
  real, parameter :: zcoh2 = 0.40
  real, parameter :: zcoq2 = 0.62


contains

!#######################################################################

 subroutine compute_ocean_roughness ( ocean, tau_x, tau_y,  &
                                      rough_mom, rough_heat, rough_moist )

 logical, intent(in)  :: ocean(:,:)
 real,    intent(in)  :: tau_x(:,:), tau_y(:,:)
 real,    intent(out) :: rough_mom(:,:), rough_heat(:,:), rough_moist(:,:)

!-----------------------------------------------------------------------
!  computes ocean roughness for momentum using wind stress
!  and sets roughness for heat/moisture using namelist value
!-----------------------------------------------------------------------

   real, dimension(size(ocean,1),size(ocean,2)) :: ustar2, xx1, xx2

   if (do_init) call ocean_rough_init


   if (trim(rough_scheme) == 'fixed') then

!  --- set roughness for momentum and heat/moisture ---

      call fixed_ocean_roughness ( ocean, rough_mom, rough_heat, &
                                          rough_moist )


!  --- compute roughness for momentum, heat, moisture ---

   else if (trim(rough_scheme) == 'beljaars' .or. &
            trim(rough_scheme) == 'charnock') then

      where (ocean)
          ustar2(:,:) = sqrt(tau_x*tau_x + tau_y*tau_y) / rho_atm
          ustar2(:,:) = max(gnu*gnu,ustar2(:,:))
          xx1(:,:) = gnu / sqrt(ustar2(:,:))
          xx2(:,:) = ustar2(:,:) / grav
      elsewhere
          rough_mom   = 0.0
          rough_heat  = 0.0
          rough_moist = 0.0
      endwhere

      if (trim(rough_scheme) == 'charnock') then
          where (ocean)
              rough_mom  (:,:) = charnock * xx2(:,:)
              rough_mom  (:,:) = max( rough_mom(:,:), roughness_min )
              rough_heat (:,:) = rough_mom  (:,:)
              rough_moist(:,:) = rough_mom  (:,:)
          endwhere
      else if (trim(rough_scheme) == 'beljaars') then
!     --- Beljaars scheme ---
          where (ocean)
              rough_mom  (:,:) = zcom1 * xx2(:,:) + zcom2 * xx1(:,:)
              rough_heat (:,:) = zcoh1            + zcoh2 * xx1(:,:)
              rough_moist(:,:) = zcoq1            + zcoq2 * xx1(:,:)
!             --- lower limit on roughness? ---
              rough_mom  (:,:) = max( rough_mom  (:,:), roughness_min )
              rough_heat (:,:) = max( rough_heat (:,:), roughness_min )
              rough_moist(:,:) = max( rough_moist(:,:), roughness_min )
          endwhere
      endif

   endif

!-----------------------------------------------------------------------

 end subroutine compute_ocean_roughness

!#######################################################################

 subroutine fixed_ocean_roughness ( ocean, rough_mom, rough_heat, rough_moist )

 logical, intent(in)  :: ocean(:,:)
 real,    intent(out) :: rough_mom(:,:), rough_heat(:,:), rough_moist(:,:)

   if (do_init) call ocean_rough_init

    where (ocean)
       rough_mom   = roughness_mom
       rough_heat  = roughness_heat
       rough_moist = roughness_moist
    elsewhere
       rough_mom   = 0.0
       rough_heat  = 0.0
       rough_moist = 0.0
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
         write (unit,11)
    endif
    call close_file (unit)

!------ constants -----

    roughness_moist = max (roughness_moist, roughness_min)
    roughness_heat  = max (roughness_heat , roughness_min)
    roughness_mom   = max (roughness_mom  , roughness_min)

    do_init = .false.

11 format (/,'namelist option USE_FIXED_ROUGH is no longer supported', &
           /,'use variable ROUGH_SCHEME instead')

 end subroutine ocean_rough_init

!#######################################################################

end module ocean_rough_mod

