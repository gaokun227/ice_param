
module ice_albedo_mod

!=======================================================================
!
!                     ICE SURFACE ALBEDO MODULE
!
!   Routines for computing the surface albedo over ice 
!
!=======================================================================

use      utilities_mod, only:  error_mesg, file_exist,  &
                               check_nml_error, open_file,  &
                               FATAL, get_my_pe, close_file

implicit none
private

!======= public interface =============================================

public  ice_albedo, ice_albedo_init

!--------------------- version number ----------------------------------

character(len=128) :: version = '$Id: ice_albedo.F90,v 1.2 2000/07/28 20:17:11 fms Exp $'
character(len=128) :: tag = '$Name: galway $'

!=======================================================================

!     DEFAULT VALUES OF NAMELIST PARAMETERS:

real :: crit_thickness       = 1.00   
real :: t_range              = 10.0

real :: min_ice_alb          = 0.50
real :: max_ice_alb          = 0.80

real :: const_alb            = 0.65


namelist /ice_albedo_nml/  crit_thickness      , &
                           t_range             , &
                           min_ice_alb         , &
                           max_ice_alb         , &
                           const_alb

!=======================================================================

!  OTHER MODULE VARIABLES

real :: temp_ice_freeze

logical :: do_init = .true.

CONTAINS

!#######################################################################

subroutine ice_albedo_init ( t_freeze )

real, intent(in) :: t_freeze

integer  unit, io, ierr

temp_ice_freeze = t_freeze

!------------------- read namelist input -------------------------------

      if (file_exist('input.nml')) then
         unit = open_file ('input.nml', action='read')
         ierr=1; do while (ierr /= 0)
            read  (unit, nml=ice_albedo_nml, iostat=io, end=10)
            ierr = check_nml_error(io,'ice_albedo_nml')
         enddo
  10     call close_file (unit)
      endif

!---------- output namelist to log-------------------------------------

      unit = open_file ('logfile.out', action='append')
      if ( get_my_pe() == 0 ) then
           write (unit, '(/,80("="),/(a))') trim(version),trim(tag)
           write (unit, nml=ice_albedo_nml)
      endif
      call close_file (unit)

  do_init = .false.

end subroutine ice_albedo_init

!#######################################################################

subroutine ice_albedo (ice, thickness, temp, albedo, i1, j1)

!-----------------------------------------------------------------------
!----------- CALCULATE SURFACE ALBEDO OVER ICE ------------------------
!-----------------------------------------------------------------------
!
! INPUT
!     (i1,j1)    = postition in global ice grid of first elements in
!                       the input
!     thickness  = thickness of ice (in meters)
!     alb_ocean  - albedo of open ocean (on ice grid)
!     temp       = surface temperature (in degrees kelvin)
!
!  OUTPUT
!     albedo = surface albedo

logical, intent(in),    dimension(:,:,:) :: ice
real,    intent(in),    dimension(:,:,:) :: temp, thickness
real,    intent(inout), dimension(:,:,:) :: albedo
integer, intent(in), optional :: i1, j1

!-----------------------------------------------------------------------

integer :: is, js, ie, je, n
real :: tcrit
real, dimension(size(temp,1),size(temp,2)) :: thick_ice_alb, alb_ocean

if (do_init) call error_mesg ('ice_albedo',  &
                              'initialization not called', FATAL)

is = 1; if (present(i1)) is = i1
js = 1; if (present(j1)) js = j1

ie = is + size(temp,1) - 1
je = js + size(temp,2) - 1

tcrit = temp_ice_freeze - t_range

!--- first partition assumed to be open ocean ----

 alb_ocean = albedo(:,:,1)

 do n = 2, size(temp,3)

  where(ice(:,:,n) .and. temp(:,:,n) <= tcrit) &
             thick_ice_alb = max_ice_alb
  where(ice(:,:,n) .and. temp(:,:,n) >= temp_ice_freeze) &
             thick_ice_alb = min_ice_alb
  where(ice(:,:,n) .and. temp(:,:,n) <  temp_ice_freeze &
                               .and. temp(:,:,n) > tcrit)  &
     thick_ice_alb = max_ice_alb +  &
           (min_ice_alb - max_ice_alb)*(temp(:,:,n) - tcrit)/t_range


  where (ice(:,:,n).and.thickness(:,:,n) >= crit_thickness)  &
        albedo(:,:,n) = thick_ice_alb
  where (ice(:,:,n).and.thickness(:,:,n) < crit_thickness) &
        albedo(:,:,n) =  alb_ocean + &
      (thick_ice_alb - alb_ocean) *sqrt(thickness(:,:,n)/crit_thickness)

end do

end subroutine ice_albedo

!######################################################################

end module ice_albedo_mod

