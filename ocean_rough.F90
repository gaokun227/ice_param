!***********************************************************************
!*                   GNU Lesser General Public License
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* FMS is free software: you can redistribute it and/or modify it under
!* the terms of the GNU Lesser General Public License as published by
!* the Free Software Foundation, either version 3 of the License, or (at
!* your option) any later version.
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
!* FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
!* for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with FMS.  If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************

module ocean_rough_mod

!-----------------------------------------------------------------------

use          mpp_mod, only: input_nml_file

use       fms_mod, only: error_mesg, FATAL,  mpp_error, &
                         check_nml_error, mpp_pe, mpp_root_pe, &
                         write_version_number, stdlog
use constants_mod, only: grav, vonkarm

implicit none
private

public :: compute_ocean_roughness, fixed_ocean_roughness
public :: cal_z0_hwrf17, cal_zt_hwrf17, cal_z0_hwrf15, cal_zt_hwrf15

!-----------------------------------------------------------------------
character(len=256) :: version = '$Id$'
character(len=256) :: tagname = '$Name$'
!-----------------------------------------------------------------------
!----- namelist -----

  real    :: roughness_init = 0.00044   ! not used in this version
  real    :: roughness_min  = 1.e-6
  real    :: charnock       = 0.032
  
  real    :: roughness_mom   = 5.8e-5
  real    :: roughness_heat  = 5.8e-5   ! was 4.00e-4
  real    :: roughness_moist = 5.8e-5
!  real, parameter :: zcoh1 = 0.0       ! Beljaars 1994 values
!  real, parameter :: zcoq1 = 0.0
! real, parameter :: zcoh1 = 1.4e-5
! real, parameter :: zcoq1 = 1.3e-4
  real            :: zcoh1 = 0.0 !miz
  real            :: zcoq1 = 0.0 !miz
  real            :: v10m  = 32.5 !jhc
  real            :: v10n  = 17.5 !jhc
  logical :: do_highwind     = .false.
  logical :: do_cap40        = .false.

  character(len=32) :: rough_scheme = 'fixed'   ! possible values:
                                                !   'fixed'
                                                !   'charnock'
                                                !   'beljaars'

namelist /ocean_rough_nml/ roughness_init, roughness_heat,  &
                           roughness_mom,  roughness_moist, &
                           roughness_min,                   &
                           charnock,                        &
                           rough_scheme, do_highwind,       &!miz
                           do_cap40, zcoh1, zcoq1,          &!sjl
                           v10m, v10n                        !jhc


!-----------------------------------------------------------------------

  logical :: do_init = .true.

!-----------------------------------------------------------------------
! ---- constants ----

! ..... high wind speed - rough sea
  real, parameter :: zcom1 = 1.8e-2    ! Charnock's constant
! ..... low wind speed - smooth sea
  real, parameter :: gnu   = 1.5e-5
  real, parameter :: zcom2 = 0.11
  real, parameter :: zcoh2 = 0.40
  real, parameter :: zcoq2 = 0.62


contains

!#######################################################################

 subroutine compute_ocean_roughness ( ocean, u_star,  &
                                      rough_mom, rough_heat, rough_moist )

 logical, intent(in)  :: ocean(:,:)
 real,    intent(in)  :: u_star(:,:)
 real,    intent(out) :: rough_mom(:,:), rough_heat(:,:), rough_moist(:,:)

!-----------------------------------------------------------------------
!  computes ocean roughness for momentum using wind stress
!  and sets roughness for heat/moisture using namelist value
!-----------------------------------------------------------------------

   real, dimension(size(ocean,1),size(ocean,2)) :: ustar2, xx1, xx2, w10 !miz
   real, dimension(size(ocean,1),size(ocean,2)) :: ustar, xx3, u10n, z0, zt, z1, alpha_v, reynolds_rough
   real:: zt1
   integer :: i, j, n, iter
   real :: ustar_min, m, b, u_max, rough_mom_init

   if (do_init) call ocean_rough_init


   if (trim(rough_scheme) == 'fixed') then

!  --- set roughness for momentum and heat/moisture ---

      call fixed_ocean_roughness ( ocean, rough_mom, rough_heat, &
                                          rough_moist )


!  --- compute roughness for momentum, heat, moisture ---

   else if (trim(rough_scheme) == 'beljaars' .or. &
            trim(rough_scheme) == 'charnock') then

      where (ocean)
          ustar2(:,:) = max(gnu*gnu,u_star(:,:)*u_star(:,:))          
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
! --- SJL ---- High Wind correction following Moon et al 2007 ------
          if (do_highwind) then       !  Moon et al. formular
              do j=1,size(ocean,2)
                 do i=1,size(ocean,1)
                    if ( ocean(i,j) ) then
                      w10(i,j) = 2.458 + u_star(i,j)*(20.255-0.56*u_star(i,j))  ! Eq(7) Moon et al.
                      if ( w10(i,j) > 12.5 ) then
                           rough_mom(i,j) = 0.001*(0.085*w10(i,j) - 0.58)    ! Eq(8b) Moon et al.
! SJL mods: cap the growth of z0 with w10 up to 40 m/s
! z0 (w10=40) = 2.82E-3
                           if(do_cap40) rough_mom(i,j) = min( rough_mom(i,j), 2.82E-3)
                      else
                           rough_mom(i,j) = 0.0185/grav*u_star(i,j)**2  ! (8a) Moon et al.
                      endif
                           zt1 = min( 1., max(0., (w10(i,j)-v10n)/(v10m-v10n)) )
                           rough_moist(i,j) = zcoq1*zt1*xx2(i,j) + zcoq2 * xx1(i,j)
                           rough_heat (i,j) = zcoh1*zt1*xx2(i,j) + zcoh2 * xx1(i,j)
!                 --- lower limit on roughness? ---
                      rough_mom  (i,j) = max( rough_mom  (i,j), roughness_min )
                      rough_heat (i,j) = max( rough_heat (i,j), roughness_min )
                      rough_moist(i,j) = max( rough_moist(i,j), roughness_min )
                    endif
                 enddo
              enddo
! SJL -----------------------------------------------------------------------------------
          else
!     --- Beljaars scheme ---
          where (ocean)
              rough_mom  (:,:) = zcom1 * xx2(:,:) + zcom2 * xx1(:,:)
              rough_heat (:,:) = zcoh1 * xx2(:,:) + zcoh2 * xx1(:,:)
              rough_moist(:,:) = zcoq1 * xx2(:,:) + zcoq2 * xx1(:,:)
!             --- lower limit on roughness? ---
              rough_mom  (:,:) = max( rough_mom  (:,:), roughness_min )
              rough_heat (:,:) = max( rough_heat (:,:), roughness_min )
              rough_moist(:,:) = max( rough_moist(:,:), roughness_min )
          endwhere
          endif
      endif

   else if (trim(rough_scheme) == 'coare3.5') then    !IH -ref Edson et al, JPO, Aug 2013 

      !if rough_mom is available in input (rough_mon_in) then set rough_mom1_init = rough_mom_in
      rough_mom_init = 1.e-03
      ustar_min = 1.e-05 ! added by IH just in case ustar is unphysically small 
      m = 0.0017         !Edson 2013
      b = -0.005         !Edson 2013
      u_max = 18.0       !Edson 2013
      ustar(:,:)=u_star(:,:)
      where (ocean)
          ustar(:,:)  = max(ustar(:,:), ustar_min)  ! IH
	  ustar2(:,:) = ustar(:,:)*ustar(:,:)
          xx1(:,:)    = gnu/ustar(:,:)
          xx2(:,:)    = ustar2(:,:)/grav
	  xx3(:,:)    = ustar(:,:)/vonkarm 
      endwhere	  
      !if rough_mom is available in input then set z0(:,:) to this input value
      z0(:,:) = rough_mom_init
      iter = 5      ! IH - overkill
      do j=1,size(ocean,2)
     	do i=1,size(ocean,1)
           if ( ocean(i,j) ) then
              do n = 1, iter
	      	 u10n(i,j) = xx3(i,j)*log(10/z0(i,j))             ! "neutral" 10m wind
	    	 alpha_v(i,j) = m*min(u10n(i,j),u_max) + b;     ! Charnock coefficient, Edson 2013
	    	 z1(i,j) = zcom2*xx1(i,j) + alpha_v(i,j)*xx2(i,j) ! Edson 2013
	    	 z0(i,j) = z1(i,j);                               !IH -- iteration
	      enddo
           endif
        enddo
      enddo

      where (ocean)
      	  rough_mom  (:,:) = z0(:,:)
      	  rough_mom  (:,:) = max( rough_mom  (:,:), roughness_min )

	  reynolds_rough(:,:) = ustar(:,:)*rough_mom(:,:)/gnu
	  rough_heat (:,:)    = 5.5e-05*(reynolds_rough(:,:)**(-0.6))
	  rough_heat (:,:)    = min(1.1e-04, rough_heat(:,:))
	  rough_moist(:,:)    = rough_heat (:,:)
      elsewhere
          rough_mom   = 0.0
          rough_heat  = 0.0
          rough_moist = 0.0
      endwhere

   ! KGao - introduce z0/zt options as in HWRF 2015 and 2017 versions
   ! this may cross over Baoqiang's work
   ! the whole point is to obtain u10n, z0, zt based on ustar via iterations
   else if (trim(rough_scheme) == 'hwrf15' .or. trim(rough_scheme) == 'hwrf17') then

      rough_mom_init = 1.e-03
      ustar_min = 1.e-05
      ustar(:,:) = u_star(:,:)
      where (ocean)
          ustar(:,:)  = max(ustar(:,:), ustar_min)
          ustar2(:,:) = ustar(:,:)*ustar(:,:)
          xx1(:,:)    = gnu/ustar(:,:)
          xx2(:,:)    = ustar2(:,:)/grav
          xx3(:,:)    = ustar(:,:)/vonkarm
      endwhere

      z0(:,:) = rough_mom_init
      iter = 5

      do j = 1, size(ocean, 2)
        do i = 1, size(ocean, 1)
           if ( ocean(i,j) ) then
              !z0(i,j) = 0.005 * xx2(i,j) ! KGao test - use charnock coeff to get initial guess for z0
              do n = 1, iter
                 u10n(i,j) = xx3(i,j) * log(10 / z0(i,j))
                 if (trim(rough_scheme) == 'hwrf15') then
                   call cal_z0_hwrf15(u10n(i,j), z0(i,j))
                 else if (trim(rough_scheme) == 'hwrf17') then
                   call cal_z0_hwrf17(u10n(i,j), z0(i,j))
                 endif
              enddo
              if (trim(rough_scheme) == 'hwrf15') then
                 call cal_zt_hwrf15(u10n(i,j), zt(i,j))
              else if (trim(rough_scheme) == 'hwrf17') then
                 call cal_zt_hwrf17(u10n(i,j), zt(i,j))
              endif
           endif
        enddo
      enddo

      where (ocean)
          rough_mom  (:,:) = z0(:,:)
          rough_mom  (:,:) = max( rough_mom  (:,:), roughness_min )
          rough_heat (:,:) = min(1.1e-04, zt(:,:))
          rough_moist(:,:) = rough_heat (:,:)
      elsewhere
          rough_mom   = 0.0
          rough_heat  = 0.0
          rough_moist = 0.0
      endwhere

   else
      call mpp_error(FATAL, '==>Error from ocean_rough_mod(compute_ocean_roughness): '//&
            'Unknown roughness scheme (case sensitive): ' //trim(rough_scheme))
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
    endwhere

 end subroutine fixed_ocean_roughness

!#######################################################################

 subroutine ocean_rough_init

   integer :: unit, ierr, io

!   ----- read and write namelist -----

  read (input_nml_file, nml=ocean_rough_nml, iostat=io)
  ierr = check_nml_error(io, 'ocean_rough_nml')

!------- write version number and namelist ---------

    if ( mpp_pe() == mpp_root_pe() ) then
         call write_version_number(version, tagname)
         unit = stdlog()
         write (unit,nml=ocean_rough_nml)
         write (unit,11)
    endif

!------ constants -----

    roughness_moist = max (roughness_moist, roughness_min)
    roughness_heat  = max (roughness_heat , roughness_min)
    roughness_mom   = max (roughness_mom  , roughness_min)

    do_init = .false.

11 format (/,'namelist option USE_FIXED_ROUGH is no longer supported', &
           /,'use variable ROUGH_SCHEME instead')

 end subroutine ocean_rough_init

!#######################################################################
      
 subroutine cal_z0_hwrf15(ws10m, z0)
      ! coded by Kun Gao (Kun.Gao@noaa.gov)
      ! originally developed by URI/GFDL
      real, intent (in) :: ws10m
      real, intent (out):: z0
      real, parameter ::         &
      a0=-8.367276172397277e-12, &          
      a1=1.7398510865876079e-09, &
      a2=-1.331896578363359e-07, &
      a3=4.507055294438727e-06,  &
      a4=-6.508676881906914e-05, &
      a5=0.00044745137674732834, &
      a6=-0.0010745704660847233, &
      b0=2.1151080765239772e-13, &
      b1=-3.2260663894433345e-11,&
      b2=-3.329705958751961e-10, &
      b3=1.7648562021709124e-07, &
      b4=7.107636825694182e-06,  &
      b5=-0.0013914681964973246, &
      b6=0.0406766967657759
               
      if (ws10m <= 5.0) then
         z0 = 0.0185/9.8*(7.59e-4*ws10m**2+2.46e-2*ws10m)**2
      elseif (ws10m > 5.0 .and. ws10m <= 10.) then
         z0 = 0.00000235*(ws10m**2-25.)+3.805129199617346e-05
      elseif (ws10m > 10.0 .and. ws10m <= 60.) then
         z0 = a6 + a5*ws10m + a4*ws10m**2 + a3*ws10m**3&
            + a2*ws10m**4 + a1*ws10m**5 + a0*ws10m**6
      else
         z0 = b6 + b5*ws10m + b4*ws10m**2 + b3*ws10m**3&
            + b2*ws10m**4 + b1*ws10m**5 + b0*ws10m**6
      endif

 end subroutine cal_z0_hwrf15

 subroutine cal_zt_hwrf15(ws10m, zt)
      ! coded by Kun Gao (Kun.Gao@noaa.gov)
      ! originally developed by URI/GFDL
      real, intent (in) :: ws10m
      real, intent (out):: zt
      real, parameter ::      &
      a0 = 2.51715926619e-09, &
      a1 = -1.66917514012e-07,&
      a2 = 4.57345863551e-06, &
      a3 = -6.64883696932e-05,&
      a4 = 0.00054390175125,  &
      a5 = -0.00239645231325, &
      a6 = 0.00453024927761,  &
      b0 = -1.72935914649e-14,&
      b1 = 2.50587455802e-12, &
      b2 = -7.90109676541e-11,&
      b3 = -4.40976353607e-09,&
      b4 = 3.68968179733e-07, &
      b5 = -9.43728336756e-06,&
      b6 = 8.90731312383e-05, &
      c0 = 4.68042680888e-14, &
      c1 = -1.98125754931e-11,&
      c2 = 3.41357133496e-09, &
      c3 = -3.05130605309e-07,&
      c4 = 1.48243563819e-05, &
      c5 = -0.000367207751936,&
      c6 = 0.00357204479347

      if (ws10m <= 7.0) then
         zt = 0.0185/9.8*(7.59e-4*ws10m**2+2.46e-2*ws10m)**2
      elseif (ws10m > 7.0 .and. ws10m <= 15.) then
         zt = a6 + a5*ws10m + a4*ws10m**2 + a3*ws10m**3&
            + a2*ws10m**4 + a1*ws10m**5 + a0*ws10m**6
      elseif (ws10m > 15.0 .and. ws10m <= 60.) then
         zt = b6 + b5*ws10m + b4*ws10m**2 + b3*ws10m**3&
            + b2*ws10m**4 + b1*ws10m**5 + b0*ws10m**6
      else
         zt = c6 + c5*ws10m + c4*ws10m**2 + c3*ws10m**3&
            + c2*ws10m**4 + c1*ws10m**5 + c0*ws10m**6
      endif
 end subroutine cal_zt_hwrf15

 subroutine cal_z0_hwrf17(ws10m, z0)
      ! coded by Kun Gao (Kun.Gao@noaa.gov)
      real, intent (in) :: ws10m
      real, intent (out):: z0
      real, parameter ::          &
      p13=-1.296521881682694e-02, &
      p12= 2.855780863283819e-01, &
      p11=-1.597898515251717e+00, &
      p10=-8.396975715683501e+00, &
      p25= 3.790846746036765e-10, &
      p24= 3.281964357650687e-09, &
      p23= 1.962282433562894e-07, &
      p22=-1.240239171056262e-06, &
      p21=1.739759082358234e-07,  &
      p20=2.147264020369413e-05,  &
      p35=1.840430200185075e-07,  &
      p34=-2.793849676757154e-05, &
      p33=1.735308193700643e-03,  &
      p32=-6.139315534216305e-02, &
      p31=1.255457892775006e+00,  &
      p30=-1.663993561652530e+01, &
      p40=4.579369142033410e-04

      if (ws10m <= 6.5) then
         z0 = exp( p10 + p11*ws10m + p12*ws10m**2 + p13*ws10m**3)
      elseif (ws10m > 6.5 .and. ws10m <= 15.7) then
         z0 = p25*ws10m**5 + p24*ws10m**4 + p23*ws10m**3&
            + p22*ws10m**2 + p21*ws10m + p20
      elseif (ws10m > 15.7 .and. ws10m <= 53.) then
         z0 = exp( p35*ws10m**5 + p34*ws10m**4 + p33*ws10m**3&
            + p32*ws10m**2 + p31*ws10m + p30 )
      else
         z0 = p40
      endif
      !if (ws10m > 30.) z0 = max(z0,0.003)

 end subroutine cal_z0_hwrf17

 subroutine cal_zt_hwrf17(ws10m, zt)
      ! coded by Kun Gao (Kun.Gao@noaa.gov)
      real, intent (in) :: ws10m
      real, intent (out):: zt
      real, parameter  :: p00 =  1.100000000000000e-04, &
           p15 = -9.144581627678278e-10, p14 =  7.020346616456421e-08,&
           p13 = -2.155602086883837e-06, p12 =  3.333848806567684e-05,&
           p11 = -2.628501274963990e-04, p10 =  8.634221567969181e-04,&
           p25 = -8.654513012535990e-12, p24 =  1.232380050058077e-09,&
           p23 = -6.837922749505057e-08, p22 =  1.871407733439947e-06,&
           p21 = -2.552246987137160e-05, p20 =  1.428968311457630e-04,&
           p35 =  3.207515102100162e-12, p34 = -2.945761895342535e-10,&
           p33 =  8.788972147364181e-09, p32 = -3.814457439412957e-08,&
           p31 = -2.448983648874671e-06, p30 =  3.436721779020359e-05,&
           p45 = -3.530687797132211e-11, p44 =  3.939867958963747e-09,&
           p43 = -1.227668406985956e-08, p42 = -1.367469811838390e-05,&
           p41 =  5.988240863928883e-04, p40 = -7.746288511324971e-03,&
           p56 = -1.187982453329086e-13, p55 =  4.801984186231693e-11,&
           p54 = -8.049200462388188e-09, p53 =  7.169872601310186e-07,&
           p52 = -3.581694433758150e-05, p51 =  9.503919224192534e-04,&
           p50 = -1.036679430885215e-02,&
           p60 =  4.751256171799112e-05

      if (ws10m >= 0.0 .and. ws10m < 5.9 ) then
         zt = p00
      elseif (ws10m >= 5.9 .and. ws10m <= 15.4) then
         zt = p10 + ws10m * (p11 + ws10m * (p12 + ws10m * (p13&
                  + ws10m * (p14 + ws10m * p15))))
      elseif (ws10m > 15.4 .and. ws10m <= 21.6) then
         zt = p20 + ws10m * (p21 + ws10m * (p22 + ws10m * (p23&
                  + ws10m * (p24 + ws10m * p25))))
      elseif (ws10m > 21.6 .and. ws10m <= 42.2) then
         zt = p30 + ws10m * (p31 + ws10m * (p32 + ws10m * (p33&
                  + ws10m * (p34 + ws10m * p35))))
      elseif ( ws10m > 42.2 .and. ws10m <= 53.3) then
         zt = p40 + ws10m * (p41 + ws10m * (p42 + ws10m * (p43&
                  + ws10m * (p44 + ws10m * p45))))
      elseif ( ws10m > 53.3 .and. ws10m <= 80.0) then
         zt = p50 + ws10m * (p51 + ws10m * (p52 + ws10m * (p53&
                  + ws10m * (p54 + ws10m * (p55 + ws10m * p56)))))
      elseif ( ws10m > 80.0) then
         zt = p60
      endif
 end subroutine cal_zt_hwrf17

end module ocean_rough_mod

