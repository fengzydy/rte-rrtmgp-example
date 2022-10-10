program stratospheric_adjustments

use, intrinsic :: iso_fortran_env, only: error_unit, output_unit, real32
use netcdf

use mo_fluxes, only: ty_fluxes_broadband
use mo_gas_concentrations, only: ty_gas_concs
use mo_gas_optics_rrtmgp, only: ty_gas_optics_rrtmgp
use mo_optical_props, only: ty_optical_props_1scl, ty_optical_props_2str
use mo_rfmip_io, only: determine_gas_names, read_and_block_gases_ty, read_and_block_lw_bc, &
                       read_and_block_pt, read_and_block_sw_bc, read_size
use mo_rte_kind, only: wp
use mo_rte_lw, only: rte_lw
use mo_rte_sw, only: rte_sw
use mo_source_functions, only: ty_source_func_lw
use mo_load_coefficients, only: load_and_init

implicit none


character(len=132) :: rfmip_file
character(len=64) :: name_of_program
character(len=64) :: buffer
integer :: nargs
integer :: ncol
integer :: nlay
integer :: nexp
integer :: forcing_index
integer :: physics_index
integer :: n_quad_angles
logical :: top_at_1
integer :: nbnd
integer :: ngpt
character(len=4) :: forcing_index_char
character(len=4) :: physics_index_char
real(kind=wp), dimension(:,:,:), allocatable :: p_lay !Atmospheric layer pressure [Pa] (ncol, nlay, nexp).
real(kind=wp), dimension(:,:,:), allocatable :: p_lev !Atmospheric level pressure [Pa] (ncol, nlev, nexp).
real(kind=wp), dimension(:,:,:), allocatable :: t_lay !Atmospheric layer temperature [K] (ncol, nlay, nexp).
real(kind=wp), dimension(:,:,:), allocatable :: t_lev !Atmoshperic level temperature [K] (ncol, nlev, nexp).
real(kind=wp), dimension(:,:), allocatable :: sfc_emis !Surface emissivity (ncol, nexp).
real(kind=wp), dimension(:,:), allocatable :: sfc_t !Surface temperature [K] (ncol, nexp).
real(kind=wp), dimension(:,:), allocatable :: sfc_emis_spec !Surface emissivity (nbnd, ncol).
type(ty_source_func_lw) :: source !Source function object.
real(kind=wp), dimension(:,:), allocatable :: surface_albedo !Surface albedo (ncol, nexp).
real(kind=wp), dimension(:,:), allocatable :: total_solar_irradiance !Solar irradiance [W/m^2] (ncol, nexp).
real(kind=wp), dimension(:,:), allocatable :: solar_zenith_angle !Solar zenith angle (ncol, nexp).
real(kind=wp), dimension(:,:), allocatable :: sfc_alb_spec !Surface albedo (nbnd, ncol).
real(kind=wp), dimension(:,:), allocatable :: toa_flux !Shortwave flux at top-of-atmosphere [W/m^2] (ncol, ngpt).
real(kind=wp), dimension(:), allocatable :: def_tsi !(ncol)
real(kind=wp), dimension(:), allocatable :: mu0 !Cosine of solar zenith angle (ncol).
logical, dimension(:,:), allocatable :: usecol !Is column daytime? (ncol, nexp).

character(len=132) :: lw_kdist_file !K-distribution configuration for the longwave.
character(len=32), dimension(:), allocatable :: lw_kdist_gas_names !Gas names used in the model.
character(len=32), dimension(:), allocatable :: lw_rfmip_gas_games !Gas names used in the input atmosphere file.
real(kind=wp), dimension(:,:,:), allocatable, target :: lw_flux_up !Upwelling longwave flux [W/m^2] (ncol, nlev, nexp).
real(kind=wp), dimension(:,:,:), allocatable, target :: lw_flux_dn !Downwelling longwave flux [W/m^2] (ncol, nlev, nexp).
type(ty_gas_optics_rrtmgp) :: lw_k_dist !Longwave k-distribution object.
type(ty_optical_props_1scl) :: lw_optical_props !Longwave optics object.
type(ty_fluxes_broadband) :: lw_fluxes !Longwave fluxes object.
type(ty_gas_concs), dimension(:), allocatable :: lw_gas_conc_array !Longwave gas concentrations object.

character(len=132) :: sw_kdist_file !K-distribution configuration for the shortwave.
character(len=32), dimension(:), allocatable :: sw_kdist_gas_names !Gas names used in the model.
character(len=32), dimension(:), allocatable :: sw_rfmip_gas_games !Gas names used in the input atmosphere file.
real(kind=wp), dimension(:,:,:), allocatable, target :: sw_flux_up !Upwelling shortwave flux [W/m^2] (ncol, nlev, nexp).
real(kind=wp), dimension(:,:,:), allocatable, target :: sw_flux_dn !Downwelling shortwave flux [W/m^2] (ncol, nlev, nexp).
type(ty_gas_optics_rrtmgp) :: sw_k_dist !Shortwave k-distribution object.
type(ty_optical_props_2str) :: sw_optical_props !Shortwave optics object.
type(ty_fluxes_broadband) :: sw_fluxes !Shortwave fluxes object.
type(ty_gas_concs), dimension(:), allocatable :: sw_gas_conc_array !Shortwave gas concentrations object.


real(kind=wp), parameter :: stratosphere_max_pressure = 20000._wp !Bottom of the stratosphere [Pa].
integer, dimension(:,:), allocatable :: stratosphere_starting_index !Index corresponding to the top, bottom of
                                                                    !the stratosphere if the column is oriented
                                                                    !toward, away from the surface (ncol, nexp).
integer, dimension(:,:), allocatable :: stratosphere_ending_index !Index corresponding to the bottom, top of
                                                                  !the stratosphere if the column is oriented
                                                                  !toward, away from the surface (ncol, nexp).
real(kind=wp), dimension(:,:,:), allocatable :: lw_heating_rate !Longwave heating rate [K/day] (ncol, nlay, nexp).
real(kind=wp), dimension(:,:,:), allocatable :: sw_heating_rate !Shortwave heating rate [K/day] (ncol, nlay, nexp).
real(kind=wp), dimension(:,:), allocatable :: dynamic_heating_rate !Dynamic heating rate [K/day] (ncol, nlay).
real(kind=wp) :: dq
real(kind=wp) :: d_hr
real(kind=wp) :: lw_hr
real(kind=wp) :: sw_hr
integer :: num_iterations = 5
integer :: b
integer :: i
integer :: j
integer :: k


enum, bind(c)
  enumerator :: expt = 1
  enumerator :: site
  enumerator :: layer
  enumerator :: level
  enumerator :: time
  enumerator :: num_dims = 5
end enum
enum, bind(c)
  enumerator :: rld = 1
  enumerator :: rlu
  enumerator :: rsd
  enumerator :: rsu
  enumerator :: rlhr
  enumerator :: rshr
  enumerator :: tlevid
  enumerator :: tlayid
  enumerator :: num_vars = 8
end enum
type :: OutputFile
  integer :: ncid
  integer, dimension(num_dims) :: dimid
  integer, dimension(num_vars) :: varid
end type
type(OutputFile), dimension(18) :: output


!Parse command line arguments.
call get_command_argument(0, name_of_program)
nargs = command_argument_count()
do i = 1, nargs
  call get_command_argument(i, buffer)
  if (trim(buffer) .eq. "-h" .or. trim(buffer) .eq. "--help") then
    write(output_unit, "(a)") "Usage: "//trim(name_of_program)//" [rfmip_file]" &
         //" [lw k-distribution_file] [sw k-distribution_file]" &
         //" [forcing_index (1,2,3)] [physics_index (1,2)]"
    stop
  endif
enddo

rfmip_file = "multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc"
if (nargs .ge. 1) call get_command_argument(1, rfmip_file)
call read_size(rfmip_file, ncol, nlay, nexp)

forcing_index_char = "1"
if (nargs .ge. 4) call get_command_argument(4, forcing_index_char)
read(forcing_index_char, '(i4)') forcing_index
if (forcing_index .lt. 1 .or. forcing_index .gt. 3) then
  call stop_on_err("Forcing index is invalid (must be 1,2 or 3).")
endif

physics_index_char = "1"
n_quad_angles = 1
if (nargs .ge. 5) call get_command_argument(5, physics_index_char)
read(physics_index_char, "(i4)") physics_index
if (physics_index .lt. 1 .or. physics_index .gt. 2) then
  call stop_on_err("Physics index is invalid (must be 1 or 2).")
endif
if (physics_index .eq. 2) n_quad_angles = 3


!Read input atmospheric data.
call read_and_block_pt(rfmip_file, ncol, p_lay, p_lev, t_lay, t_lev)
top_at_1 = p_lay(1,1,1) .lt. p_lay(1,nlay,1)
call read_and_block_lw_bc(rfmip_file, ncol, sfc_emis, sfc_t)
call read_and_block_sw_bc(rfmip_file, ncol, surface_albedo, &
                          total_solar_irradiance, solar_zenith_angle)


!Find pressure layer indices for the stratosphere.
allocate(stratosphere_starting_index(ncol, nexp))
allocate(stratosphere_ending_index(ncol, nexp))
if (top_at_1) then
  stratosphere_starting_index(:,:) = 1
else
  stratosphere_ending_index(:,:) = nlay
endif
do k = 1, nexp
  do i = 1, ncol
    do j = 1, nlay
      if (top_at_1) then
        if (p_lay(i,j,k) .gt. stratosphere_max_pressure) then
          !Leaving the stratosphere during downward sweep through the layers.
          stratosphere_ending_index(i,k) = j - 1
          exit
        endif
      else
        if (p_lay(i,j,k) .lt. stratosphere_max_pressure) then
          !Entering the stratosphere during upward sweep through the layers.
          stratosphere_starting_index(i,k) = j - 1
          exit
        endif
      endif
    enddo
  enddo
enddo


!Initialize longwave.
lw_kdist_file = "coefficients_lw.nc"
if (nargs .ge. 2) call get_command_argument(2, lw_kdist_file)
call determine_gas_names(rfmip_file, lw_kdist_file, forcing_index, &
                         lw_kdist_gas_names, lw_rfmip_gas_games)
print *, "Lw calculation uses RFMIP gases: ", &
         (trim(lw_rfmip_gas_games(b))//" ", b = 1, size(lw_rfmip_gas_games))
call read_and_block_gases_ty(rfmip_file, ncol, lw_kdist_gas_names, &
                             lw_rfmip_gas_games, lw_gas_conc_array)
call load_and_init(lw_k_dist, trim(lw_kdist_file), lw_gas_conc_array(1))
if (.not. lw_k_dist%source_is_internal()) then
  call stop_on_err("k-distribution file isn't LW.")
endif
if (top_at_1) then
  p_lev(:,1,:) = lw_k_dist%get_press_min() + epsilon(lw_k_dist%get_press_min())
else
  p_lev(:,nlay+1,:) = lw_k_dist%get_press_min() + epsilon(lw_k_dist%get_press_min())
endif
allocate(lw_flux_up(ncol, nlay+1, nexp))
allocate(lw_flux_dn(ncol, nlay+1, nexp))
call stop_on_err(source%alloc(ncol, nlay, lw_k_dist))
call stop_on_err(lw_optical_props%alloc_1scl(ncol, nlay, lw_k_dist))
nbnd = lw_k_dist%get_nband()
allocate(sfc_emis_spec(nbnd, ncol))
allocate(lw_heating_rate(ncol, nlay, nexp))


!Initialize shortwave.
sw_kdist_file = "coefficients_sw.nc"
if (nargs .ge. 3) call get_command_argument(3, sw_kdist_file)
call determine_gas_names(rfmip_file, sw_kdist_file, forcing_index, &
                         sw_kdist_gas_names, sw_rfmip_gas_games)
print *, "Sw calculation uses RFMIP gases: ", &
         (trim(sw_rfmip_gas_games(b))//" ", b = 1, size(sw_rfmip_gas_games))
call read_and_block_gases_ty(rfmip_file, ncol, sw_kdist_gas_names, &
                             sw_rfmip_gas_games, sw_gas_conc_array)
call load_and_init(sw_k_dist, trim(sw_kdist_file), sw_gas_conc_array(1))
if (.not. sw_k_dist%source_is_external()) then
  call stop_on_err("k-distribution file isn't SW.")
endif
nbnd = sw_k_dist%get_nband()
ngpt = sw_k_dist%get_ngpt()
allocate(toa_flux(ncol, ngpt))
allocate(def_tsi(ncol))
allocate(usecol(ncol, nexp))
if (top_at_1) then
  p_lev(:,1,:) = sw_k_dist%get_press_min() + epsilon(sw_k_dist%get_press_min())
else
  p_lev(:,nlay+1,:) = sw_k_dist%get_press_min() + epsilon(sw_k_dist%get_press_min())
endif
do b = 1, nexp
  usecol(1:ncol,b) = solar_zenith_angle(1:ncol,b) .lt. &
                           90._wp - 2._wp * spacing(90._wp)
enddo
allocate(sw_flux_up(ncol, nlay+1, nexp))
allocate(sw_flux_dn(ncol, nlay+1, nexp))
allocate(mu0(ncol))
allocate(sfc_alb_spec(nbnd, ncol))
call stop_on_err(sw_optical_props%alloc_2str(ncol, nlay, sw_k_dist))
allocate(sw_heating_rate(ncol, nlay, nexp))

!Create output files
do i = 1, nexp
  buffer = ""
  if (i-1 .lt. 10) then
    write(buffer, '(a,i1,a)') "fluxes_experiment_", i-1, ".nc"
  else
    write(buffer, '(a,i2,a)') "fluxes_experiment_", i-1, ".nc"
  endif
  call create_output_file(trim(buffer), ncol, nlay+1, output(i))
enddo

!Calculate present-day longwave heating rates.
call calculate_lw_fluxes(1, lw_fluxes, lw_flux_up, lw_flux_dn, lw_k_dist, sfc_emis_spec, &
                         sfc_emis, p_lay, p_lev, t_lay, sfc_t, lw_gas_conc_array, &
                         lw_optical_props, source, t_lev, top_at_1, ncol, n_quad_angles)
lw_heating_rate(:,:,:) = 0.0
call calculate_heating_rate(lw_flux_up(:,:,1), lw_flux_dn(:,:,1), p_lev(:,:,1), top_at_1, &
                            stratosphere_starting_index(:,1), stratosphere_ending_index(:,1), &
                            lw_heating_rate(:,:,1))

!Calculate present-day shortwave heating rates.
call calculate_sw_fluxes(1, sw_fluxes, sw_flux_up, sw_flux_dn, sw_k_dist, p_lay, p_lev, &
                         t_lay, sw_gas_conc_array, sw_optical_props, toa_flux, ncol, &
                         def_tsi, total_solar_irradiance, sfc_alb_spec, surface_albedo, &
                         mu0, solar_zenith_angle, usecol, top_at_1)
sw_heating_rate(:,:,:) = 0.0
call calculate_heating_rate(sw_flux_up(:,:,1), sw_flux_dn(:,:,1), p_lev(:,:,1), top_at_1, &
                            stratosphere_starting_index(:,1), stratosphere_ending_index(:,1), &
                            sw_heating_rate(:,:,1))
call write_output(output(1), lw_flux_dn(:,:,1), lw_flux_up(:,:,1), sw_flux_dn(:,:,1), &
                  sw_flux_up(:,:,1), lw_heating_rate(:,:,1), sw_heating_rate(:,:,1), &
                  t_lev(:,:,1), t_lay(:,:,1), 1)
call catch_netcdf_error(nf90_close(output(1)%ncid))

!Calculate the "dynamic" heating rate.
allocate(dynamic_heating_rate(ncol, nlay))
dynamic_heating_rate(:,:) = 0._wp
do i = 1, ncol
  do j = stratosphere_starting_index(i,1), stratosphere_ending_index(i,1)
    dynamic_heating_rate(i,j) = -1._wp*(lw_heating_rate(i,j,1) + sw_heating_rate(i,j,1))
  enddo
enddo

num_iterations = 250
do i = 1, num_iterations
  write(output_unit, "(a,i0.3,a,i0.3)") "Executing iteration ", i, "/", num_iterations
  do b = 2, nexp

    !Calculate longwave fluxes and heating rates.
    call calculate_lw_fluxes(b, lw_fluxes, lw_flux_up, lw_flux_dn, lw_k_dist, sfc_emis_spec, &
                             sfc_emis, p_lay, p_lev, t_lay, sfc_t, lw_gas_conc_array, &
                             lw_optical_props, source, t_lev, top_at_1, ncol, n_quad_angles)
    call calculate_heating_rate(lw_flux_up(:,:,b), lw_flux_dn(:,:,b), p_lev(:,:,b), top_at_1, &
                                stratosphere_starting_index(:,b), stratosphere_ending_index(:,b), &
                                lw_heating_rate(:,:,b))

    !Calculate shortwave fluxes and heating rates.
    call calculate_sw_fluxes(b, sw_fluxes, sw_flux_up, sw_flux_dn, sw_k_dist, p_lay, p_lev, &
                             t_lay, sw_gas_conc_array, sw_optical_props, toa_flux, ncol, &
                             def_tsi, total_solar_irradiance, sfc_alb_spec, surface_albedo, &
                             mu0, solar_zenith_angle, usecol, top_at_1)
    call calculate_heating_rate(sw_flux_up(:,:,b), sw_flux_dn(:,:,b), p_lev(:,:,b), top_at_1, &
                                stratosphere_starting_index(:,b), stratosphere_ending_index(:,b), &
                                sw_heating_rate(:,:,b))

    !Adjust temperatures.
    do j = 1, ncol
      do k = stratosphere_starting_index(j,b), stratosphere_ending_index(j,b)
        dq = dynamic_heating_rate(j,k) + lw_heating_rate(j,k,b) + sw_heating_rate(j,k,b)
        t_lay(j,k,b) = t_lay(j,k,b) + dq
      enddo

      !Interpolate heating rates to atmospheric levels.
      k = stratosphere_starting_index(j,b)
      t_lev(j,k,b) = t_lev(j,k,b) + dynamic_heating_rate(j,k) + lw_heating_rate(j,k,b) + &
                     sw_heating_rate(j,k,b)
      k = stratosphere_ending_index(j,b)
      t_lev(j,k+1,b) = t_lev(j,k+1,b) + dynamic_heating_rate(j,k) + lw_heating_rate(j,k,b) + &
                     sw_heating_rate(j,k,b)
      do k = stratosphere_starting_index(j,b), stratosphere_ending_index(j,b)-1
        call interpolate(p_lay(j,k,b), p_lay(j,k+1,b), dynamic_heating_rate(j,k), &
                         dynamic_heating_rate(j,k+1), p_lev(j,k+1,b), d_hr)
        call interpolate(p_lay(j,k,b), p_lay(j,k+1,b), lw_heating_rate(j,k,b), &
                         lw_heating_rate(j,k+1,b), p_lev(j,k+1,b), lw_hr)
        call interpolate(p_lay(j,k,b), p_lay(j,k+1,b), sw_heating_rate(j,k,b), &
                         sw_heating_rate(j,k+1,b), p_lev(j,k+1,b), sw_hr)
        t_lev(j,k+1,b) = t_lev(j,k+1,b) + lw_hr + sw_hr + d_hr
      enddo
    enddo

    !Write out data.
    call write_output(output(b), lw_flux_dn(:,:,b), lw_flux_up(:,:,b), sw_flux_dn(:,:,b), &
                      sw_flux_up(:,:,b), lw_heating_rate(:,:,b), sw_heating_rate(:,:,b), &
                      t_lev(:,:,b), t_lay(:,:,b), i)
  enddo
enddo

do i = 2, nexp
  call catch_netcdf_error(nf90_close(output(i)%ncid))
enddo

!Free memory.
deallocate(p_lay)
deallocate(p_lev)
deallocate(t_lay)
deallocate(t_lev)
deallocate(sfc_emis)
deallocate(sfc_t)
deallocate(surface_albedo)
deallocate(total_solar_irradiance)
deallocate(solar_zenith_angle)
deallocate(lw_kdist_gas_names)
deallocate(lw_rfmip_gas_games)
deallocate(lw_gas_conc_array)
deallocate(sw_kdist_gas_names)
deallocate(sw_rfmip_gas_games)
deallocate(sw_gas_conc_array)
deallocate(lw_flux_up)
deallocate(lw_flux_dn)
deallocate(sfc_emis_spec)
deallocate(toa_flux)
deallocate(def_tsi)
deallocate(usecol)
deallocate(sw_flux_up)
deallocate(sw_flux_dn)
deallocate(mu0)
deallocate(sfc_alb_spec)
deallocate(stratosphere_starting_index)
deallocate(stratosphere_ending_index)
deallocate(lw_heating_rate)
deallocate(sw_heating_rate)
deallocate(dynamic_heating_rate)


contains


subroutine stop_on_err(error_msg)
  character(len=*), intent(in) :: error_msg
  if (error_msg .ne. "") then
    write(error_unit, *) "Error: "//trim(error_msg)
    stop 1
  end if
end subroutine stop_on_err


subroutine interpolate(x1, x2, y1, y2, x, y)

  real(kind=wp), intent(in) :: x1
  real(kind=wp), intent(in) :: x2
  real(kind=wp), intent(in) :: y1
  real(kind=wp), intent(in) :: y2
  real(kind=wp), intent(in) :: x
  real(kind=wp), intent(out) :: y

  real(kind=wp) :: m
  real(kind=wp) :: b

  m = (y2 - y1)/(x2 - x1)
  b = y1 - m*x1
  y = m*x + b
end subroutine interpolate


subroutine calculate_heating_rate(flux_up, flux_down, p, top_at_1, layer_starting_index, &
                                  layer_ending_index, heating_rate)

  real(kind=wp), dimension(:,:), intent(in) :: flux_up !Upwelling flux [W/m^2]. (column, level).
  real(kind=wp), dimension(:,:), intent(in) :: flux_down !Downwelling flux [W/m^2]. (column, level).
  real(kind=wp), dimension(:,:), intent(in) :: p !Pressure at atmospheric levels [Pa]. (column, level).
  logical, intent(in) :: top_at_1 !Is toa the first level?.
  integer, dimension(:), intent(in) :: layer_starting_index !Starting layer that heating rates
                                                            !will be calculated for.
  integer, dimension(:), intent(in) :: layer_ending_index !Ending layer that heating rates
                                                          !will be calculated for.
  real(kind=wp), dimension(:,:), intent(inout) :: heating_rate ![K/day].

  integer :: ncol
  integer :: nlay
  integer :: i
  integer :: j
  real(kind=wp) :: dp ![Pa].
  real(kind=wp) :: rho ![kg/m^2].
  real(kind=wp) :: net_flux ![W/m^2].
  real(kind=wp), parameter :: cp = 1004.6_wp !Specific heat of air at constant pressure [J/(kg*K)].
  real(kind=wp), parameter :: g = 9.8_wp !Acceleration due to gravity [m/s].
  real(kind=wp), parameter :: spd = 86400._wp !Seconds per day [s/day].

  !Calculate the integrated number density across each layer.
  ncol = size(flux_up, 1)
  nlay = size(flux_up, 2) - 1
  do i = 1, ncol
    do j = layer_starting_index(i), layer_ending_index(i)
      dp = abs(p(i,j) - p(i,j+1))
      rho = dp/g
      net_flux = flux_down(i,j) - flux_up(i,j) - flux_down(i,j+1) + flux_up(i,j+1)
      if (.not. top_at_1) then
        net_flux = -1._wp*net_flux
      endif
      heating_rate(i,j) = (net_flux*spd)/(cp*rho)
    enddo
  enddo
end subroutine calculate_heating_rate


subroutine calculate_lw_fluxes(experiment, fluxes, flux_up, flux_dn, k_dist, sfc_emis_spec, &
                               sfc_emis, p_lay, p_lev, t_lay, sfc_t, gas_conc_array, &
                               optical_props, source, t_lev, top_at_1, ncol, n_quad_angles)

  integer, intent(in) :: experiment
  type(ty_fluxes_broadband), intent(inout) :: fluxes
  real(wp), dimension(:,:,:), intent(inout), target :: flux_up
  real(wp), dimension(:,:,:), intent(inout), target :: flux_dn
  type(ty_gas_optics_rrtmgp), intent(inout) :: k_dist
  real(wp), dimension(:,:), intent(inout) :: sfc_emis_spec
  real(wp), dimension(:,:), intent(in) :: sfc_emis
  real(wp), dimension(:,:,:), intent(in) :: p_lay
  real(wp), dimension(:,:,:), intent(in) :: p_lev
  real(wp), dimension(:,:,:), intent(in) :: t_lay
  real(wp), dimension(:,:), intent(in) :: sfc_t
  type(ty_gas_concs), dimension(:), intent(in) :: gas_conc_array
  type(ty_optical_props_1scl), intent(inout) :: optical_props
  type(ty_source_func_lw), intent(inout) :: source
  real(wp), dimension(:,:,:), intent(in) :: t_lev
  logical, intent(in) :: top_at_1
  integer, intent(in) :: ncol
  integer, intent(in) :: n_quad_angles

  integer :: b
  integer :: nbnd
  integer :: icol
  integer :: ibnd

  b = experiment
  nbnd = k_dist%get_nband()
  fluxes%flux_up => flux_up(:,:,b)
  fluxes%flux_dn => flux_dn(:,:,b)
  do icol = 1, ncol
    do ibnd = 1, nbnd
      sfc_emis_spec(ibnd,icol) = sfc_emis(icol,b)
    enddo
  enddo
  call stop_on_err(k_dist%gas_optics(p_lay(:,:,b), p_lev(:,:,b), t_lay(:,:,b), sfc_t(:,b), &
                                     gas_conc_array(b), optical_props, source, &
                                     tlev=t_lev(:,:,b)))
  call stop_on_err(rte_lw(optical_props, top_at_1, source, sfc_emis_spec, fluxes, &
                          n_gauss_angles=n_quad_angles))
end subroutine calculate_lw_fluxes


subroutine calculate_sw_fluxes(experiment, fluxes, flux_up, flux_dn, k_dist, p_lay, p_lev, &
                               t_lay, gas_conc_array, optical_props, toa_flux, ncol, &
                               def_tsi, total_solar_irradiance, sfc_alb_spec, surface_albedo, &
                               mu0, solar_zenith_angle, usecol, top_at_1)

  integer, intent(in) :: experiment
  type(ty_fluxes_broadband), intent(inout) :: fluxes
  real(wp), dimension(:,:,:), intent(inout), target :: flux_up
  real(wp), dimension(:,:,:), intent(inout), target :: flux_dn
  type(ty_gas_optics_rrtmgp), intent(inout) :: k_dist
  real(wp), dimension(:,:,:), intent(in) :: p_lay
  real(wp), dimension(:,:,:), intent(in) :: p_lev
  real(wp), dimension(:,:,:), intent(in) :: t_lay
  type(ty_gas_concs), dimension(:), intent(in) :: gas_conc_array
  type(ty_optical_props_2str), intent(inout) :: optical_props
  real(wp), dimension(:,:), intent(inout) :: toa_flux
  integer, intent(in) :: ncol
  real(wp), dimension(:), intent(inout) :: def_tsi
  real(wp), dimension(:,:), intent(in) :: total_solar_irradiance
  real(wp), dimension(:,:), intent(inout) :: sfc_alb_spec
  real(wp), dimension(:,:), intent(in) :: surface_albedo
  real(wp), dimension(:), intent(inout) :: mu0
  real(wp), dimension(:,:), intent(in) :: solar_zenith_angle
  logical, dimension(:,:), intent(in) :: usecol
  logical, intent(in) :: top_at_1

  integer :: b
  integer :: nbnd
  integer :: icol
  integer :: igpt
  integer :: ibnd
  integer :: ngpt
  real(wp), parameter :: deg_to_rad = acos(-1._wp)/180._wp

  b = experiment
  nbnd = k_dist%get_nband()
  ngpt = k_dist%get_ngpt()
  fluxes%flux_up => flux_up(:,:,b)
  fluxes%flux_dn => flux_dn(:,:,b)
  call stop_on_err(k_dist%gas_optics(p_lay(:,:,b), p_lev(:,:,b), t_lay(:,:,b), &
                                     gas_conc_array(b), optical_props, toa_flux))
  do icol = 1, ncol
    def_tsi(icol) = toa_flux(icol, 1)
  enddo
  do igpt = 1, ngpt
    do icol = 1, ncol
      def_tsi(icol) = def_tsi(icol) + toa_flux(icol, igpt)
    enddo
  enddo
  do igpt = 1, ngpt
    do icol = 1, ncol
      toa_flux(icol,igpt) = toa_flux(icol,igpt) * total_solar_irradiance(icol,b)/def_tsi(icol)
    enddo
  enddo
  do icol = 1, ncol
    do ibnd = 1, nbnd
      sfc_alb_spec(ibnd,icol) = surface_albedo(icol,b)
    enddo
  enddo
  do icol = 1, ncol
    mu0(icol) = merge(cos(solar_zenith_angle(icol,b)*deg_to_rad), 1._wp, usecol(icol,b))
  end do
  call stop_on_err(rte_sw(optical_props, top_at_1, mu0, toa_flux, sfc_alb_spec, &
                          sfc_alb_spec, fluxes))
  do icol = 1, ncol
    if (.not. usecol(icol,b)) then
      flux_up(icol,:,b) = 0._wp
      flux_dn(icol,:,b) = 0._wp
    endif
  enddo
end subroutine calculate_sw_fluxes


subroutine catch_netcdf_error(code)

  integer, intent(in) :: code

  character(len=80) :: buffer

  if (code .ne. nf90_noerr) then
    buffer = nf90_strerror(code)
    call stop_on_err(trim(buffer))
  endif
end subroutine catch_netcdf_error


subroutine switch_col_and_z_2d(in_var, out_var)

  real(kind=wp), dimension(:,:), intent(in) :: in_var
  real(kind=wp), dimension(:,:), allocatable, intent(inout) :: out_var

  integer :: i
  integer :: j

  if (allocated(out_var)) deallocate(out_var)
  allocate(out_var(size(in_var, 2), size(in_var, 1)))
  do j = 1, size(in_var, 1)
    do i = 1, size(in_var, 2)
      out_var(i,j) = in_var(j,i)
    enddo
  enddo
end subroutine switch_col_and_z_2d


function add_variable(ncid, name, standard_name, units, dimid) result(varid)

  integer, intent(in) :: ncid
  character(len=*), intent(in) :: name
  character(len=*), intent(in) :: standard_name
  character(len=*), intent(in) :: units
  integer, dimension(:), intent(in) :: dimid
  integer :: varid

  call catch_netcdf_error(nf90_def_var(ncid, trim(name), nf90_double, dimid, varid))
  call catch_netcdf_error(nf90_put_att(ncid, varid, "units", trim(units)))
  call catch_netcdf_error(nf90_put_att(ncid, varid, "standard_name", trim(standard_name)))
end function add_variable


subroutine create_output_file(path, ncol, nlev, output)

  character(len=*), intent(in) :: path
  integer, intent(in) :: ncol
  integer, intent(in) :: nlev
  type(OutputFile), intent(inout) :: output

  integer, dimension(4) :: dimid

  call catch_netcdf_error(nf90_create(path, ior(nf90_clobber, nf90_netcdf4), output%ncid))
  call catch_netcdf_error(nf90_def_dim(output%ncid, "expt", 1, output%dimid(expt)))
  call catch_netcdf_error(nf90_def_dim(output%ncid, "site", ncol, output%dimid(site)))
  call catch_netcdf_error(nf90_def_dim(output%ncid, "layer", nlev-1, output%dimid(layer)))
  call catch_netcdf_error(nf90_def_dim(output%ncid, "level", nlev, output%dimid(level)))
  call catch_netcdf_error(nf90_def_dim(output%ncid, "time", nf90_unlimited, output%dimid(time)))
  dimid(1) = output%dimid(level)
  dimid(2) = output%dimid(site)
  dimid(3) = output%dimid(expt)
  dimid(4) = output%dimid(time)
  output%varid(rlu) = add_variable(output%ncid, "rlu", "upwelling_longwave_flux_in_air", &
                                   "W m-2", dimid)
  output%varid(rld) = add_variable(output%ncid, "rld", "downwelling_longwave_flux_in_air", &
                                   "W m-2", dimid)
  output%varid(rsu) = add_variable(output%ncid, "rsu", "upwelling_shortwave_flux_in_air", &
                                   "W m-2", dimid)
  output%varid(rsd) = add_variable(output%ncid, "rsd", "downwelling_shortwave_flux_in_air", &
                                   "W m-2", dimid)
  output%varid(tlevid) = add_variable(output%ncid, "tlev", "level_temperature", "K", dimid)
  dimid(1) = output%dimid(layer)
  output%varid(rlhr) = add_variable(output%ncid, "rlhr", "longwave_heating_rate_in_air", &
                                    "K day-1", dimid)
  output%varid(rshr) = add_variable(output%ncid, "rshr", "shortwave_heating_rate_in_air", &
                                    "K day-1", dimid)
  output%varid(tlayid) = add_variable(output%ncid, "tlay", "layer_temperature", "K", dimid)
end subroutine create_output_file


subroutine write_output(output, lw_flux_dn, lw_flux_up, sw_flux_dn, sw_flux_up, &
                        lw_heating_rate, sw_heating_rate, t_lev, t_lay, iter)

  type(OutputFile), intent(in) :: output
  real(kind=wp), dimension(:,:), intent(in) :: lw_flux_dn
  real(kind=wp), dimension(:,:), intent(in) :: lw_flux_up
  real(kind=wp), dimension(:,:), intent(in) :: sw_flux_dn
  real(kind=wp), dimension(:,:), intent(in) :: sw_flux_up
  real(kind=wp), dimension(:,:), intent(in) :: lw_heating_rate
  real(kind=wp), dimension(:,:), intent(in) :: sw_heating_rate
  real(kind=wp), dimension(:,:), intent(in) :: t_lev
  real(kind=wp), dimension(:,:), intent(in) :: t_lay
  integer, intent(in) :: iter

  real(kind=wp), dimension(:,:), allocatable :: buffer
  integer, dimension(4) :: start
  integer, dimension(4) :: counts
  integer :: ncol
  integer :: nlev

  ncol = size(lw_flux_dn, 1)
  nlev = size(lw_flux_dn, 2)
  start(1:3) = 1
  start(4) = iter
  counts(1) = nlev
  counts(2) = ncol
  counts(3) = 1
  counts(4) = 1
  call switch_col_and_z_2d(lw_flux_dn, buffer)
  call catch_netcdf_error(nf90_put_var(output%ncid, output%varid(rld), buffer, &
                                       start=start, count=counts))
  call switch_col_and_z_2d(lw_flux_up, buffer)
  call catch_netcdf_error(nf90_put_var(output%ncid, output%varid(rlu), buffer, &
                                       start=start, count=counts))
  call switch_col_and_z_2d(sw_flux_dn, buffer)
  call catch_netcdf_error(nf90_put_var(output%ncid, output%varid(rsd), buffer, &
                                       start=start, count=counts))
  call switch_col_and_z_2d(sw_flux_up, buffer)
  call catch_netcdf_error(nf90_put_var(output%ncid, output%varid(rsu), buffer, &
                                       start=start, count=counts))
  call switch_col_and_z_2d(t_lev, buffer)
  call catch_netcdf_error(nf90_put_var(output%ncid, output%varid(tlevid), buffer, &
                                       start=start, count=counts))
  counts(1) = nlev - 1
  call switch_col_and_z_2d(lw_heating_rate, buffer)
  call catch_netcdf_error(nf90_put_var(output%ncid, output%varid(rlhr), buffer, &
                                       start=start, count=counts))
  call switch_col_and_z_2d(sw_heating_rate, buffer)
  call catch_netcdf_error(nf90_put_var(output%ncid, output%varid(rshr), buffer, &
                                       start=start, count=counts))
  call switch_col_and_z_2d(t_lay, buffer)
  call catch_netcdf_error(nf90_put_var(output%ncid, output%varid(tlayid), buffer, &
                                       start=start, count=counts))
  deallocate(buffer)
end subroutine write_output


end program stratospheric_adjustments
