program clearsky_rrtmgp

use, intrinsic :: iso_fortran_env, only: error_unit, output_unit, real32
use netcdf

use mo_fluxes, only: ty_fluxes_broadband
use mo_gas_concentrations, only: ty_gas_concs
use mo_gas_optics_rrtmgp, only: ty_gas_optics_rrtmgp
use mo_optical_props, only: ty_optical_props_1scl, ty_optical_props_2str
use mo_rte_kind, only: wp
use mo_rte_lw, only: rte_lw
use mo_source_functions, only: ty_source_func_lw
use mo_load_coefficients, only: load_and_init
use mo_clr_io, only: read_atmos, write_lw_fluxes
use mo_rfmip_io, only: determine_gas_names, read_and_block_gases_ty, read_and_block_lw_bc, &
                       read_and_block_pt, read_and_block_sw_bc, read_size
implicit none


character(len=132) :: input_file
character(len=64) :: name_of_program
character(len=64) :: buffer
integer :: nargs
integer :: ncol
integer :: nlay, nlev
integer :: n_quad_angles
logical :: top_at_1
integer :: nbnd
integer :: ngpt
real(kind=wp), dimension(:,:), allocatable :: p_lay !Atmospheric layer pressure [Pa] (ncol, nlay).
real(kind=wp), dimension(:,:), allocatable :: p_lev !Atmospheric level pressure [Pa] (ncol, nlev).
real(kind=wp), dimension(:,:), allocatable :: t_lay !Atmospheric layer temperature [K] (ncol, nlay).
real(kind=wp), dimension(:,:), allocatable :: t_lev !Atmoshperic level temperature [K] (ncol, nlev).
real(kind=wp), dimension(:,:), allocatable :: col_dry
real(kind=wp), dimension(:), allocatable :: sfc_t !Surface temperature [K] (ncol).
real(kind=wp), dimension(:,:), allocatable :: sfc_emis_spec !Surface emissivity (nbnd, ncol).
type(ty_source_func_lw) :: source,  source_trop_off, source_strat_off, source_air_off !Source function object.
real(kind=wp), dimension(:,:), allocatable :: toa_flux !Shortwave flux at top-of-atmosphere [W/m^2] (ncol, ngpt).
integer, parameter :: ngas = 8
real(kind=wp), parameter :: sfc_emis = 1 !Surface emissivity
real(kind=wp), parameter :: stratosphere_max_pressure = 20000._wp !Bottom of the stratosphere [Pa].
integer, dimension(:), allocatable :: stratosphere_starting_index !Index corresponding to the top, bottom of
                                                                    !the stratosphere if the column is oriented
                                                                    !toward, away from the surface (ncol, nexp).
integer, dimension(:), allocatable :: stratosphere_ending_index !Index corresponding to the bottom, top of
                                                                  !the stratosphere if the column is oriented
                                                                  !toward, away from the surface (ncol, nexp).
character(len=3), dimension(ngas) &
                     :: gas_names = ['h2o', 'co2', 'o3 ', 'n2o', 'co ', 'ch4', 'o2 ','n2 ']

character(len=132) :: lw_kdist_file !K-distribution configuration for the longwave.
real(kind=wp), dimension(:,:), allocatable, target :: lw_flux_up, lw_flux_up_trop_off, lw_flux_up_strat_off, lw_flux_up_air_off !Upwelling longwave flux [W/m^2] (ncol, nlev).
real(kind=wp), dimension(:,:), allocatable, target :: lw_flux_dn !Downwelling longwave flux [W/m^2] (ncol, nlev).
type(ty_gas_optics_rrtmgp) :: lw_k_dist !Longwave k-distribution object.
type(ty_optical_props_1scl) :: lw_optical_props !Longwave optics object.
type(ty_fluxes_broadband) :: lw_fluxes, lw_fluxes_trop_off, lw_fluxes_strat_off, lw_fluxes_air_off !Longwave fluxes object.
type(ty_gas_concs) :: gas_conc !Longwave gas concentrations object.
real(kind=wp), dimension(:,:), allocatable :: lw_heating_rate !Longwave heating rate [K/day] (ncol, nlay).
real(kind=wp) :: dq
real(kind=wp) :: d_hr
real(kind=wp) :: lw_hr
integer :: num_iterations = 5
integer :: b
integer :: i
integer :: j
integer :: k
integer :: igas

!Parse command line arguments.
call get_command_argument(0, name_of_program)
nargs = command_argument_count()
do i = 1, nargs
  call get_command_argument(i, buffer)
  if (trim(buffer) .eq. "-h" .or. trim(buffer) .eq. "--help") then
    write(output_unit, "(a)") "Usage: "//trim(name_of_program)//" [input_file]" &
         //" [lw k-distribution_file]"
    stop
  endif
enddo

input_file = "input.nc"
if (nargs .ge. 1) then
  call get_command_argument(1, input_file)
endif
call read_atmos(input_file, ncol, nlay, nlev, &
                p_lay, t_lay, p_lev, t_lev, sfc_t, &
                gas_conc, col_dry)
n_quad_angles = 3
top_at_1 = p_lay(1,1) .lt. p_lay(1,nlay)

!Find pressure layer indices for the stratosphere.
allocate(stratosphere_starting_index(ncol))
allocate(stratosphere_ending_index(ncol))
if (top_at_1) then
  stratosphere_starting_index(:) = 1
else
  stratosphere_starting_index(:) = nlay
endif

do i = 1, ncol
  if (top_at_1) then
    !Leaving the stratosphere during downward sweep through the layers.
    stratosphere_ending_index(i) = minloc(t_lev(i,5:), dim=1)
  else
    !Entering the stratosphere during upward sweep through the layers.
    stratosphere_ending_index(i) = minloc(t_lev(i,:nlev-4), dim=1)
  endif
enddo

!Initialize longwave.
lw_kdist_file = "coefficients_lw.nc"
if (nargs .ge. 2) then
  call get_command_argument(2, lw_kdist_file)
endif
call load_and_init(lw_k_dist, trim(lw_kdist_file), gas_conc)
if (.not. lw_k_dist%source_is_internal()) then
  call stop_on_err("k-distribution file isn't LW.")
endif
if (top_at_1) then
  p_lev(:,1) = lw_k_dist%get_press_min() + epsilon(lw_k_dist%get_press_min())
else
  p_lev(:,nlay+1) = lw_k_dist%get_press_min() + epsilon(lw_k_dist%get_press_min())
endif
allocate(lw_flux_up(ncol, nlay+1))
allocate(lw_flux_up_trop_off(ncol, nlay+1))
allocate(lw_flux_up_strat_off(ncol, nlay+1))
allocate(lw_flux_up_air_off(ncol, nlay+1))
allocate(lw_flux_dn(ncol, nlay+1))
call stop_on_err(source%alloc(ncol, nlay, lw_k_dist))
call stop_on_err(source_trop_off%alloc(ncol, nlay, lw_k_dist))
call stop_on_err(source_strat_off%alloc(ncol, nlay, lw_k_dist))
call stop_on_err(source_air_off%alloc(ncol, nlay, lw_k_dist))
call stop_on_err(lw_optical_props%alloc_1scl(ncol, nlay, lw_k_dist))
nbnd = lw_k_dist%get_nband()
allocate(sfc_emis_spec(nbnd, ncol))
allocate(lw_heating_rate(ncol, nlay))

!Calculate present-day longwave heating rates.
call calculate_lw_fluxes(lw_fluxes, lw_fluxes_trop_off, lw_fluxes_strat_off, lw_fluxes_air_off, &
                         lw_flux_up, lw_flux_up_trop_off, lw_flux_up_strat_off, lw_flux_up_air_off, &
                         lw_flux_dn, lw_k_dist, sfc_emis_spec, &
                         sfc_emis, p_lay, p_lev, t_lay, sfc_t, gas_conc, &
                         lw_optical_props, source, source_trop_off, source_strat_off, source_air_off, &
                         t_lev, top_at_1, &
                         stratosphere_starting_index(:), stratosphere_ending_index(:), &
                         ncol, n_quad_angles)
lw_heating_rate(:,:) = 0.0
call calculate_heating_rate(lw_flux_up(:,:), lw_flux_dn(:,:), p_lev(:,:), &
                            lw_heating_rate(:,:))
call write_lw_fluxes(input_file, lw_flux_up, lw_flux_dn, lw_heating_rate, &
                     lw_flux_up_trop_off, &
                     lw_flux_up_strat_off, &
                     lw_flux_up_air_off)
!Free memory.
deallocate(p_lay)
deallocate(p_lev)
deallocate(t_lay)
deallocate(t_lev)
deallocate(sfc_t)
deallocate(lw_flux_up)
deallocate(lw_flux_up_trop_off)
deallocate(lw_flux_up_strat_off)
deallocate(lw_flux_up_air_off)
deallocate(lw_flux_dn)
deallocate(sfc_emis_spec)
deallocate(stratosphere_starting_index)
deallocate(stratosphere_ending_index)
deallocate(lw_heating_rate)


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


subroutine calculate_heating_rate(flux_up, flux_down, p, heating_rate)

  real(kind=wp), dimension(:,:), intent(in) :: flux_up !Upwelling flux [W/m^2]. (column, level).
  real(kind=wp), dimension(:,:), intent(in) :: flux_down !Downwelling flux [W/m^2]. (column, level).
  real(kind=wp), dimension(:,:), intent(in) :: p !Pressure at atmospheric levels [Pa]. (column, level).
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
    do j = 1, nlay
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


subroutine calculate_lw_fluxes(fluxes, fluxes_trop_off, fluxes_strat_off, fluxes_air_off, &
                               flux_up, flux_up_trop_off, flux_up_strat_off, flux_up_air_off, &
                               flux_dn, k_dist, sfc_emis_spec, &
                               sfc_emis, p_lay, p_lev, t_lay, sfc_t, gas_conc, &
                               optical_props, source, source_trop_off, source_strat_off, source_air_off, &
                               t_lev, top_at_1, layer_starting_index, layer_ending_index, ncol, n_quad_angles)

  type(ty_fluxes_broadband), intent(inout) :: fluxes, fluxes_trop_off, fluxes_strat_off, fluxes_air_off
  real(wp), dimension(:,:), intent(inout), target :: flux_up, flux_up_trop_off, flux_up_strat_off, flux_up_air_off
  real(wp), dimension(:,:), intent(inout), target :: flux_dn
  type(ty_gas_optics_rrtmgp), intent(inout) :: k_dist
  real(wp), dimension(:,:), intent(inout) :: sfc_emis_spec
  real(wp), intent(in) :: sfc_emis
  real(wp), dimension(:,:), intent(in) :: p_lay
  real(wp), dimension(:,:), intent(in) :: p_lev
  real(wp), dimension(:,:), intent(in) :: t_lay
  real(wp), dimension(:), intent(in) :: sfc_t
  type(ty_gas_concs), intent(in) :: gas_conc
  type(ty_optical_props_1scl), intent(inout) :: optical_props
  type(ty_source_func_lw), intent(inout) :: source, source_trop_off, source_strat_off, source_air_off
  real(wp), dimension(:,:), intent(in) :: t_lev
  logical, intent(in) :: top_at_1
  integer, intent(in) :: ncol
  integer, intent(in) :: n_quad_angles
  integer, dimension(:), intent(in) :: layer_starting_index !Starting layer that heating rates
                                                            !will be calculated for.
  integer, dimension(:), intent(in) :: layer_ending_index !Ending layer that heating rates
                                                          !will be calculated for.
  integer :: nbnd
  integer :: icol
  integer :: ibnd

  nbnd = k_dist%get_nband()
  fluxes%flux_up => flux_up
  fluxes_trop_off%flux_up => flux_up_trop_off
  fluxes_strat_off%flux_up => flux_up_strat_off
  fluxes_air_off%flux_up => flux_up_air_off
  fluxes%flux_dn => flux_dn
  do icol = 1, ncol
    do ibnd = 1, nbnd
      sfc_emis_spec(ibnd,icol) = sfc_emis
    enddo
  enddo
  call stop_on_err(k_dist%gas_optics(p_lay, p_lev, t_lay, sfc_t, &
                                     gas_conc, optical_props, source, &
                                     tlev=t_lev))
  source_trop_off%sfc_source = source%sfc_source
  source_trop_off%lev_source_inc = source%lev_source_inc
  source_trop_off%lev_source_dec = source%lev_source_dec
  source_trop_off%lay_source = source%lay_source

  source_strat_off%sfc_source = source%sfc_source
  source_strat_off%lev_source_inc = source%lev_source_inc
  source_strat_off%lev_source_dec = source%lev_source_dec
  source_strat_off%lay_source = source%lay_source

  source_air_off%sfc_source = source%sfc_source
  source_air_off%lev_source_inc = 0.
  source_air_off%lev_source_dec = 0.
  source_air_off%lay_source = 0.

  do icol = 1, ncol
    source_strat_off%lev_source_inc(icol, layer_starting_index(icol):layer_ending_index(icol)+1,:) = 0.
    source_strat_off%lev_source_dec(icol, layer_starting_index(icol):layer_ending_index(icol)+1,:) = 0.
    source_strat_off%lay_source(icol, layer_starting_index(icol):layer_ending_index(icol),:) = 0.                 
    if (top_at_1) then
      source_trop_off%lev_source_inc(icol, layer_ending_index(icol)+1:size(t_lev, 2),:) = 0.
      source_trop_off%lev_source_dec(icol, layer_ending_index(icol)+1:size(t_lev, 2),:) = 0.
      source_trop_off%lay_source(icol, layer_ending_index(icol)+1:size(t_lay, 2),:) = 0.
    else
      source_trop_off%lev_source_inc(icol, 1:layer_starting_index(icol),:) = 0.
      source_trop_off%lev_source_dec(icol, 1:layer_starting_index(icol),:) = 0.
      source_trop_off%lay_source(icol, 1:layer_starting_index(icol)-1,:) = 0.
    endif
  enddo

  call stop_on_err(rte_lw(optical_props, top_at_1, source_trop_off, sfc_emis_spec, fluxes_trop_off, &
                          n_gauss_angles=n_quad_angles))
  call stop_on_err(rte_lw(optical_props, top_at_1, source_strat_off, sfc_emis_spec, fluxes_strat_off, &
                          n_gauss_angles=n_quad_angles))
  call stop_on_err(rte_lw(optical_props, top_at_1, source_air_off, sfc_emis_spec, fluxes_air_off, &
                          n_gauss_angles=n_quad_angles))
  call stop_on_err(rte_lw(optical_props, top_at_1, source, sfc_emis_spec, fluxes, &
                          n_gauss_angles=n_quad_angles))
end subroutine calculate_lw_fluxes

subroutine catch_netcdf_error(code)

  integer, intent(in) :: code

  character(len=80) :: buffer

  if (code .ne. nf90_noerr) then
    buffer = nf90_strerror(code)
    call stop_on_err(trim(buffer))
  endif
end subroutine catch_netcdf_error


end program clearsky_rrtmgp
