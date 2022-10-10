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
use mo_clr_io_kernel, only: read_atmos, write_lw_fluxes, write_kernels
use mo_rfmip_io, only: determine_gas_names, read_and_block_gases_ty, read_and_block_lw_bc, &
                       read_and_block_pt, read_and_block_sw_bc, read_size
implicit none


character(len=132) :: input_file, input_file_pert_lay, input_file_pert_layp1
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
real(kind=wp), dimension(:,:), allocatable :: t_lay, t_lay_pert !Atmospheric layer temperature [K] (ncol, nlay).
real(kind=wp), dimension(:,:), allocatable :: t_lev, t_lev_pert_lay, t_lev_pert_layp1 !Atmoshperic level temperature [K] (ncol, nlev).
real(kind=wp), dimension(:,:), allocatable :: col_dry
real(kind=wp), dimension(:), allocatable :: sfc_t, sfc_t_pert !Surface temperature [K] (ncol).
real(kind=wp), dimension(:,:), allocatable :: sfc_emis_spec !Surface emissivity (nbnd, ncol).
type(ty_source_func_lw) :: source !Source function object.
real(kind=wp), dimension(:,:), allocatable :: toa_flux !Shortwave flux at top-of-atmosphere [W/m^2] (ncol, ngpt).
integer, parameter :: ngas = 8
real(kind=wp), parameter :: sfc_emis = 1 !Surface emissivity
character(len=3), dimension(ngas) &
                     :: gas_names = ['h2o', 'co2', 'o3 ', 'n2o', 'co ', 'ch4', 'o2 ','n2 ']

character(len=132) :: lw_kdist_file !K-distribution configuration for the longwave.
real(kind=wp), dimension(:,:), allocatable, target :: lw_flux_up
real(kind=wp), dimension(:,:), allocatable :: lw_flux_up_surf_pert !Upwelling longwave flux [W/m^2] (ncol, nlev).
real(kind=wp), dimension(:,:), allocatable, target :: lw_flux_dn 
real(kind=wp), dimension(:,:), allocatable :: lw_flux_dn_surf_pert !Downwelling longwave flux [W/m^2] (ncol, nlev).
real(kind=wp), dimension(:,:,:), allocatable :: lw_flux_up_gas_pert, lw_flux_up_T_pert !Upwelling longwave flux [W/m^2] (ncol, nlev, nlay).
real(kind=wp), dimension(:,:,:), allocatable :: lw_flux_dn_gas_pert, lw_flux_dn_T_pert !Downwelling longwave flux [W/m^2] (ncol, nley, nlay).
real(kind=wp), dimension(:,:), allocatable :: kernel_srf_up
real(kind=wp), dimension(:,:,:), allocatable :: kernel_q_up, kernel_t_up !kernel (ncol, nlev, nlay).
real(kind=wp), dimension(:,:,:), allocatable :: kernel_q_dn, kernel_t_dn !
type(ty_gas_optics_rrtmgp) :: lw_k_dist !Longwave k-distribution object.
type(ty_optical_props_1scl) :: lw_optical_props !Longwave optics object.
type(ty_fluxes_broadband) :: lw_fluxes !Longwave fluxes object.
type(ty_gas_concs) :: gas_conc, gas_conc_pert, gas_conc_1d !Longwave gas concentrations object.
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
    write(output_unit, "(a)") "Usage: "//trim(name_of_program)//" [input_file]"//" [input_file_pert]" &
         //" [lw k-distribution_file]"
    stop
  endif
enddo

input_file = "input.nc"
if (nargs .ge. 1) call get_command_argument(1, input_file)
if (nargs .ge. 2) call get_command_argument(2, input_file_pert_lay)
if (nargs .ge. 3) call get_command_argument(3, input_file_pert_layp1)
  write(output_unit, "(a)") "Reading: "//trim(input_file)
  call read_atmos(input_file, ncol, nlay, nlev, &
                  p_lay, t_lay, p_lev, t_lev, sfc_t, &
                  gas_conc, gas_conc_1d, col_dry)
  write(output_unit, "(a)") "Reading: "//trim(input_file_pert_lay)
  call read_atmos(input_file_pert_lay, ncol, nlay, nlev, &
                  p_lay, t_lay_pert, p_lev, t_lev_pert_lay, sfc_t_pert, &
                  gas_conc_pert, gas_conc_1d, col_dry)
  write(output_unit, "(a)") "Reading: "//trim(input_file_pert_layp1)
  call read_atmos(input_file_pert_layp1, ncol, nlay, nlev, &
                  p_lay, t_lay_pert, p_lev, t_lev_pert_layp1, sfc_t_pert, &
                  gas_conc_pert, gas_conc_1d, col_dry)
n_quad_angles = 3
top_at_1 = p_lay(1,1) .lt. p_lay(1,nlay)

!Initialize longwave.
lw_kdist_file = "coefficients_lw.nc"
if (nargs .ge. 4) call get_command_argument(4, lw_kdist_file)
call load_and_init(lw_k_dist, trim(lw_kdist_file), gas_conc_1d)
if (.not. lw_k_dist%source_is_internal()) then
  call stop_on_err("k-distribution file isn't LW.")
endif
if (top_at_1) then
  p_lev(:,1) = lw_k_dist%get_press_min() + epsilon(lw_k_dist%get_press_min())
else
  p_lev(:,nlay+1) = lw_k_dist%get_press_min() + epsilon(lw_k_dist%get_press_min())
endif
allocate(lw_flux_up(ncol, nlay+1))
allocate(lw_flux_dn(ncol, nlay+1))
allocate(lw_flux_up_surf_pert(ncol, nlay+1))
allocate(lw_flux_dn_surf_pert(ncol, nlay+1))
allocate(lw_flux_dn_T_pert(ncol, nlay+1, nlay))
allocate(lw_flux_up_T_pert(ncol, nlay+1, nlay))
allocate(lw_flux_dn_gas_pert(ncol, nlay+1, nlay))
allocate(lw_flux_up_gas_pert(ncol, nlay+1, nlay))
allocate(kernel_srf_up(ncol, nlay+1))
allocate(kernel_t_up(ncol, nlay+1, nlay))
allocate(kernel_t_dn(ncol, nlay+1, nlay))
allocate(kernel_q_up(ncol, nlay+1, nlay))
allocate(kernel_q_dn(ncol, nlay+1, nlay))

call stop_on_err(source%alloc(ncol, nlay, lw_k_dist))
call stop_on_err(lw_optical_props%alloc_1scl(ncol, nlay, lw_k_dist))
nbnd = lw_k_dist%get_nband()
allocate(sfc_emis_spec(nbnd, ncol))
allocate(lw_heating_rate(ncol, nlay))

!Calculate present-day longwave heating rates.
call calculate_lw_fluxes(lw_fluxes, &
                         lw_flux_up, lw_flux_up_surf_pert, lw_flux_up_gas_pert, lw_flux_up_T_pert, &
                         lw_flux_dn, lw_flux_dn_surf_pert, lw_flux_dn_gas_pert, lw_flux_dn_T_pert, &
                         lw_k_dist, sfc_emis_spec, &
                         sfc_emis, p_lay, p_lev, t_lay, t_lay_pert, sfc_t,sfc_t_pert, gas_conc, gas_conc_pert, &
                         lw_optical_props, source, &
                         t_lev, t_lev_pert_lay, t_lev_pert_layp1, top_at_1, &
                         ncol, n_quad_angles)                                                                                               
lw_heating_rate(:,:) = 0.0

    kernel_srf_up = lw_flux_up_surf_pert - lw_flux_up
  do i = 1, nlay
    kernel_t_up(:, :, i)   = lw_flux_up_T_pert(:, :, i)   - lw_flux_up(:,:)
    kernel_t_dn(:, :, i)   = lw_flux_dn_T_pert(:, :, i)   - lw_flux_dn(:,:)
    kernel_q_up(:, :, i)   = lw_flux_up_gas_pert(:, :, i) - lw_flux_up(:,:)
    kernel_q_dn(:, :, i)   = lw_flux_dn_gas_pert(:, :, i) - lw_flux_dn(:,:)
enddo

call calculate_heating_rate(lw_flux_up(:,:), lw_flux_dn(:,:), p_lev(:,:), &
                            lw_heating_rate(:,:))
call write_lw_fluxes(input_file, lw_flux_up, lw_flux_dn)
call write_kernels(input_file_pert_lay, kernel_srf_up, kernel_t_up, kernel_q_up, kernel_t_dn, kernel_q_dn)
!Free memory.
deallocate(p_lay)
deallocate(p_lev)
deallocate(t_lay, t_lay_pert)
deallocate(t_lev, t_lev_pert_lay, t_lev_pert_layp1)
deallocate(sfc_t, sfc_t_pert)
deallocate(lw_flux_up)
deallocate(lw_flux_dn)
deallocate(sfc_emis_spec)
deallocate(lw_heating_rate)
deallocate(kernel_t_dn)
deallocate(kernel_q_dn)
deallocate(kernel_t_up)
deallocate(kernel_q_up)
deallocate(kernel_srf_up)
deallocate(lw_flux_up_surf_pert)
deallocate(lw_flux_dn_surf_pert)
deallocate(lw_flux_dn_T_pert)
deallocate(lw_flux_up_T_pert)
deallocate(lw_flux_dn_gas_pert)
deallocate(lw_flux_up_gas_pert)

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


subroutine calculate_heating_rate(flux_up, flux_down, p, &
                                  heating_rate)

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


subroutine calculate_lw_fluxes(fluxes, &
                               flux_up, flux_up_surf_pert, flux_up_gas_pert, flux_up_T_pert, &
                               flux_dn, flux_dn_surf_pert, flux_dn_gas_pert, flux_dn_T_pert, &
                               k_dist, sfc_emis_spec, &
                               sfc_emis, p_lay, p_lev, t_lay, t_lay_pert, sfc_t, sfc_t_pert, gas_conc, gas_conc_pert, &
                               optical_props, source, &
                               t_lev, t_lev_pert_lay, t_lev_pert_layp1, &
                               top_at_1, ncol, n_quad_angles)

  type(ty_fluxes_broadband), intent(inout) :: fluxes
  real(wp), dimension(:,:), intent(inout), target :: flux_up, flux_dn
  real(wp), dimension(:, :), intent(inout) :: flux_up_surf_pert
  real(wp), dimension(:, :), intent(inout) :: flux_dn_surf_pert
  real(wp), dimension(:,:,:), intent(inout):: flux_up_gas_pert, flux_up_T_pert
  real(wp), dimension(:,:,:), intent(inout):: flux_dn_gas_pert, flux_dn_T_pert
  type(ty_gas_optics_rrtmgp), intent(inout) :: k_dist
  real(wp), dimension(:,:), intent(inout) :: sfc_emis_spec
  real(wp), intent(in) :: sfc_emis
  real(wp), dimension(:,:), intent(in) :: p_lay
  real(wp), dimension(:,:), intent(in) :: p_lev
  real(wp), dimension(:,:), intent(in) :: t_lay, t_lay_pert
  real(wp), dimension(:), intent(in) :: sfc_t, sfc_t_pert
  type(ty_gas_concs), intent(in) :: gas_conc, gas_conc_pert
  type(ty_optical_props_1scl), intent(inout) :: optical_props
  type(ty_source_func_lw), intent(inout) :: source
  real(wp), dimension(:,:), intent(in) :: t_lev, t_lev_pert_lay, t_lev_pert_layp1
  logical, intent(in) :: top_at_1
  integer, intent(in) :: ncol
  integer, intent(in) :: n_quad_angles
  type(ty_fluxes_broadband) :: fluxes_tmp  
  real(wp), dimension(size(flux_up,1),size(flux_up,2)), target :: flux_up_tmp, flux_dn_tmp                                                   
  type(ty_optical_props_1scl) :: optical_props_pert, optical_props_tmp        
  type(ty_source_func_lw) :: source_lay, source_layp1, source_pert                                                  
  type(ty_gas_concs) :: gas_conc_1d                                           
  integer :: nbnd
  integer :: icol, ilay
  integer :: ibnd

  call stop_on_err(source_lay%alloc(size(p_lay,1), size(p_lay,2), k_dist))
  call stop_on_err(source_layp1%alloc(size(p_lay,1), size(p_lay,2), k_dist))
  call stop_on_err(source_pert%alloc(size(p_lay,1), size(p_lay,2), k_dist))
  
  call stop_on_err(optical_props_pert%alloc_1scl(size(p_lay,1), size(p_lay,2), k_dist))
  call stop_on_err(optical_props_tmp%alloc_1scl(size(p_lay,1), size(p_lay,2), k_dist))

  nbnd = k_dist%get_nband()
  fluxes%flux_up => flux_up
  fluxes%flux_dn => flux_dn
  fluxes_tmp%flux_up => flux_up_tmp
  fluxes_tmp%flux_dn => flux_dn_tmp
    do icol = 1, ncol
      do ibnd = 1, nbnd
        sfc_emis_spec(ibnd,icol) = sfc_emis
      enddo
    enddo                 
                                         
  call stop_on_err(k_dist%gas_optics_int(p_lay, p_lev, t_lay, sfc_t, &
                                     gas_conc, optical_props, source, &
                                     tlev=t_lev))
    
  call stop_on_err(k_dist%gas_optics_int(p_lay, p_lev, t_lay_pert, sfc_t_pert, &
                                     gas_conc, optical_props, source_lay, &
                                     tlev=t_lev_pert_lay))
  call stop_on_err(k_dist%gas_optics_int(p_lay, p_lev, t_lay_pert, sfc_t_pert, &
                                     gas_conc_pert, optical_props_pert, source_layp1, &
                                     tlev=t_lev_pert_layp1))   

  call stop_on_err(rte_lw(optical_props, top_at_1, source, sfc_emis_spec, fluxes, &
                          n_gauss_angles=n_quad_angles))

optical_props_tmp%tau = optical_props%tau
    do ilay = 1, size(p_lay,2)
      optical_props_tmp%tau(:,ilay,:) = optical_props_pert%tau(:,ilay,:)
          call stop_on_err(rte_lw(optical_props_tmp, top_at_1, source, sfc_emis_spec, fluxes_tmp, &
                          n_gauss_angles=n_quad_angles))    
            flux_up_gas_pert(:,:,ilay)   = fluxes_tmp%flux_up(:,:)
            flux_dn_gas_pert(:,:,ilay)   = fluxes_tmp%flux_dn(:,:)
      optical_props_tmp%tau(:,ilay,:) = optical_props%tau(:,ilay,:)
    enddo      

source_pert%sfc_source = source%sfc_source
source_pert%lev_source_inc = source%lev_source_inc
source_pert%lev_source_dec = source%lev_source_dec
source_pert%lay_source     = source%lay_source  
      do ilay = 1, size(p_lay,2)  
        source_pert%lev_source_inc(:,ilay,:) = source_lay%lev_source_inc(:,ilay,:) 
        source_pert%lev_source_dec(:,ilay,:) = source_lay%lev_source_dec(:,ilay,:) 
        if (ilay .NE. size(p_lay,2)) then
            source_pert%lev_source_inc(:,ilay+1,:)= source_layp1%lev_source_inc(:,ilay,:)
            source_pert%lev_source_dec(:,ilay+1,:)= source_layp1%lev_source_dec(:,ilay,:)
        endif
        source_pert%lay_source(:,ilay,:)     = source_lay%lay_source(:,ilay,:)   

        call stop_on_err(rte_lw(optical_props, top_at_1, source_pert, sfc_emis_spec, fluxes_tmp, &
                          n_gauss_angles=n_quad_angles))   
        flux_up_T_pert(:,:, ilay)   = fluxes_tmp%flux_up(:,:)
        flux_dn_T_pert(:,:, ilay)   = fluxes_tmp%flux_dn(:,:)

        source_pert%lev_source_inc(:,ilay,:) = source%lev_source_inc(:,ilay,:)
        source_pert%lev_source_dec(:,ilay,:) = source%lev_source_dec(:,ilay,:)
        source_pert%lay_source(:,ilay,:)     = source%lay_source(:,ilay,:)  
      enddo

    source_pert%sfc_source = source_lay%sfc_source
    source_pert%lev_source_inc = source%lev_source_inc
    source_pert%lev_source_dec = source%lev_source_dec
    source_pert%lay_source     = source%lay_source       
    call stop_on_err(rte_lw(optical_props, top_at_1, source_pert, sfc_emis_spec, fluxes_tmp, &
                          n_gauss_angles=n_quad_angles))   
    flux_up_surf_pert   = fluxes_tmp%flux_up
    flux_dn_surf_pert   = fluxes_tmp%flux_dn                                                                                
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

