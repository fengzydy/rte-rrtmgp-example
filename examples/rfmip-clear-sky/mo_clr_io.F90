! This code is part of Radiative Transfer for Energetics (RTE) and
!   RRTM for GCM Applications - Parallel (RRTMGP)
!
! Contacts: Robert Pincus and Eli Mlawer
! email:  rrtmgp@aer.com
!
! Copyright 2015-2018,  Atmospheric and Environmental Research and
! Regents of the University of Colorado.  All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
! -------------------------------------------------------------------------------------------------
!
! This module reads and writes to netCDF files (nominally for the Garand atmospheres) using
!   serial I/O. The files follow arbitrary conventions adopted by RTE+RRTMGP developers.
!   It may be useful as an example for other formats/conventions.
!
! Reading routines use "allocation on assignment," a feature of Fortran 2003 that may require
!   particular compilation flags.
!
module mo_clr_io
  !
  ! RTE+RRTMGP modules
  !
  use mo_rte_kind,           only: wp
  use mo_gas_concentrations, only: ty_gas_concs
  use mo_optical_props,      only: ty_optical_props
  !
  ! NetCDF I/O routines, shared with other RTE+RRTMGP examples
  !
  use mo_simple_netcdf,      only: read_field, read_string, var_exists, get_dim_size, &
                                   write_field, create_dim, create_var

  use netcdf
  implicit none
  private

  public :: read_atmos, write_lw_fluxes, write_sw_fluxes
contains
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Read profiles for all columns  -- T, p, and gas concentrations
  !   Allocation occurs on assignments (says the F2003 standard)
  !
  subroutine read_atmos(fileName,  ncol, nlay, nlev,    &
                        p_lay, t_lay, p_lev, t_lev, sfc_t, gas_concs, col_dry)
    character(len=*),   intent(in   ) :: fileName
    real(wp), dimension(:,:), allocatable,                 &
                        intent(inout) :: p_lay, t_lay, p_lev, t_lev, col_dry
                        type(ty_gas_concs), intent(inout) :: gas_concs
    integer, intent(inout) :: ncol, nlay, nlev
    real(wp), dimension(:), allocatable, intent(inout) :: sfc_t
    ! -------------------
    real(wp), dimension(:,:), allocatable:: gas_temp
    integer :: igas
    integer :: ncid, ndims, varid
    integer, parameter :: ngas = 8
    character(len=3), dimension(ngas) &
                       :: gas_names = ['h2o', 'co2', 'o3 ', 'n2o', 'co ', 'ch4', 'o2 ', 'n2 ']
    character(len=7) :: vmr_name
    ! -------------------

    if(nf90_open(trim(fileName), NF90_NOWRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("read_atmos: can't find file " // trim(fileName))

    ncol = get_dim_size(ncid, 'col')
    nlay = get_dim_size(ncid, 'lay')
    nlev = get_dim_size(ncid, 'lev')
    if(nlev /= nlay+1) call stop_on_err("read_atmos: nlev should be nlay+1")

    !
    ! These lines assume that compilers follow the Fortran 2003 standard for
    !   allocating on assignment. This may require explicit compiler support
    !   e.g. -assume realloc_lhs flag for Intel
    !
    sfc_t = read_field(ncid, 'sfc_t', ncol)
    p_lay = read_field(ncid, 'p_lay', ncol, nlay)
    t_lay = read_field(ncid, 't_lay', ncol, nlay)
    p_lev = read_field(ncid, 'p_lev', ncol, nlev)
    t_lev = read_field(ncid, 't_lev', ncol, nlev)

    allocate(gas_temp(ncol,nlay))

    call stop_on_err(gas_concs%init(gas_names))

    do igas = 1, ngas
      vmr_name = 'vmr_' // trim(gas_names(igas))
      if (.not. var_exists(ncid, trim(vmr_name))) then
        call stop_on_err("read_atmos: can't read concentration of " // trim(gas_names(igas)))
      endif
      if (nf90_inq_varid(ncid, trim(vmr_name), varid) .ne. nf90_noerr) then
        call stop_on_err("read_atmos: error trying to get variable id.")
      endif
      if (nf90_inquire_variable(ncid, varid, ndims=ndims) .ne. nf90_noerr) then
        call stop_on_err("read_atmos: error trying to get number of dimensions.")
      endif
      if (ndims .eq. 1) then
      gas_temp(:,1) = read_field(ncid, trim(vmr_name), ncol)
        call stop_on_err(gas_concs%set_vmr(trim(gas_names(igas)), gas_temp(1,1)))
      else
        gas_temp = read_field(ncid, trim(vmr_name), ncol, nlay)
        call stop_on_err(gas_concs%set_vmr(trim(gas_names(igas)), gas_temp))
      endif
    end do

    ! col_dry has unchanged allocation status on return if the variable isn't present in the netCDF file
    if(var_exists(ncid, 'col_dry')) col_dry = read_field(ncid, 'col_dry', ncol, nlay)

    ncid = nf90_close(ncid)
    deallocate(gas_temp)
  end subroutine read_atmos
  !--------------------------------------------------------------------------------------------------------------------
  subroutine write_sw_fluxes(fileName, flux_up, flux_dn, flux_dir)
    character(len=*),         intent(in) :: fileName
    real(wp), dimension(:,:), intent(in) :: flux_up, flux_dn, flux_dir
    ! -------------------
    integer :: ncid, ncol, nlev
    ! -------------------
    if(nf90_open(trim(fileName), NF90_WRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("write_fluxes: can't open file " // trim(fileName))
    !
    ! At present these dimension sizes aren't used
    !   We could certainly check the array sizes against these dimension sizes
    !
    ncol  = size(flux_up, dim=1)
    nlev  = get_dim_size(ncid, 'lev')
    call create_dim(ncid, "col_flx", ncol)

    call create_var(ncid, "sw_flux_up",  ["col_flx",  "lev    "], [ncol, nlev])
    call create_var(ncid, "sw_flux_dn",  ["col_flx",  "lev    "], [ncol, nlev])
    call create_var(ncid, "sw_flux_dir", ["col_flx",  "lev    "], [ncol, nlev])

    call stop_on_err(write_field(ncid, "sw_flux_up",  flux_up ))
    call stop_on_err(write_field(ncid, "sw_flux_dn",  flux_dn ))
    call stop_on_err(write_field(ncid, "sw_flux_dir", flux_dir))

    ncid = nf90_close(ncid)
  end subroutine write_sw_fluxes
  !--------------------------------------------------------------------------------------------------------------------
  subroutine write_lw_fluxes(fileName, flux_up, flux_dn, heating_rate, &
                     flux_up_trop_off, &
                     flux_up_strat_off, &
                     flux_up_air_off)
    character(len=*),         intent(in) :: fileName
    real(wp), dimension(:,:), intent(in) :: flux_up, flux_dn, flux_up_trop_off, flux_up_strat_off, flux_up_air_off, heating_rate
    ! -------------------
    integer :: ncid, ncol, nlev, nlay
    ! -------------------
    if(nf90_open(trim(fileName), NF90_WRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("write_fluxes: can't open file " // trim(fileName))

    ncol  = size(flux_up, dim=1)
    nlev  = get_dim_size(ncid, 'lev')
    nlay  = nlev - 1
    call create_dim(ncid, "col_flx", ncol)
    call create_dim(ncid, "lay    ", nlay)

    call create_var(ncid, "lw_flux_up",  ["col_flx",  "lev    "], [ncol, nlev])
    call create_var(ncid, "lw_flux_up_trop_off",  ["col_flx",  "lev    "], [ncol, nlev])
    call create_var(ncid, "lw_flux_up_strat_off",  ["col_flx",  "lev    "], [ncol, nlev])
    call create_var(ncid, "lw_flux_up_air_off",  ["col_flx",  "lev    "], [ncol, nlev])
    call create_var(ncid, "lw_flux_dn",  ["col_flx",  "lev    "], [ncol, nlev])
    call create_var(ncid, "lw_heating_rate", ["col_flx",  "lay    "], [ncol, nlay])

    call stop_on_err(write_field(ncid, "lw_flux_up",  flux_up ))
    call stop_on_err(write_field(ncid, "lw_flux_up_trop_off",  flux_up_trop_off ))
    call stop_on_err(write_field(ncid, "lw_flux_up_strat_off",  flux_up_strat_off ))
    call stop_on_err(write_field(ncid, "lw_flux_up_air_off",  flux_up_air_off ))
    call stop_on_err(write_field(ncid, "lw_flux_dn",  flux_dn ))
    call stop_on_err(write_field(ncid, "lw_heating_rate",  heating_rate ))

    ncid = nf90_close(ncid)
  end subroutine write_lw_fluxes
  !--------------------------------------------------------------------------------------------------------------------
  subroutine stop_on_err(msg)
    !
    ! Print error message and stop
    !
    use iso_fortran_env, only : error_unit
        character(len=*), intent(in) :: msg
    if(len_trim(msg) > 0) then
      write(error_unit,*) trim(msg)
      stop
    end if
  end subroutine
  !--------------------------------------------------------------------------------------------------------------------

subroutine reorder123x312(array, array_out)
  real(kind=wp), dimension(:,:,:), intent(in) :: array
  real(kind=wp), dimension(:,:,:), intent(out) :: array_out

  integer :: i, j, k
  do j = 1, size(array, 2)
    do i = 1, size(array, 1)
      do k = 1, size(array, 3)
        array_out(k, i, j) = array(i, j, k)
      enddo
    enddo
  enddo
end subroutine reorder123x312


end module mo_clr_io