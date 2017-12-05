subroutine stop_on_err(error_msg)
  use iso_fortran_env, only : error_unit
  character(len=*), intent(in) :: error_msg

  if(error_msg /= "") then
    write (error_unit,*) trim(error_msg)
    write (error_unit,*) "test_flux_compute stopping"
    stop
  end if

end subroutine stop_on_err
program rrtmgp_rfmip_lw
  use mo_rte_kind,           only: wp
  use mo_gas_optics,         only: ty_gas_optics_specification
  use mo_gas_concentrations, only: ty_gas_concs
  use mo_optical_props,      only: ty_optical_props, ty_optical_props_1scl
  use mo_fluxes,             only: ty_fluxes
  ! ---- RRTMPG driver
  use mo_rte_lw,          only: rte_lw, rte_lw_init

  ! ---- I/O
  use mo_rfmip_io,           only: read_size, read_and_block_pt, read_and_block_gases_ty, unblock_and_write, &
                                   read_and_block_lw_bc, read_kdist_gas_names
  use mo_load_coefficients,  only: load_and_init
  implicit none
  ! --------------------------------------------------
  character(len=132)         :: fileName    = 'multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-0-4_none.nc'
  character(len=132)         :: k_dist_file = 'coefficients_lw.nc'
  logical                    :: top_at_1
  integer                    :: ncol, nlay, nexp, ngpt, nblocks, block_size = 4
  character(len=32 ), &
            dimension(:), &
                 allocatable :: kdist_gas_names, gases_to_use
  integer                    :: b
  real(wp), dimension(:,:,:), &
                 allocatable  :: p_lay, p_lev, t_lay, t_lev ! block_size, nlay, nblocks
  real(wp), dimension(:,:,:), target,  &
                 allocatable  :: flux_up, flux_dn
  real(wp), dimension(:,:  ), &
                 allocatable  :: sfc_emis, sfc_t            ! block_size, nblocks

  real(wp), dimension(:,:,:), &
                 allocatable  :: lay_src, lev_src_dec, lev_src_inc
  real(wp), dimension(:,:  ), &
                 allocatable  :: sfc_src

  type(ty_gas_concs), dimension(:), &
                       allocatable  :: gas_conc_array
  type(ty_gas_optics_specification) :: k_dist
  type(ty_optical_props_1scl)       :: optical_props
  type(ty_fluxes)                   :: fluxes
  ! --------------------------------------------------
  !
  ! Update file names, block size
  !
  ! How big is the problem? Does it fit into blocks of the size we've specified?
  !
  call read_size(fileName, ncol, nlay, nexp)
  if(mod(ncol*nexp, block_size) /= 0 ) call stop_on_err("rrtmgp_rfmip_lw: number of columns doesn't fit evenly into blocks.")
  nblocks = (ncol*nexp)/block_size

  !
  ! Names of gases known to the k-distribution - default is all.
  !
  call read_kdist_gas_names(k_dist_file, kdist_gas_names)
  !
  ! Here could provide variants i.e. using equivalent concentrations
  !
  gases_to_use = kdist_gas_names
  print *, "Radiation calculation uses gases "
  print *, "  ", [(trim(gases_to_use(b)) // " ", b = 1, size(gases_to_use))]

  !
  ! Allocation on assignment within reading routines
  !
  call read_and_block_pt(   fileName, block_size, p_lay, p_lev, t_lay, t_lev)
  top_at_1 = p_lay(1, 1, 1) < p_lay(1, nlay, 1)
  !
  ! Read the gas concentrations
  !
  call read_and_block_gases_ty(fileName, block_size, gases_to_use, gas_conc_array)

  call read_and_block_lw_bc(fileName, block_size, sfc_emis, sfc_t)

  ! Read k-distribution
  !
  call load_and_init(k_dist, trim(k_dist_file), gas_conc_array(1))
  if(.not. k_dist%is_internal_source_present()) &
    stop "rrtmgp_rfmip_lw: k-distribution file isn't LW"
  ngpt = k_dist%get_ngpt()

  !
  ! RRTMGP won't run with pressure less than its minimum. The top level in the RFMIP file
  !   is set to 10^-3 Pa. Here we pretend the layer is just a bit less deep.
  !
  if(top_at_1) then
    p_lev(:,1,:) = k_dist%get_press_ref_min() + epsilon(k_dist%get_press_ref_min())
  else
    p_lev(:,nlay+1,:) &
                 = k_dist%get_press_ref_min() + epsilon(k_dist%get_press_ref_min())
  end if

  allocate(flux_up(    block_size, nlay+1, nblocks), &
           flux_dn(    block_size, nlay+1, nblocks))
  allocate(lay_src(    block_size, nlay,   ngpt), &
           lev_src_inc(block_size, nlay,   ngpt), &
           lev_src_dec(block_size, nlay,   ngpt), &
           sfc_src(    block_size,         ngpt))
  call stop_on_err(optical_props%init_1scl(block_size, nlay, ngpt))
  !
  ! Loop over blocks -- use OpenMP?
  !   Would need private copies of source arrays, optical props, fluxes type
  !   Maybe the latter would need to be initialized inside the loop
  !
  do b = 1, nblocks
    call stop_on_err(k_dist%gas_optics(p_lay(:,:,b), &
                  									   p_lev(:,:,b),       &
                  									   t_lay(:,:,b),       &
                  									   sfc_t(:  ,b),       &
                  									   gas_conc_array(b),  &
                  									   optical_props,      &
                  									   lay_src,            &
                  									   lev_src_inc,        &
                  									   lev_src_dec,        &
                  									   sfc_src,            &
                  									   tlev = t_lev(:,:,b)))

    fluxes%flux_up => flux_up(:,:,b)
    fluxes%flux_dn => flux_dn(:,:,b)
	  call stop_on_err(rte_lw(optical_props,   &
            							  top_at_1,        &
            							  k_dist,          &
            							  lay_src,         &
            							  lev_src_inc,     &
            							  lev_src_dec,     &
            							  spread(sfc_emis(:,b), 1, ncopies = k_dist%get_nband()), &
            							  sfc_src,         &
            							  fluxes))
  end do

  call unblock_and_write('rlu_template.nc', 'rlu', flux_up)
  call unblock_and_write('rld_template.nc', 'rld', flux_dn)
end program rrtmgp_rfmip_lw
