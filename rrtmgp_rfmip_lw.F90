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
  use mo_rrtmgp_kind,        only: wp
  use mo_gas_optics_specification, &
                             only: ty_gas_optics_specification
  use mo_gas_concentrations, only: ty_gas_concs
  use mo_optical_props,      only: ty_optical_props, ty_optical_props_1scl
  use mo_fluxes,             only: ty_fluxes
  ! ---- RRTMPG driver
  use mo_rrtmgp_lw,          only: rrtmgp_lw, rrtmgp_lw_init
                        
  ! ---- I/O 
  use mo_rfmip_io,           only: read_size, read_and_block_pt, read_and_block_gases_ty, unblock_and_write, 
                                   read_and_block_lw_bc
  use mo_load_coefficients,  only: load_and_init
  implicit none 
  ! --------------------------------------------------
  character(len=132)         :: input_file  = 'multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-0-4_none.nc'
  character(len=132)         :: k_dist_file = 'coefficients_lw.nc'
  logical                    :: top_at_1
  integer                    :: ncol, nlay, nexp, ngpt, nblocks, block_size = 8 
  integer                    :: b 
  real(wp), dimension(:,:,:), & 
                 allocatable  :: p_lay, p_lev, t_lay, t_lev ! block_size, nlay, nblocks 
  real(wp), dimension(:,:,:), & 
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
  type(ty_optical_props_1scl), & 
                      dimension(1)  :: optical_props
  type(ty_fluxes),    dimension(1)  :: fluxes
  ! --------------------------------------------------
  ! 
  ! Update file names, block size 
  !
  
  call load_and_init(k_dist, trim(k_dist_file))
  if(k_dist%is_internal_source_present()) & 
    stop "rrtmgp_rfmip_lw: k-distribution file isn't LW"   
  ngpt = k_dist%get_ngpt() 
  
  ! 
  ! How big is the problem? Does it fit into blocks of the size we've specified? 
  ! 
  call read_size(fileName, ncol, nlay, nexp) 
  if(mod(ncol*nexp, block_size) /= 0 ) call stop_on_err("rrtmgp_rfmip_lw: number of columns doesn't fit evenly into blocks.") 
  nblocks = (ncol*nexp)/block_size
  
  allocate(flux_up(    block_size, nlay, nblocks), & 
           flux_dn(    block_size, nlay, nblocks))
  allocate(lay_src(    block_size, nlay, ngpt), & 
           lev_src_inc(block_size, nlay, ngpt), & 
           lev_src_dec(block_size, nlay, ngpt), & 
           sfc_src(    block_size,       ngpt))
  call stop_on_err(optical_props(1)%init_1scl(block_size, nlay, ngpt))
  
  !
  ! Allocation on assignment within reading routines 
  !
  call read_and_block_pt(   fileName, block_size, p_lay, p_lev, t_lay, t_lev)
  call read_and_block_lw_bc(fileName, block_size, sfc_emis, sfc_t) 
  !
  ! Need to provide variants i.e. using equivalent concentrations 
  !
  call read_and_block_gases(fileName, block_size, k_dist%get_gases, gas_conc_array)
  top_at_1 = p_lay(1, 1) < p_lay(1, nlay)
  
  ! 
  ! Loop over blocks -- use OpenMP? 
  !   Would need private copies of source arrays, optical props, fluxes type  
  !   Maybe the latter would need to be initialized inside the loop 
  !
  do b = 1, nblocks
    call stop_on_err(k_dist%gas_optics(p_lay(:,:,b)/100._wp, &
									   p_lev(:,:,b)/100._wp, &
									   t_lay(:,:,b),       & 
									   sfc_t(:,:,b),       & 
									   gas_conc_array(b),  &
									   optical_props(1),   &
									   lay_src,            & 
									   lev_src_inc,        & 
									   lev_src_dec,        &
									   sfc_src,            & 
									   tlev = t_lev(:,:,b)))

    fluxes(1)%flux_up => flux_up(:,:,b) 
    fluxes(1)%flux_dn => flux_dn(:,:,b) 
    
	call stop_on_err(rrtmgp_lw(optical_props,   &
							   top_at_1,        &
							   k_dist,          &
							   lay_src,         & 
							   lev_src_inc,     & 
							   lev_src_dec,     &
							   sfc_emis(:,:,b), &
							   sfc_src,         & 
							   fluxes))
  end do 
  
  call unblock_and_write(fileName, 'rlu', flux_up)                      
  call unblock_and_write(fileName, 'rld', flux_dn)                      
end program rrtmgp_rfmip_lw