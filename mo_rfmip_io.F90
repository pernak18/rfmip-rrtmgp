module mo_rfmip_io
  use mo_rte_kind,   only: wp
  use mo_gas_concentrations, &
                        only: ty_gas_concs
  use mo_util_string,   only: lower_case, string_in_array, string_loc_in_array
  use netcdf
  implicit none
  private
  public :: read_kdist_gas_names, read_size, read_and_block_pt, &
            read_and_block_sw_bc, read_and_block_lw_bc, read_and_block_gases_ty, &
            unblock_and_write

  interface read_field
    module procedure read_scalar, read_1d_field, read_2d_field, read_3d_field, read_4d_field
  end interface

  integer :: ncol_l = 0, nlay_l = 0, nexp_l = 0 ! Local copies

! Return T, p, surface albedos, temperatures,
! Gas concentrations can be scalars except where experiment boundaries are crossed
! Options: all gases, (CO2, CH4, N2O) + {CFC11eq; CFC12eq + HFC-134eq}. Always include ozone
!   What about other gases eg. CO that aren't part of the specification?

! Does output come blocked? Yes -- so block size must be known ahead of time.

! Routines --
! read_and_block_gases
! unblock_and_write flux

! Pack so that all mus > 0 are contiguous (offline, have an upacking routine too)
contains
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! This routine reads the names of the gases known to the k-distribution 
  !   
  !
  subroutine read_kdist_gas_names(fileName, kdist_gas_names)
    character(len=*),          intent(in   ) :: fileName
    character(len=32), dimension(:), allocatable, & 
                               intent(  out) :: kdist_gas_names
    ! ---------------------------
    integer :: ncid, varid
    character(len=8), parameter :: varName = "gas_names" 
    ! ---------------------------
    if(nf90_open(trim(fileName), NF90_NOWRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("read_kdist_gas_names: can't find file " // trim(fileName))

    allocate(kdist_gas_names(get_dim_length(ncid, 'absorber')))
    
    if(nf90_inq_varid(ncid, trim(varName), varid) /= NF90_NOERR) &
      call stop_on_err("read_kdist_gas_names: can't find variable" // trim(varName))
    if(nf90_get_var(ncid, varid, kdist_gas_names)  /= NF90_NOERR) &
      call stop_on_err("read_kdist_gas_names: can't read variable" // trim(varName))
    
    ncid = nf90_close(ncid)
  end subroutine read_kdist_gas_names
  !--------------------------------------------------------------------------------------------------------------------
  subroutine read_size(fileName, ncol, nlay, nexp)
    character(len=*),          intent(in   ) :: fileName
    integer,         optional, intent(  out) :: ncol, nlay, nexp
    ! ---------------------------
    integer :: ncid
    ! ---------------------------
    if(nf90_open(trim(fileName), NF90_NOWRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("read_atmos: can't find file " // trim(fileName))

    ncol = get_dim_length(ncid, 'site')
    nlay = get_dim_length(ncid, 'layer')
    nexp = get_dim_length(ncid, 'expt')
    if(get_dim_length(ncid, 'level') /= nlay+1) call stop_on_err("read_size: number of levels should be nlay+1")
    ncid = nf90_close(ncid)

    ncol_l = ncol
    nlay_l = nlay
    nexp_l = nexp
  end subroutine read_size
  !--------------------------------------------------------------------------------------------------------------------
  subroutine read_and_block_pt(fileName, blocksize, &
                               p_lay, p_lev, t_lay, t_lev)
    character(len=*),           intent(in   ) :: fileName
    integer,                    intent(in   ) :: blocksize
    real(wp), dimension(:,:,:), allocatable, & ! [blocksize, nlay/+1, nblocks]
                                intent(  out) :: p_lay, p_lev, t_lay, t_lev
    ! ---------------------------
    integer :: ncid
    integer :: nblocks
    ! ---------------------------
    if(any([ncol_l, nlay_l, nexp_l]  == 0)) call stop_on_err("read_and_block_pt: Haven't read problem size yet.")
    if(mod(ncol_l*nexp_l, blocksize) /= 0 ) call stop_on_err("read_and_block_pt: number of columns doesn't fit evenly into blocks.")
    nblocks = (ncol_l*nexp_l)/blocksize
    !
    ! Check that output arrays are sized correctly : blocksize, nlay, (ncol * nexp)/blocksize
    !

    if(nf90_open(trim(fileName), NF90_NOWRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("read_and_block_pt: can't find file " // trim(fileName))

    p_lay = reshape(read_field(ncid, "press_layer", ncol_l, nlay_l,   nexp_l), &
                    shape = [blocksize, nblocks, nlay_l], order = [1, 3, 2])
    t_lay = reshape(read_field(ncid, "temp_layer",  ncol_l, nlay_l,   nexp_l), &
                    shape = [blocksize, nblocks, nlay_l], order = [1, 3, 2])
    p_lay = reshape(read_field(ncid, "press_level", ncol_l, nlay_l+1, nexp_l), &
                    shape = [blocksize, nblocks, nlay_l], order = [1, 3, 2])
    t_lay = reshape(read_field(ncid, "temp_level",  ncol_l, nlay_l+1, nexp_l), &
                    shape = [blocksize, nblocks, nlay_l], order = [1, 3, 2])

    ncid = nf90_close(ncid)
  end subroutine read_and_block_pt
  !--------------------------------------------------------------------------------------------------------------------
  subroutine read_and_block_sw_bc(fileName, blocksize, &
                               surface_albedo, total_solar_irradiance, solar_zenith_angle)
    character(len=*),           intent(in   ) :: fileName
    integer,                    intent(in   ) :: blocksize
    real(wp), dimension(:,:), allocatable, &
                                intent(  out) :: surface_albedo, total_solar_irradiance, solar_zenith_angle
    ! ---------------------------
    integer :: ncid
    integer :: nblocks
    ! ---------------------------
    if(any([ncol_l, nlay_l, nexp_l]  == 0)) call stop_on_err("read_and_block_sw_bc: Haven't read problem size yet.")
    if(mod(ncol_l*nexp_l, blocksize) /= 0 ) call stop_on_err("read_and_block_sw_bc: number of columns doesn't fit evenly into blocks.")
    nblocks = (ncol_l*nexp_l)/blocksize
    !
    ! Check that output arrays are sized correctly : blocksize, nlay, (ncol * nexp)/blocksize
    !

    if(nf90_open(trim(fileName), NF90_NOWRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("read_and_block_sw_bc: can't find file " // trim(fileName))

    surface_albedo         = reshape(spread(read_field(ncid, "surface_albedo",         ncol_l), dim=1, ncopies=nexp_l), &
                                     shape = [blocksize, nblocks])
    total_solar_irradiance = reshape(spread(read_field(ncid, "total_solar_irradiance", ncol_l), dim=1, ncopies=nexp_l), &
                                    shape = [blocksize, nblocks])
    solar_zenith_angle     = reshape(spread(read_field(ncid, "solar_zenith_angle",     ncol_l), dim=1, ncopies=nexp_l), &
                                     shape = [blocksize, nblocks])

    ncid = nf90_close(ncid)
  end subroutine read_and_block_sw_bc
  !--------------------------------------------------------------------------------------------------------------------
  subroutine read_and_block_lw_bc(fileName, blocksize, &
                                  surface_emissivity, surface_temperature)
    character(len=*),           intent(in   ) :: fileName
    integer,                    intent(in   ) :: blocksize
    real(wp), dimension(:,:), allocatable, &
                                intent(  out) :: surface_emissivity, surface_temperature
    ! ---------------------------
    integer :: ncid
    integer :: nblocks
    ! ---------------------------
    if(any([ncol_l, nlay_l, nexp_l]  == 0)) &
      call stop_on_err("read_and_block_lw_bc: Haven't read problem size yet.")
    if(mod(ncol_l*nexp_l, blocksize) /= 0 ) &
      call stop_on_err("read_and_block_lw_bc: number of columns doesn't fit evenly into blocks.")
    nblocks = (ncol_l*nexp_l)/blocksize
    !
    ! Check that output arrays are sized correctly : blocksize, nlay, (ncol * nexp)/blocksize
    !

    if(nf90_open(trim(fileName), NF90_NOWRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("read_and_block_lw_bc: can't find file " // trim(fileName))

    surface_emissivity  = reshape(spread(read_field(ncid, "surface_emissivity",  ncol_l), dim=1, ncopies=nexp_l), &
                                  shape = [blocksize, nblocks])
    surface_temperature = reshape(spread(read_field(ncid, "surface_temperature", ncol_l), dim=1, ncopies=nexp_l), &
                                  shape = [blocksize, nblocks])

    ncid = nf90_close(ncid)
  end subroutine read_and_block_lw_bc
  !--------------------------------------------------------------------------------------------------------------------
  subroutine read_and_block_gases_ty(fileName, blocksize, gas_names, gas_conc_array)
    character(len=*),           intent(in   ) :: fileName
    integer,                    intent(in   ) :: blocksize
    character(len=*),  dimension(:), &
                                intent(in   ) :: gas_names ! Provided by gas_optics -- which gases do we want to read
    type(ty_gas_concs), dimension(:), allocatable, &
                                intent(  out) :: gas_conc_array

    ! ---------------------------
    integer :: ncid
    integer :: nblocks
    integer :: b, g
    integer,  dimension(:,:),   allocatable :: exp_num
    real(wp), dimension(:),     allocatable :: gas_conc_temp_1d
    real(wp), dimension(:,:,:), allocatable :: gas_conc_temp_3d
    character(len=32)                       :: gas_name_in_file
    character(len=32), dimension(10) :: &
      chem_name = ['co   ', &
                   'ch4  ', &
          			   'o2   ', &
          			   'n2o  ', &
          			   'n2   ', &
          			   'co2  ', &
          			   'CCl4 ', &
          			   'ch4  ', &
          			   'CH3Br', &
			   'CH3Cl'], &
      desc_name = ['carbon_monoxide     ', &
                   'methane             ', &
                   'oxygen              ', &
      			       'nitrous_oxide       ', &
      			       'nitrogen            ', &
      			       'carbon_dioxide      ', &
        				   'carbon_tetrachloride', &
        				   'methane             ', &
        				   'methyl_bromide      ', &
        				   'methyl_chloride     ']
    ! ---------------------------
    if(any([ncol_l, nlay_l, nexp_l]  == 0)) &
      call stop_on_err("read_and_block_lw_bc: Haven't read problem size yet.")
    if(mod(ncol_l*nexp_l, blocksize) /= 0 ) &
      call stop_on_err("read_and_block_lw_bc: number of columns doesn't fit evenly into blocks.")
    nblocks = (ncol_l*nexp_l)/blocksize
    ! Experiment index for each colum
    exp_num = reshape(spread([(b, b = 1, nexp_l)], 1, ncopies = ncol_l), shape = [blocksize, nblocks])

    if(nf90_open(trim(fileName), NF90_NOWRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("read_and_block_gases_ty: can't find file " // trim(fileName))


    allocate(gas_conc_array(nblocks))

    !
    ! Water vapor and ozone depend on col, lay, exp: look just like other fields
    !
    gas_conc_temp_3d = reshape(read_field(ncid, "water_vapor", ncol_l, nlay_l, nexp_l), &
                               shape = [blocksize, nblocks, nlay_l], order = [1, 3, 2])
    do b = 1, nblocks
      call stop_on_err(gas_conc_array(b)%set_vmr('h2o', gas_conc_temp_3d(:,:,b)))
    end do
    gas_conc_temp_3d = reshape(read_field(ncid, "ozone", ncol_l, nlay_l, nexp_l), &
                               shape = [blocksize, nblocks, nlay_l], order = [1, 3, 2])
    do b = 1, nblocks
      call stop_on_err(gas_conc_array(b)%set_vmr('o3', gas_conc_temp_3d(:,:,b)))
    end do

    !
    ! All other gases are a function of experiment only
    !
    do g = 1, size(gas_names)
      gas_name_in_file = trim(lower_case(gas_names(g)))
      if(gas_name_in_file == 'h2o' .or. gas_name_in_file == 'o3') cycle
      !
      ! Use a mapping between chemical formula and name if it exists
      !
      if(string_in_array(gas_name_in_file, chem_name)) &
        gas_name_in_file = desc_name(string_loc_in_array(gas_name_in_file, chem_name))
      gas_name_in_file = gas_name_in_file // "_GM"

      ! Read the values as a function of experiment
      gas_conc_temp_1d = read_field(ncid, gas_name_in_file, nexp_l)

	  do b = 1, nblocks
        ! Does every value in this block belong to the same experiment?
	    if(all(exp_num(1,b) == exp_num(2:,b))) then
	      ! Provide a scalar value
		    call stop_on_err(gas_conc_array(b)%set_vmr(gas_names(g), gas_conc_temp_1d(exp_num(1,b))))
		  else
		  ! Create 2D field, blocksize x nlay, with scalar values from each experiment
		  call stop_on_err(gas_conc_array(b)%set_vmr(gas_names(g), &
		                                             spread(gas_conc_temp_1d(exp_num(:,b)), 2, ncopies = nlay_l)))
		  end if
	  end do

    end do
    ncid = nf90_close(ncid)
  end subroutine read_and_block_gases_ty
  !--------------------------------------------------------------------------------------------------------------------
  subroutine unblock_and_write(fileName, varName, values)
    character(len=*),           intent(in   ) :: fileName, varName
    real(wp), dimension(:,:,:),  & ! [blocksize, nlay/+1, nblocks]
                                intent(in   ) :: values
    ! ---------------------------
    integer :: ncid
    integer :: nblocks
    ! ---------------------------
    if(any([ncol_l, nlay_l, nexp_l]  == 0)) call stop_on_err("read_and_block_pt: Haven't read problem size yet.")
    !
    ! Check that output arrays are sized correctly : blocksize, nlay, (ncol * nexp)/blocksize
    !

    if(nf90_open(trim(fileName), NF90_WRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("unblock_and_write: can't find file " // trim(fileName))

    call stop_on_err(write_3d_field(ncid, varName,                      &
                     reshape(reshape(values,                            &
                                     shape = [ncol_l, nexp_l, nlay_l]), &
                             shape = [ncol_l, nlay_l, nexp_l], order = [1, 3, 2])))

    ncid = nf90_close(ncid)
  end subroutine unblock_and_write
  !
  !
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Ancillary functions
  !
  !--------------------------------------------------------------------------------------------------------------------
  function read_scalar(ncid, varName)
    integer,          intent(in) :: ncid
    character(len=*), intent(in) :: varName
    real(wp)                     :: read_scalar

    integer :: varid

    if(nf90_inq_varid(ncid, trim(varName), varid) /= NF90_NOERR) &
      call stop_on_err("read_field: can't find variable " // trim(varName))
    if(nf90_get_var(ncid, varid, read_scalar)  /= NF90_NOERR) &
      call stop_on_err("read_field: can't read variable " // trim(varName))

  end function read_scalar
  !--------------------------------------------------------------------------------------------------------------------
  function read_1d_field(ncid, varName, nx)
    integer,          intent(in) :: ncid
    character(len=*), intent(in) :: varName
    integer,          intent(in) :: nx
    real(wp), dimension(nx)      :: read_1d_field

    integer :: varid

    if(nf90_inq_varid(ncid, trim(varName), varid) /= NF90_NOERR) &
      call stop_on_err("read_field: can't find variable " // trim(varName))
    if(nf90_get_var(ncid, varid, read_1d_field)  /= NF90_NOERR) &
      call stop_on_err("read_field: can't read variable " // trim(varName))

  end function read_1d_field
  !--------------------------------------------------------------------------------------------------------------------
  function read_2d_field(ncid, varName, nx, ny)
    integer,          intent(in) :: ncid
    character(len=*), intent(in) :: varName
    integer,          intent(in) :: nx, ny
    real(wp), dimension(nx, ny)  :: read_2d_field

    integer :: varid

    if(nf90_inq_varid(ncid, trim(varName), varid) /= NF90_NOERR) &
      call stop_on_err("read_field: can't find variable " // trim(varName))
    if(nf90_get_var(ncid, varid, read_2d_field)  /= NF90_NOERR) &
      call stop_on_err("read_field: can't read variable " // trim(varName))

  end function read_2d_field
  !--------------------------------------------------------------------------------------------------------------------
  function read_3d_field(ncid, varName, nx, ny, nz)
    integer,          intent(in) :: ncid
    character(len=*), intent(in) :: varName
    integer,          intent(in) :: nx, ny, nz
    real(wp), dimension(nx, ny, nz)  :: read_3d_field

    integer :: varid

    if(nf90_inq_varid(ncid, trim(varName), varid) /= NF90_NOERR) &
      call stop_on_err("read_field: can't find variable " // trim(varName))
    if(nf90_get_var(ncid, varid, read_3d_field)  /= NF90_NOERR) &
      call stop_on_err("read_field: can't read variable " // trim(varName))

  end function read_3d_field
  !--------------------------------------------------------------------------------------------------------------------
  function read_4d_field(ncid, varName, nw, nx, ny, nz)
    integer,          intent(in) :: ncid
    character(len=*), intent(in) :: varName
    integer,          intent(in) :: nw, nx, ny, nz
    real(wp), dimension(nw, nx, ny, nz)  :: read_4d_field

    integer :: varid

    if(nf90_inq_varid(ncid, trim(varName), varid) /= NF90_NOERR) &
      call stop_on_err("read_field: can't find variable " // trim(varName))
    if(nf90_get_var(ncid, varid, read_4d_field)  /= NF90_NOERR) &
      call stop_on_err("read_field: can't read variable " // trim(varName))

  end function read_4d_field
  !--------------------------------------------------------------------------------------------------------------------
  function get_dim_length(ncid, dimname)
    !
    ! Get the length of a dimension from an open netCDF file
    !  This is unfortunately a two-step process
    !
    integer,          intent(in) :: ncid
    character(len=*), intent(in) :: dimname
    integer :: get_dim_length

    integer :: dimid

    if(nf90_inq_dimid(ncid, trim(dimname), dimid) == NF90_NOERR) then
      if(nf90_inquire_dimension(ncid, dimid, len=get_dim_length) /= NF90_NOERR) get_dim_length = 0
    else
      get_dim_length = 0
    end if

  end function get_dim_length
  !--------------------------------------------------------------------------------------------------------------------
  function var_exists(ncid, varName)
    !
    ! Does this variable exist (have a valid var_id) in the open netCDF file?
    !
    integer,          intent(in) :: ncid
    character(len=*), intent(in) :: varName
    logical :: var_exists

    integer :: varId
    var_exists = nf90_inq_varid(ncid, trim(varName), varid) == NF90_NOERR
  end function var_exists
  !--------------------------------------------------------------------------------------------------------------------
  function dim_exists(ncid, dimName)
    !
    ! Does this dimension exist (have a valid dim_id) in the open netCDF file?
    !
    integer,          intent(in) :: ncid
    character(len=*), intent(in) :: dimName
    logical :: dim_exists

    integer :: dimid
    dim_exists = nf90_inq_dimid(ncid, trim(dimName), dimid) == NF90_NOERR
  end function dim_exists
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
  function write_3d_field(ncid, varName, var) result(err_msg)
    integer,                    intent(in) :: ncid
    character(len=*),           intent(in) :: varName
    real(wp), dimension(:,:,:), intent(in) :: var
    character(len=128)                     :: err_msg

    integer :: varid

    err_msg = ""
    if(nf90_inq_varid(ncid, trim(varName), varid) /= NF90_NOERR) then
      err_msg = "write_field: can't find variable " // trim(varName)
      return
    end if
    if(nf90_put_var(ncid, varid, var)  /= NF90_NOERR) &
      err_msg = "write_field: can't write variable " // trim(varName)

  end function write_3d_field

end module mo_rfmip_io
