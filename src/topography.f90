module topography
  use iso_fortran_env
  use netcdf
  use utils
  implicit none

  type topography_t
    ! Dimensions
    integer(int32) :: nxt = 0
    integer(int32) :: nyt = 0
    ! Depth variable and attributes
    real(real32), allocatable :: depth(:,:)
    character(len=3) :: lakes_removed = "no "
    ! Fraction of cell covered by water variable
    real(real32), allocatable :: frac(:,:)
    ! Global attributes
    character(len=:), allocatable :: original_file
  contains
    procedure :: write => topography_write
    procedure :: copy => topography_copy
    generic   :: assignment(=) => copy
  end type topography_t

  interface topography_t
    module procedure topography_constructor
  end interface topography_t

contains

  type(topography_t) function topography_constructor(filename) result(topog)
    character(len=*), intent(in) :: filename

    integer(int32) :: ncid, depth_id, frac_id, dids(2) ! NetCDF ids

    write(*,*) "Reading topography from file '", trim(filename), "'"

    topog%original_file = filename

    ! Open file
    call handle_error(nf90_open(trim(filename), nf90_nowrite, ncid))

    ! Get dimensions
    call handle_error(nf90_inq_dimid(ncid, 'nx', dids(1)))
    call handle_error(nf90_inq_dimid(ncid, 'ny', dids(2)))
    call handle_error(nf90_inquire_dimension(ncid, dids(1), len=topog%nxt))
    call handle_error(nf90_inquire_dimension(ncid, dids(2), len=topog%nyt))

    ! Get depth
    allocate(topog%depth(topog%nxt, topog%nyt))
    call handle_error(nf90_inq_varid(ncid, 'depth', depth_id))
    call handle_error(nf90_get_var(ncid, depth_id, topog%depth))
    call handle_error(nf90_get_att(ncid, depth_id, 'lakes_removed', topog%lakes_removed), isfatal=.false.)

    ! Get fraction of water
    call handle_error(nf90_inq_varid(ncid, 'frac', frac_id))
    allocate(topog%frac(topog%nxt, topog%nyt))
    call handle_error(nf90_get_var(ncid, frac_id, topog%frac))

    ! Close file
    call handle_error(nf90_close(ncid))

  end function topography_constructor

  subroutine topography_copy(topog_out, topog_in)
    class(topography_t), intent(out) :: topog_out
    class(topography_t), intent(in)  :: topog_in

    ! Dimensions
    topog_out%nxt = topog_in%nxt
    topog_out%nyt = topog_in%nyt

    ! Depth variable and attributes
    allocate(topog_out%depth, source=topog_in%depth)
    topog_out%lakes_removed = topog_in%lakes_removed

    ! Fraction of cell covered by water variable
    allocate(topog_out%frac, source=topog_in%frac)

    ! Global attributes
    topog_out%original_file = topog_in%original_file

  end subroutine topography_copy

  subroutine topography_write(this, filename)
    class(topography_t), intent(in) :: this
    character(len=*), intent(in) :: filename

    integer(int32) :: ncid, depth_id, frac_id, dids(2) ! NetCDF ids

    write(*,*) "Writing topography to file '", trim(filename), "'"

    ! Open file
    call handle_error(nf90_create(trim(filename), ior(nf90_netcdf4, nf90_clobber), ncid))

    ! Write dimensions
    call handle_error(nf90_def_dim(ncid, 'nx', this%nxt, dids(1)))
    call handle_error(nf90_def_dim(ncid, 'ny', this%nyt, dids(2)))

    ! Write depth
    call handle_error(nf90_def_var(ncid, 'depth', nf90_float, dids, depth_id, chunksizes=[this%nxt/10, this%nyt/10], &
      deflate_level=1, shuffle=.true.))
    call handle_error(nf90_def_var_fill(ncid, depth_id, 0, MISSING_VALUE))
    call handle_error(nf90_put_att(ncid, depth_id, 'long_name', 'depth'))
    call handle_error(nf90_put_att(ncid, depth_id, 'units', 'm'))
    call handle_error(nf90_put_att(ncid, depth_id, 'lakes_removed', this%lakes_removed))
    call handle_error(nf90_put_var(ncid, depth_id, this%depth))

    ! Write frac
    call handle_error(nf90_def_var(ncid, 'frac', nf90_float, dids, frac_id, chunksizes=[this%nxt/10, this%nyt/10], &
      deflate_level=1, shuffle=.true.))
    call handle_error(nf90_def_var_fill(ncid, frac_id, 0, MISSING_VALUE))
    call handle_error(nf90_put_var(ncid, frac_id, this%frac))

    ! Write global attributes
    call handle_error(nf90_put_att(ncid, nf90_global, 'original_file', trim(this%original_file)))

    ! Close file
    call handle_error(nf90_enddef(ncid))

  end subroutine topography_write

end module topography
