program topogtools
  use, intrinsic :: iso_fortran_env
  use M_CLI2
  use check_nonadvective_m
  use deseas_m
  use fix_nonadvective_m
  use gen_topo_m
  use min_max_depth_m
  use topog2mask_m
  implicit none

  character(len=:), allocatable :: name
  character(len=:), allocatable :: help_text(:)
  character(len=:), allocatable :: file_in, file_out, hgrid, vgrid
  integer :: ii

  ! Read command line
  name = get_subcommand()
  select case (name)
  case ('gen_topo')
    call set_args('--input:i "unset" --output:o "unset" --hgrid:h "unset" --tripolar:t F --longitude-offset 0.0')

  case ('deseas')
    call set_args('--input:i "unset" --output:o "unset"')

  case ('min_max_depth')
    call set_args('--input:i "unset" --output:o "unset" --vgrid:v "ocean_vgrid.nc" --level:l 0')

  case ('fix_nonadvective')
    call set_args('--input:i "unset" --output:o "unset" --vgrid:v "ocean_vgrid.nc"')

  case ('check_nonadvective')
    call set_args('--input:i "unset" --vgrid:v "ocean_vgrid.nc"')

  case ('topog2mask')
    call set_args('--input:i "unset" --output:o "unset" --fraction:f 0.0')

  case ('')
    ! general help for "topogtools"
    help_text=[character(len=80) :: &
      "usage: topogtools [--help] <command> [<args>]", &
      "", &
      "Collection of tools to edit and manipulate ocean model topographies.", &
      "See 'topogtools --help <command>' to read about a specific subcommand.", &
      "", &
      "Available coomands:", &
      "  gen_topo - Generate a new topography file from a bathymetry", &
      "  deseas - Remove enclosed seas", &
      "  min_max_depth - Set minimum and maximum depth", &
      "  check_nonadvective - Check for cells that are nonadvective", &
      "  fix_nonadvective - Fix cells that are non-advective", &
      "  topog2mask - Generate mask", &
      "" ]

    ! Print help in case the user specified the --help flag
    call set_args(' ', help_text)

    ! Also print even if the user did not specify --help
    if (.not. lget('help')) then
      write(output_unit,'(g0)') (trim(help_text(ii)), ii=1, size(help_text))
      stop
    end if

  case default
    write(*,'(3a)') "topogtools: '", trim(name), "' is not a topogtools subcommand. See 'topogtools --help'."
    stop
  end select

  ! Sanity checks for common arguments
  file_in = sget('input')
  if (file_in == 'unset') then
    write(*,*) 'ERROR: no input file specified'
    stop
  end if
  call check_file_exist(file_in)

  select case (name)
  case ('gen_topo', 'deseas', 'min_max_depth', 'fix_nonadvective', 'topog2mask')
    file_out = sget('output')
    if (file_out == 'unset') then
      write(*,*) 'ERROR: no output file specified'
      stop
    end if
  end select

  select case (name)
  case ('min_max_depth', 'fix_nonadvective', 'check_nonadvective')
    vgrid = sget('vgrid')
    call check_file_exist(vgrid)
  end select

  ! Run subcommand
  select case (name)
  case ('gen_topo')
    hgrid = sget('hgrid')
    if (hgrid == 'unset') then
      write(*,*) 'ERROR: no hgrid file specified'
      stop
    end if
    call check_file_exist(hgrid)
    call gen_topo(file_in, file_out, hgrid, lget('tripolar'), rget('longitude-offset'))

  case ('deseas')
    call deseas(file_in, file_out)

  case ('min_max_depth')
    call min_max_depth(file_in, file_out, vgrid, iget('level'))

  case ('fix_nonadvective')
    call fix_nonadvective(file_in, file_out, vgrid)

  case ('check_nonadvective')
    call check_nonadvective(file_in, vgrid)

  case ('topog2mask')
    call topog2mask(file_in, file_out, rget('fraction'))

  end select

end program topogtools