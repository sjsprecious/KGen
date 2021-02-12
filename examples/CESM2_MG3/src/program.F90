PROGRAM main
    USE micro_mg_cam, only : micro_mg_cam_tend_pack
    USE shr_kind_mod, only : r8 => SHR_KIND_R8
    USE netcdf

    REAL(r8)           :: dtime = 1800._r8
    INTEGER            :: ncid, varid, nlev, mgncol
    CHARACTER(LEN=128) :: FILE_NAME
    CHARACTER(LEN=5)   :: col, lev 

    ! see examples at:
    !     https://stackoverflow.com/questions/1262695/convert-integers-to-strings-to-create-output-filenames-at-run-time
    write(col, '(I5)') PCOLS
    write(lev, '(I5)') PVER
    FILE_NAME = "./input_data/PCOLS"//trim(adjustl(col))//"_PVER"//trim(adjustl(lev))//"/mg3_kgen_input.nc"

    ! open an existing netcdf file
    CALL check( nf90_open(FILE_NAME, NF90_WRITE, ncid) )

    ! see examples at:
    !     https://www.nag.com/market/training/fortran-workshop/netcdf-f90.pdf
    CALL check( nf90_inq_dimid(ncid, "mgncol", varid) )
    CALL check( nf90_inquire_dimension(ncid, varid, len=mgncol) )

    CALL check( nf90_inq_dimid(ncid, "nlev", varid) )
    CALL check( nf90_inquire_dimension(ncid, varid, len=nlev) )

    CALL micro_mg_cam_tend_pack ( dtime, mgncol, nlev, ncid )

    print *, "Program done..."

contains

subroutine check(status)

   integer, intent(in) :: status

   if (status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      stop "Stopped"
   end if

end subroutine check

END PROGRAM
