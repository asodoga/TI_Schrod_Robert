MODULE NDindex_m
  USE NumParameters_m
!  USE Basis_m, only : Basis_t
  USE UtilLib_m
  IMPLICIT NONE
  TYPE :: Tab_t

    integer, allocatable           :: Tab_i(:)

  END TYPE Tab_t

  PRIVATE
  PUBLIC :: tab_ind,Int_tab_ind,alloc_Tab

 contains

 SUBROUTINE alloc_Tab(Tab,nb)
    TYPE(Tab_t),  intent(inout) :: Tab
    integer,     intent(in)    :: nb
    !logical,    parameter      :: debug = .true.
    logical,     parameter      ::debug = .false.

    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING alloc_Tab'
      !call write_basis(Basis)
      flush(out_unitp)
    END IF

    IF (nb < 1) STOP 'ERROR in init_Tab: nb < 1!'

    allocate(Tab%Tab_i(nb))
    IF (debug) THEN
      write(out_unitp,*) 'END alloc_Tab'
      flush(out_unitp)
    END IF
  END SUBROUTINE alloc_Tab

  SUBROUTINE Int_tab_ind(Tab)
    USE UtilLib_m
    IMPLICIT NONE

    TYPE(Tab_t),  intent(inout) :: Tab
    !logical,    parameter      :: debug = .true.
    logical,     parameter      ::debug = .false.

    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING Int_tab_ind'
      !call write_basis(Basis)
      flush(out_unitp)
    END IF

    Tab%Tab_i(1)=0
    Tab%Tab_i(2)=1

    IF (debug) THEN
      write(out_unitp,*) 'END Int_tab_ind'
      flush(out_unitp)
    END IF

  END SUBROUTINE Int_tab_ind

  SUBROUTINE Tab_ind(Tab,n1)
    USE UtilLib_m
    IMPLICIT NONE

      TYPE(Tab_t),  intent(inout) :: Tab
      integer,     intent(in)     :: n1
      !logical,    parameter      :: debug = .true.
      logical,     parameter      ::debug = .false.

      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING Tab_ind'
        !call write_basis(Basis)
        flush(out_unitp)
      END IF

      IF(Tab%Tab_i(1) == n1) THEN
        Tab%Tab_i(2)=Tab%Tab_i(2)+1
        Tab%Tab_i(1)=1
      ELSE
        Tab%Tab_i(1)=Tab%Tab_i(1)+1
      END IF

      IF (debug) THEN
        write(out_unitp,*) 'END Tab_ind'
        flush(out_unitp)
      END IF

   END SUBROUTINE Tab_ind


END MODULE NDindex_m
