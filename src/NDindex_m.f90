MODULE NDindex_m
  USE NumParameters_m
  USE UtilLib_m
  IMPLICIT NONE

  TYPE :: NDindex_t
    integer                       :: Ndim            = 0
    integer                       :: Tab0_i(2)       = [0,1]
    integer, allocatable          :: NDsize(:)
    integer, allocatable          :: NDend(:)
  END TYPE NDindex_t

  PRIVATE
  PUBLIC :: increase_NDindex,Init_tab_ind,Testindex,NDindex_t!,alloc_Tab

CONTAINS

  !SUBROUTINE alloc_Tab(Tab_i,NDindex)
  !  TYPE(NDindex_t),  intent(in) :: NDindex
  !  integer,     intent(inout)   :: Tab_i(:)
   !logical,    parameter      :: debug = .true.
  !  logical,     parameter      :: debug = .false.

  !  IF (debug) THEN
  !    write(out_unitp,*) 'BEGINNING alloc_Tab'
  !    flush(out_unitp)
  !  END IF

  !  IF ( NDindex%Ndim < 1) STOP 'ERROR in init_Tab: Ndim < 1!'
  !  IF (allocated(Tab_i)) deallocate(Tab_i)

  !  allocate(Tab_i(NDindex%Ndim))

  !  IF (debug) THEN
  !    write(out_unitp,*) 'END alloc_Tab'
  !    flush(out_unitp)
  !  END IF

  !END SUBROUTINE alloc_Tab

  SUBROUTINE Init_tab_ind(Tab_i,NDindex)
    USE UtilLib_m
    IMPLICIT NONE
    TYPE(NDindex_t),  intent(in) :: NDindex
    integer,     intent(inout)   :: Tab_i(:)

    !logical,    parameter       :: debug = .true.
    logical,     parameter       ::debug = .false.

    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING Init_tab_ind'
      flush(out_unitp)
    END IF
     Tab_i(:)= NDindex%Tab0_i(:)

    !Tab_i(1)=0
    !Tab_i(2)=1

    IF (debug) THEN
      write(out_unitp,*) 'END Init_tab_ind'
      flush(out_unitp)
    END IF

  END SUBROUTINE Init_tab_ind

  SUBROUTINE increase_NDindex(Tab_i,NDindex,n1)
    USE UtilLib_m
    IMPLICIT NONE

    TYPE(NDindex_t),  intent(in) :: NDindex
    integer,     intent(inout)   :: Tab_i(:)
    integer,     intent(in)      :: n1
    !logical,    parameter       :: debug = .true.
    logical,     parameter       ::debug = .false.

    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING Tab_ind'
      flush(out_unitp)
    END IF

    IF(Tab_i(1) == n1) THEN

      Tab_i(2)=Tab_i(2)+1
      Tab_i(1)=1
    ELSE
      Tab_i(1)=Tab_i(1)+1
    END IF

    IF (debug) THEN
      write(out_unitp,*) 'END Tab_ind'
      flush(out_unitp)
    END IF

  END SUBROUTINE increase_NDindex


  SUBROUTINE Testindex(NDindex)
    USE UtilLib_m
    IMPLICIT NONE

    TYPE(NDindex_t),  intent(in) :: NDindex
    integer                     :: Tab_i(2)
  !  integer,     intent(in)      :: nl
    integer                      :: i
    !logical,    parameter      :: debug = .true.
    logical,     parameter      ::debug = .false.

    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING Testindex'
      flush(out_unitp)
    END IF

    !Call alloc_Tab(Tab_i,NDindex)
    Call Init_tab_ind(Tab_i,NDindex)
    DO i= 1,100
      CALL increase_NDindex(Tab_i,NDindex,6)
      write(out_unitp,*) i, tab_i(:)
    END DO

    IF (debug) THEN
      write(out_unitp,*) 'END Testindex'
      flush(out_unitp)
    END IF

  END SUBROUTINE Testindex
END MODULE NDindex_m
