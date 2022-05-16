MODULE NDindex_m
  USE NumParameters_m
  USE UtilLib_m
  IMPLICIT NONE

  TYPE :: NDindex_t
    integer                       :: Ndim            = 0
    integer,  allocatable         :: Tab0(:)      ! = [0,1]
    integer, allocatable          :: NDsize(:)
    integer, allocatable          :: NDend(:)
  END TYPE NDindex_t

  PRIVATE
  PUBLIC :: increase_NDindex,Init_tab_ind,Testindex,NDindex_t,Init_NDindex!,alloc_Tab

CONTAINS


  SUBROUTINE Init_NDindex(NDindex,NDend,Ndim)
    TYPE(NDindex_t),intent(inout) :: NDindex
    integer ,      intent(in)     :: NDend(:)
    integer ,      intent(in)     :: Ndim
    !logical,    parameter         :: debug = .true.
    logical,     parameter        :: debug = .false.

    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING Init_NDindex'
      flush(out_unitp)
    END IF

    NDindex%Ndim = Ndim
    NDindex%NDend = NDend
    Allocate(NDindex%Tab0(NDindex%Ndim))
    NDindex%Tab0(:)   = 1
    NDindex%Tab0(1)   = 0

    IF (debug) THEN
      write(out_unitp,*) NDindex%Ndim
      write(out_unitp,*)'NDindex%Ndend', NDindex%Ndend
      write(out_unitp,*)'NDindex%Tab0', NDindex%Tab0(:)
      write(out_unitp,*) 'END Init_NDindex'

      flush(out_unitp)
    END IF

  END SUBROUTINE Init_NDindex

  SUBROUTINE Init_tab_ind(Tab_ind,NDindex)
    USE UtilLib_m
    IMPLICIT NONE
    TYPE(NDindex_t),  intent(in) :: NDindex
    integer,     intent(inout)   :: Tab_ind(:)

    !logical,    parameter       :: debug = .true.
    logical,     parameter       :: debug = .false.

    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING Init_tab_ind'
      flush(out_unitp)
    END IF

     Tab_ind(:)= NDindex%Tab0(:)

    IF (debug) THEN
      write(out_unitp,*) 'END Init_tab_ind'
      flush(out_unitp)
    END IF

  END SUBROUTINE Init_tab_ind

  SUBROUTINE increase_NDindex(Tab_ind,NDindex,Endloop)
    USE UtilLib_m
    IMPLICIT NONE

    TYPE(NDindex_t),  intent(in) :: NDindex

    integer,     intent(inout)   :: Tab_ind(:)
    !logical,    parameter       :: debug = .true.
    logical,     parameter       :: debug = .false.
    logical,     intent(inout)   :: Endloop
    integer                      :: i


    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING Tab_ind'
      write(out_unitp,*)'NDindex%Ndend', NDindex%Ndend
      write(out_unitp,*)'Tab_ind1', Tab_ind
      flush(out_unitp)
    END IF

    Tab_ind(1)=Tab_ind(1)+1
    IF (debug)  write(out_unitp,*)'Tab_indd', Tab_ind
    DO i=1,NDindex%Ndim-1

       IF(Tab_ind(i) > NDindex%Ndend(i)) THEN
         Tab_ind(i+1)=Tab_ind(i+1)+1
         Tab_ind(i)=1
       END IF
    END DO
    IF (debug) write(out_unitp,*)'Tab_indfin', Tab_ind(NDindex%Ndim)

    Endloop =(Tab_ind(NDindex%Ndim)==NDindex%Ndend(NDindex%Ndim)+1)


    !IF(Tab_ind(1) == NDindex%Ndend(1)) THEN

     ! IF(Tab_ind(2) == NDindex%Ndend(2)) THEN
     !  Tab_ind(3)=Tab_ind(3)+1
     !  Tab_ind(2)=1
     !  Tab_ind(1)=1
     ! ELSE
     !  Tab_ind(2)=Tab_ind(2)+1
     !  Tab_ind(1)=1

     ! END IF

   ! ELSE
   !   Tab_ind(1)=Tab_ind(1)+1
    !END IF

    !Endloop =(Tab_ind(3)>NDindex%Ndend(3))
    !  IF(Tab_ind(1) == NDindex%Ndend(1)) THEN
    !    Tab_ind(2)=Tab_ind(2)+1
    !    Tab_ind(1)=1
    !  ELSE
    !    Tab_ind(1)=Tab_ind(1)+1
    !  END IF
    !END DO

    !Endloop =(Tab_ind(NDindex%Ndim)>NDindex%Ndend(NDindex%Ndim))

    IF (debug) THEN
      write(out_unitp,*) 'END Tab_ind'
      flush(out_unitp)
    END IF

  END SUBROUTINE increase_NDindex

    SUBROUTINE Write_NDindex(NDindex)
      USE UtilLib_m
      IMPLICIT NONE

      TYPE(NDindex_t),  intent(in) :: NDindex
      !logical,    parameter      :: debug = .true.
      logical,     parameter      :: debug = .false.

      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING Write_NDindex'
        flush(out_unitp)
      END IF

      write(out_unitp,*)'NDindex%Ndend=', NDindex%Ndend

      write(out_unitp,*) 'NDindex%Ndim =',NDindex%Ndim
      write(out_unitp,*)'NDindex%Ndend =', NDindex%Ndend
      write(out_unitp,*)'NDindex%Tab0 =', NDindex%Tab0(:)

      IF (debug) THEN
        write(out_unitp,*) ' END Write_NDindex'
        flush(out_unitp)
      END IF

    END SUBROUTINE Write_NDindex

  SUBROUTINE Testindex(NDindex)
    USE UtilLib_m
    IMPLICIT NONE

    TYPE(NDindex_t),  intent(in) :: NDindex
    integer, allocatable         :: Tab_ind(:)
    logical                      :: Endloop
  !  integer,     intent(in)      :: nl
    integer                      :: i
    !logical,    parameter      :: debug = .true.
    logical,     parameter      :: debug = .false.

    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING Testindex'
      flush(out_unitp)
    END IF
    Allocate(Tab_ind(NDindex%Ndim))
    Call Init_tab_ind(Tab_ind,NDindex)

    i=0
    !write(out_unitp,*) i, tab_ind(:)
    DO !i= 1,100
      i=i+1
      CALL increase_NDindex(Tab_ind,NDindex,Endloop)
        IF (Endloop) exit
      write(out_unitp,*) i, tab_ind(:)
    END DO

    IF (debug) THEN
      write(out_unitp,*) 'END Testindex'
      flush(out_unitp)
    END IF

  END SUBROUTINE Testindex
END MODULE NDindex_m
