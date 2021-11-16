!###############################################################################
!#                                                                             #
!# aed_csv_reader.F90                                                          #
!#                                                                             #
!# Read csv input files for libaed2 variable spatial initialisation.           #
!#                                                                             #
!#   -----------------------------------------------------------------------   #
!#                                                                             #
!#  Developed by :                                                             #
!#      AquaticEcoDynamics (AED) Group                                         #
!#      School of Agriculture and Environment                                  #
!#      The University of Western Australia                                    #
!#                                                                             #
!#      http://aquatic.science.uwa.edu.au/                                     #
!#                                                                             #
!#  Copyright 2013 - 2021 -  The University of Western Australia               #
!#                                                                             #
!#   GLM is free software: you can redistribute it and/or modify               #
!#   it under the terms of the GNU General Public License as published by      #
!#   the Free Software Foundation, either version 3 of the License, or         #
!#   (at your option) any later version.                                       #
!#                                                                             #
!#   GLM is distributed in the hope that it will be useful,                    #
!#   but WITHOUT ANY WARRANTY; without even the implied warranty of            #
!#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             #
!#   GNU General Public License for more details.                              #
!#                                                                             #
!#   You should have received a copy of the GNU General Public License         #
!#   along with this program.  If not, see <http://www.gnu.org/licenses/>.     #
!#                                                                             #
!#   -----------------------------------------------------------------------   #
!#                                                                             #
!# Created Apr 2016                                                            #
!#                                                                             #
!###############################################################################

MODULE aed_csv_reader

   IMPLICIT NONE

!-------------------------------------------------------------------------------
!
!  PRIVATE     ! by default, make everything private
!
#ifdef SINGLE
#  define AED_REAL real(4)
#else
#  define AED_REAL real(8)
#endif

   INTEGER,PARAMETER :: bufsize=2048


   !----------------------------------------------------------------------------
   !# The AED_READER user type is used internally in the parser to track the   #
   !#  file being read.                                                        #
   !----------------------------------------------------------------------------
   TYPE AED_READER
     CHARACTER(len=bufsize) :: buf
     INTEGER                :: buf_pos
     INTEGER                :: buf_len
     INTEGER                :: lun
     INTEGER                :: n_cols
   END TYPE AED_READER

   !----------------------------------------------------------------------------
   !# The AED_SYMBOL type is a form of variable length string                  #
   !----------------------------------------------------------------------------
   TYPE AED_SYMBOL
     INTEGER           :: length
     CHARACTER,POINTER :: sym(:)
   END TYPE AED_SYMBOL

!-------------------------------------------------------------------------------

    CHARACTER(len=1), parameter :: EOLN=achar(10)
    CHARACTER(len=1), parameter :: CTAB=achar(8)
    CHARACTER(len=1), parameter :: CNUL=achar(0)
!
    CHARACTER(len=32) :: t_strs(0:32)
!
   !TYPE(AED_READER),POINTER :: units(10)
    TYPE ARP
       TYPE(AED_READER),POINTER :: p
    END TYPE ARP
    TYPE(ARP) :: units(10)
!
!-------------------------------------------------------------------------------
   PUBLIC aed_csv_read_header, aed_csv_read_row, aed_csv_close, AED_SYMBOL
   PUBLIC extract_double, extract_logical, extract_integer, extract_string
   PUBLIC copy_name

CONTAINS

!###############################################################################
SUBROUTINE init_t_strs
!-------------------------------------------------------------------------------
    t_strs( 0) = ""
    t_strs( 1) = " "
    t_strs( 2) = "  "
    t_strs( 3) = "   "
    t_strs( 4) = "    "
    t_strs( 5) = "     "
    t_strs( 6) = "      "
    t_strs( 7) = "       "
    t_strs( 8) = "        "
    t_strs( 9) = "         "
    t_strs(10) = "          "
    t_strs(11) = "           "
    t_strs(12) = "            "
    t_strs(13) = "             "
    t_strs(14) = "              "
    t_strs(15) = "               "
    t_strs(16) = "                "
    t_strs(17) = "                 "
    t_strs(18) = "                  "
    t_strs(19) = "                   "
    t_strs(20) = "                    "
    t_strs(21) = "                     "
    t_strs(22) = "                      "
    t_strs(23) = "                       "
    t_strs(24) = "                        "
    t_strs(25) = "                         "
    t_strs(26) = "                          "
    t_strs(27) = "                           "
    t_strs(28) = "                            "
    t_strs(29) = "                             "
    t_strs(30) = "                              "
    t_strs(31) = "                               "
    t_strs(32) = "                                "
END SUBROUTINE init_t_strs
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
INTEGER FUNCTION f_get_lun()
!-------------------------------------------------------------------------------
! Find the first free logical unit number
!-------------------------------------------------------------------------------
   INTEGER :: lun
   LOGICAL :: opened
!
!-------------------------------------------------------------------------------
!BEGIN
   DO lun = 10,99
      inquire(unit=lun, opened=opened)
      IF ( .not. opened ) THEN
         f_get_lun = lun
         RETURN
      ENDIF
   ENDDO
   f_get_lun = -1
END FUNCTION f_get_lun
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
INTEGER FUNCTION fopen(filename)
!-------------------------------------------------------------------------------
! open a file and return it logical unit number
!-------------------------------------------------------------------------------
!ARGUMENTS
   CHARACTER(len=*),INTENT(in) :: filename
!
!-------------------------------------------------------------------------------
!LOCALS
   INTEGER :: iostat=-1
   INTEGER :: lun=-1

!BEGIN
   lun = f_get_lun()
   IF ( lun .GT. 0 ) open(lun, status='old', file=filename, iostat=iostat)

   IF (iostat .NE. 0) lun = -1

   fopen = lun
END FUNCTION fopen
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
LOGICAL FUNCTION init_parser(filename, aedr)
!-------------------------------------------------------------------------------
! open a file and initialise the parser structure for it
!-------------------------------------------------------------------------------
!ARGUMENTS
   CHARACTER(len=*), INTENT(in) :: filename
   TYPE(AED_READER),POINTER     :: aedr

!LOCALS
   INTEGER :: lun
!
!-------------------------------------------------------------------------------
!BEGIN
   CALL init_t_strs
   lun = fopen(filename)
   IF ( lun .gt. 0 ) THEN
      ALLOCATE(aedr)
      aedr%lun=lun
      aedr%buf_pos=-1
      aedr%buf_len=0
      init_parser = .true.
   ELSE
      init_parser = .false.
   ENDIF
END FUNCTION init_parser
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
INTEGER FUNCTION char_in_str(ch, str, s)
!-------------------------------------------------------------------------------
! char_in_str returns the index of the first match of character ch in string str
!-------------------------------------------------------------------------------
!ARGUMENTS
   CHARACTER, INTENT(in) :: ch
   CHARACTER(len=*), INTENT(in) :: str
   INTEGER, optional, INTENT(in) :: s

!LOCALS
   INTEGER res
   INTEGER lnt
!
!-------------------------------------------------------------------------------
!BEGIN
   lnt=LEN_TRIM(str)
   IF (present(s)) THEN
      res=s
   ELSE
      res=1
   ENDIF
   DO WHILE ( (str(res:res) .NE. ch) .AND. (res .LE. lnt) )
      res=res+1
   ENDDO
   IF (res .GT. lnt) res = 0
   char_in_str = res
END FUNCTION char_in_str
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
LOGICAL FUNCTION next_symbol(aedr, sym)
!-------------------------------------------------------------------------------
! get the next "symbol" from the file, if there is one
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE(AED_READER),POINTER      :: aedr
   TYPE(AED_SYMBOL),INTENT(out)  :: sym

!LOCALS
   INTEGER :: iostat, s1, e1, i, j
   CHARACTER :: quot
   CHARACTER(len=4) :: term
   LOGICAL   :: endln
!
!-------------------------------------------------------------------------------
!BEGIN
   IF ( ASSOCIATED(sym%sym) ) DEALLOCATE(sym%sym)
   next_symbol = .false.
   endln = .false.
   quot=CNUL

   IF (aedr%buf_pos .GE. 0) THEN
      IF (aedr%buf_pos .GT. aedr%buf_len) THEN
         aedr%buf_pos = -1
         ALLOCATE(sym%sym(1))
         sym%length=1
         sym%sym(1) = EOLN
         next_symbol = .true.
         RETURN
      ENDIF
   ENDIF

   term='   '
   term(1:1)='"'
   term(2:2)="'"
   term(3:3)=","

   s1=bufsize+1
   DO WHILE ( s1 .GT. bufsize )
      DO WHILE ( (aedr%buf_pos .LE. 0) .OR. (aedr%buf_pos .GT. aedr%buf_len) )
         read(UNIT=aedr%lun,FMT='(A)', iostat=iostat) aedr%buf
         IF ( iostat .NE. 0) RETURN

         aedr%buf_len=LEN_TRIM(aedr%buf)
         IF ( (aedr%buf_len .GT. 0) .AND.         &
              (aedr%buf(1:1) .NE. '#') .AND.      &
              (aedr%buf(1:1) .NE. '!') ) THEN
            aedr%buf_pos=1
         ELSE
            aedr%buf_pos=-1
         ENDIF
      ENDDO

      s1=aedr%buf_pos

      !* Skip leading blanks
      DO WHILE( ( (aedr%buf(s1:s1) .EQ. ' ') .OR.         &
                  (aedr%buf(s1:s1) .EQ. CTAB) .OR.        &
                  (aedr%buf(s1:s1) .EQ. CNUL) .OR.        &
                  (aedr%buf(s1:s1) .EQ. EOLN) ) .AND.     &
                  (s1 .LE. aedr%buf_len) )
         s1=s1+1
      ENDDO

      !* If we have a # skip the rest of the line
      IF (aedr%buf(s1:s1) .EQ. '#') THEN
         s1=bufsize+1
         aedr%buf_pos=aedr%buf_len+1
      ENDIF
   ENDDO

   IF (aedr%buf(s1:s1) .EQ. '"') quot='"'
   IF (aedr%buf(s1:s1) .EQ. "'") quot="'"

   IF (quot .NE. CNUL) THEN
      !* looking for an end of quote
      s1=s1+1
      e1=s1
      DO WHILE((e1 .LE. aedr%buf_len) .AND. (aedr%buf(e1:e1) .NE. quot))
         e1=e1+1
      ENDDO
   ELSE
      e1=s1+1
      DO WHILE((e1 .LE. aedr%buf_len) .AND. (char_in_str(aedr%buf(e1:e1),term) .EQ. 0))
         e1=e1+1
      ENDDO
   ENDIF

   next_symbol = .true.
   ALLOCATE(sym%sym(e1-s1))
   sym%length=e1-s1
   j=1
   DO i=s1,e1-1
      sym%sym(j:j)=aedr%buf(i:i)
      j=j+1
   ENDDO

   IF (quot .NE. CNUL) e1=e1+1

   if (aedr%buf(e1:e1) == ",") e1 = e1+1
   aedr%buf_pos=e1
END FUNCTION next_symbol
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
LOGICAL FUNCTION end_parse(aedr)
!-------------------------------------------------------------------------------
! Close the file and free the parser storage
!-------------------------------------------------------------------------------
!ARGUMENTS
    TYPE(AED_READER),POINTER :: aedr

!LOCALS
    INTEGER :: iostat
!
!-------------------------------------------------------------------------------
!BEGIN
    close(aedr%lun, iostat=iostat)
    DEALLOCATE(aedr)
    end_parse=(iostat .eq. 0)
END FUNCTION end_parse
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION extract_double(sym) RESULT(num)
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!ARGUMENTS
    TYPE(AED_SYMBOL)  :: sym

!LOCALS
    DOUBLE PRECISION  :: num
    CHARACTER(len=80) :: tbuf
    INTEGER           :: i
!
!-------------------------------------------------------------------------------
!BEGIN
    DO i=1,sym%length
       tbuf(i:i)=sym%sym(i)
    ENDDO
    tbuf(sym%length+1:)=' '

    read(tbuf,*) num
END FUNCTION extract_double
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION extract_integer(sym) RESULT(num)
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!ARGUMENTS
    TYPE(AED_SYMBOL)  :: sym

!LOCALS
    INTEGER           :: num
    CHARACTER(len=80) :: tbuf
    INTEGER           :: i
!
!-------------------------------------------------------------------------------
!BEGIN
    DO i=1,sym%length
       tbuf(i:i)=sym%sym(i)
    ENDDO
    tbuf(sym%length+1:)=' '

    read(tbuf,*) num
END FUNCTION extract_integer
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION extract_logical(sym) RESULT(res)
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!ARGUMENTS
    TYPE(AED_SYMBOL)  :: sym

!LOCALS
    LOGICAL           :: res
    CHARACTER(len=80) :: tbuf
    INTEGER           :: i
!
!-------------------------------------------------------------------------------
!BEGIN
    DO i=1,sym%length
       tbuf(i:i)=sym%sym(i)
    ENDDO
    tbuf(sym%length+1:)=' '

    read(tbuf,*) res
END FUNCTION extract_logical
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION extract_string(sym) RESULT(str)
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!ARGUMENTS
    TYPE(AED_SYMBOL)  :: sym

!LOCALS
    CHARACTER(len=20) :: str
    INTEGER           :: i
!
!-------------------------------------------------------------------------------
!BEGIN
    DO i=1,sym%length
       str(i:i)=sym%sym(i)
    ENDDO
    str(sym%length+1:)=' '
END FUNCTION extract_string
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE copy_name(sym, name)
!-------------------------------------------------------------------------------
!ARGUMENTS
!-------------------------------------------------------------------------------
   TYPE(AED_SYMBOL),INTENT(in) :: sym
   CHARACTER(len=*),INTENT(out) :: name
!
!LOCALS
   INTEGER i
!
!-------------------------------------------------------------------------------
!BEGIN
   name = t_strs(sym%length)
   DO i=1,sym%length
      name(i:i) = sym%sym(i)
   ENDDO
END SUBROUTINE copy_name
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION scan_csv_header(aedr,titles) RESULT(count)
!-------------------------------------------------------------------------------
!ARGUMENTS
    TYPE(AED_READER),POINTER,INTENT(in) :: aedr
    CHARACTER(len=32),POINTER,INTENT(out) :: titles(:)
!
!LOCALS
    TYPE(AED_SYMBOL) :: sym
    INTEGER          :: count, i
!
!-------------------------------------------------------------------------------
!BEGIN
   NULLIFY(sym%sym)
   count = 0

   DO WHILE( next_symbol(aedr, sym) )
      IF ( sym%sym(1) .EQ. EOLN ) THEN
         REWIND(aedr%lun)
         aedr%buf_pos=-1
         aedr%buf_len=0
         i = 0;
         ALLOCATE(titles(count))
         titles = ' '
         DO WHILE( next_symbol(aedr, sym) )
            IF ( sym%sym(1) .EQ. EOLN ) RETURN
            i = i + 1
            CALL copy_name(sym,titles(i))
         ENDDO
         RETURN
      ENDIF
      count=count+1
   ENDDO
END FUNCTION scan_csv_header
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
INTEGER FUNCTION aed_csv_read_header(fname, names, ncols)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CHARACTER(*),INTENT(in) :: fname
   CHARACTER(len=32),DIMENSION(:),POINTER,INTENT(inout) :: names
   INTEGER,INTENT(out) :: ncols
!
!LOCALS
   INTEGER :: unit, i
   TYPE(AED_READER),POINTER :: aedr
   LOGICAL                  :: res
!
!-------------------------------------------------------------------------------
!BEGIN
   NULLIFY(aedr)
   unit = 0

   IF ( .NOT. init_parser(fname, aedr) ) THEN
      print*, "Failed to open file '",fname,"'"
      aed_csv_read_header = -1
      RETURN
   ENDIF
   DO i=1,10
      IF ( .NOT. ASSOCIATED(units(i)%p) ) THEN
         units(i)%p => aedr
         unit = i
         exit
      ENDIF
   ENDDO
   IF ( unit == 0 ) THEN
      res = end_parse(aedr)
   ELSE
      ncols = scan_csv_header(aedr, names)
   ENDIF

   aedr%n_cols = ncols
   aed_csv_read_header = unit
END FUNCTION aed_csv_read_header
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
LOGICAL FUNCTION aed_csv_read_row(unit, values)
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in) :: unit
   TYPE(AED_SYMBOL),DIMENSION(:),INTENT(out) :: values
!
!LOCALS
   TYPE(AED_READER),POINTER :: aedr
   INTEGER :: i, j, ncols
   TYPE(AED_SYMBOL) :: sym
!
!-------------------------------------------------------------------------------
!BEGIN
   aedr => units(unit)%p
   ncols = aedr%n_cols

   values(1:ncols)%length = 0
   i = 0
   DO WHILE ( next_symbol(aedr, sym) ) !#
      IF ( sym%sym(1) .EQ. EOLN ) EXIT
      i = i + 1
      IF ( i <= ncols ) THEN
         IF ( ASSOCIATED(values(i)%sym) ) NULLIFY( values(i)%sym )

         values(i) = sym
         NULLIFY(sym%sym)
      ENDIF
   ENDDO
   IF ( i > 0 .AND. i /= ncols ) &
      print *, "data row had ", i, " columns : expecting ", ncols

   aed_csv_read_row = (i > 0)
END FUNCTION aed_csv_read_row
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
LOGICAL FUNCTION aed_csv_close(unit)
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in) :: unit
!
!LOCALS
   TYPE(AED_READER),POINTER :: aedr
!
!-------------------------------------------------------------------------------
!BEGIN
   aed_csv_close = .FALSE.
   aedr => units(unit)%p
   IF (ASSOCIATED(aedr)) aed_csv_close = end_parse(aedr)
   NULLIFY(aedr)
   units(unit)%p => aedr
END FUNCTION aed_csv_close
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE aed_csv_reader
