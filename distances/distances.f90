program main

!*****************************************************************************80
!
!! MAIN is the main program for DISTANCES.
!
!  Discussion:
!
!    DISTANCES reads a multiple alignment file and outputs a distance matrix.
!
!    This program reads a multiple alignment file in FASTA or PIR format 
!    and outputs a distance matrix.  The distance matrix is intended as
!    input to the SPACER routine, which performs principal coordinate
!    analysis.  The sequences must be aligned first.  The gaps in the 
!    alignment are read as hyphen characters: "-".
!
!    The READ_ALIGN program may be used to convert a multiple sequence
!    data file from a related format to PIR format.
!
!    The maximum problem size is set for 5000 residues each in 650 sequences.
!    Reset these in ALL subroutines where they occur if you change them!
!
!  Author:
!
!    Des Higgins, 
!    EMBL Data Library, 
!    April 1992
!
!  Modified:
!
!    05 April 2000
!
  implicit none

  integer ( kind = 4 ), parameter :: sequence_length_max = 5000
  integer ( kind = 4 ), parameter :: sequence_max = 650

  character ( len = 20 ) aacids
  real aamat(0:20,0:20)
  integer ( kind = 4 ) conlen
  real distance(sequence_max,sequence_max)
  logical dna
  integer ( kind = 4 ) i
  character ( len = 100 ) input_name
  integer ( kind = 4 ), parameter :: input_unit = 66
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) it
  integer ( kind = 4 ) j
  character ( len = 100 ) output_name
  integer ( kind = 4 ), parameter :: output_unit = 67
  character sequence(sequence_max,sequence_length_max)
  integer ( kind = 4 ) sequence_index(sequence_max,sequence_length_max)
  integer ( kind = 4 ) sequence_length(sequence_max)
  character ( len = 16 ) sequence_name(sequence_max)
  integer ( kind = 4 ) sequence_number

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DISTANCES'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Read a multiple alignment file;'
  write ( *, '(a)' ) '  Write a distance matrix for input to SPACER.'
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Enter the name of the input alignment file:'
  read ( *, '(a)' ) input_name

  open ( unit = input_unit, file = input_name, status = 'old', iostat = ios ) 

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DISTANCES - Fatal error!'
    write ( *, '(a)' ) '  Could not open the alignment file.'
    stop
  end if

  call read_sequence ( input_unit, sequence_max, sequence_length_max, &
    sequence_number, sequence, sequence_name, sequence_length, conlen )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The sequence alignment file was read.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of sequences = ', sequence_number
  write ( *, '(a,i6)' ) '  Consensus length    = ', conlen
  write ( *, '(a)' ) ' '

  do i = 1, sequence_number
    write ( *, '(a,i4,2x,a,a,i5)' ) ' Sequence # ', i, sequence_name(i), &
      ' length = ', sequence_length(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' Type 1 to toss all gaps; '
  write ( *, '(a)' ) '      0 not to.'
  read ( *, * ) it

  if ( it == 1 ) then
    call toss_gaps ( sequence_max, sequence_length_max, sequence_number, &
      sequence, conlen )
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' Type 1 for DNA; '
  write ( *, '(a)' ) '      0 for protein:'
  read ( *, * ) it

  if ( it == 1 ) then
    dna = .true.
  end if

  if ( dna ) then
    aacids = 'ACGTU'
  else
    aacids = 'DEKRHNQSTILVFWYCMAGP'
  end if

  aamat(0:20,0:20) = 1.0E+00

  do i = 0, 20
    aamat(i,i) = 0.0E+00
  end do

  aamat(0,0) = 1.0E+00
  aamat(4,5) = 0.0E+00
  aamat(5,4) = 0.0E+00

  if ( .not. dna ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Type 1 for ID distances;'
    write ( *, '(a)' ) '     0 for Smith AA matrix:'
    read ( *, * ) it

    if ( it == 0 ) then
      call set_aamat ( aamat, aacids )
    end if

  end if
!
!  Replace the sequence character by the equivalent sequence index number.
!
  do i = 1, sequence_number
    do j = 1, conlen
      sequence_index(i,j) = index ( aacids, sequence(i,j) )
    end do
  end do

  call dist_mat_set ( sequence_max, sequence_length_max, sequence_number, &
    conlen, distance, sequence_index, aamat )
!
!  Write the distance matrix to a file.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Enter name for the output distance matrix file:'
  read ( *, '(a)' ) output_name

  open ( unit = output_unit, file = output_name, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DISTANCES - Fatal error!'
    write ( *, '(a)' ) '  Could not open the distance matrix file.'
    stop
  end if

  call dist_mat_write ( output_unit, sequence_max, sequence_number, distance )

  close ( unit = output_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DISTANCES'
  write ( *, '(a)' ) '  Normal end of execution.'
 
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine dist_mat_set ( sequence_max, sequence_length_max, sequence_number, &
  conlen, distance, sequence_index, aamat )

!*****************************************************************************80
!
!! DIST_MAT_SET calculates a matrix of distances between pairs of sequences.
!
!  Author:
!
!    Des Higgins, 
!    EMBL Data Library, 
!    April 1992
!
!  Modified:
!
!    06 April 2000
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) SEQUENCE_MAX, the maximum number of sequences.
!
!    Input, integer ( kind = 4 ) SEQUENCE_LENGTH_MAX, the maximum length of a sequence.
!
!    Input, integer ( kind = 4 ) SEQUENCE_NUMBER, the number of sequences.
!
!    Input, integer ( kind = 4 ) CONLEN, the consensus length.
!
!    Output, real DISTANCE(SEQUENCE_MAX,SEQUENCE_MAX), the distance matrix.
!
!    Input, integer ( kind = 4 ) SEQUENCE_INDEX(SEQUENCE_MAX,SEQUENCE_LENGTH_MAX), ?
!
!    Input, real AAMAT(0:20,0:20), the matching score matrix.
!
  implicit none

  integer ( kind = 4 ) sequence_length_max
  integer ( kind = 4 ) sequence_max

  real aamat(0:20,0:20)
  integer ( kind = 4 ) conlen
  real dist
  real distance(sequence_max,sequence_max)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ik
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jk
  integer ( kind = 4 ) k
  integer ( kind = 4 ) sequence_index(sequence_max,sequence_length_max)
  integer ( kind = 4 ) sequence_number

  do i = 1, sequence_number
    do j = i+1, sequence_number

      dist = 0.0E+00
  
      do k = 1, conlen

        ik = sequence_index(i,k)
        jk = sequence_index(j,k)

        if ( ik /= 0 .and. jk /= 0 ) then 
          dist = dist + aamat(ik,jk)         
        end if

      end do

      dist = sqrt ( dist )
      distance(i,j) = dist
      distance(j,i) = dist

    end do
  end do
!
!  The diagonal of the matrix is zero.
!
  do i = 1, sequence_number
    distance(i,i) = 0.0E+00
  end do

  return
end
subroutine dist_mat_write ( output_unit, sequence_max, sequence_number, &
  distance )

!*****************************************************************************80
!
!! DIST_MAT_WRITE writes the distance matrix to a file.
!
!  Author:
!
!    Des Higgins, 
!    EMBL Data Library, 
!    April 1992
!
!  Modified:
!
!    06 April 2000
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OUTPUT_UNIT, the output unit number
!
!    Input, integer ( kind = 4 ) SEQUENCE_MAX, the maximum number of sequences.
!
!    Input, integer ( kind = 4 ) SEQUENCE_NUMBER, the number of sequences.
!
!    Input, real DISTANCE(SEQUENCE_MAX,SEQUENCE_MAX), the distance matrix.
!
  implicit none

  integer ( kind = 4 ) sequence_max

  real distance(sequence_max,sequence_max)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) sequence_number
  integer ( kind = 4 ) output_unit

  write ( output_unit, * ) sequence_number

  do i = 1, sequence_number
    write ( output_unit, '(5f12.4)' ) distance(i,1:i)
  end do

  return
end
subroutine read_sequence ( input_unit, sequence_max, sequence_length_max, &
  sequence_number, sequence, sequence_name, sequence_length, conlen )

!*****************************************************************************80
!
!! READ_SEQUENCE reads a FASTA or PIR format multiple alignment file.
!
!  Author:
!
!    Des Higgins, 
!    EMBL Data Library, 
!    April 1992
!
!  Modified:
!
!    06 April 2000
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) INPUT_UNIT, the input unit number.
!
!    Input, integer ( kind = 4 ) SEQUENCE_MAX, the maximum number of sequences.
!
!    Input, integer ( kind = 4 ) SEQUENCE_LENGTH_MAX, the maximum length of a sequence.
!
!    Output, integer ( kind = 4 ) SEQUENCE_NUMBER, the number of sequences read.
!
!    Output, character SEQUENCE(SEQUENCE_MAX,SEQUENCE_LENGTH_MAX), ?
!
!    Output, character ( len = 16 ) SEQUENCE_NAME(SEQUENCE_MAX), ?
!
!    Output, integer ( kind = 4 ) SEQUENCE_LENGTH(SEQUENCE_MAX), ?
!
!    Output, integer ( kind = 4 ) CONLEN, the consensus length.
!
  implicit none

  integer ( kind = 4 ), parameter :: inlen = 250
  integer ( kind = 4 ) sequence_length_max
  integer ( kind = 4 ) sequence_max

  character ch
  integer ( kind = 4 ) conlen
  character ( len = 40 ) :: digits = '0123456789.,;:\\|/?[]{}_=+)(&^%$#@!`~'
  integer ( kind = 4 ) i
  integer ( kind = 4 ) input_unit
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) j
  character ( len = inlen ) line
  logical pir_format
  character sequence(sequence_max,sequence_length_max)
  integer ( kind = 4 ) sequence_length(sequence_max)
  character ( len = 16 ) sequence_name(sequence_max)
  integer ( kind = 4 ) sequence_number

  read ( input_unit, '(a)', iostat = ios ) line

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'READ_SEQUENCE - Fatal error!'
    write ( *, '(a)' ) '  Error while reading input from file.'
    stop
  end if
!
!  Check for a > as first character.
!
  if ( line(1:1) /= '>' ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'READ_SEQUENCE - Fatal error!'
    write ( *, '(a)' ) '  Sequences must be in FASTA or PIR format.'
    write ( *, '(a)' ) '  Sequence names must start with a > character.'
    write ( *, '(a)' ) '  EXIT'
    stop
  end if
!
!  Check for FASTA or PIR format.
!
  if ( line(4:4) == ';' ) then
    pir_format = .TRUE.
  else
    pir_format = .FALSE.
  end if

  backspace ( input_unit )

  sequence_number = 0

1 continue

  read ( input_unit, '(a)', end = 90 ) line
  sequence_number = sequence_number + 1

  if ( sequence_number > sequence_max ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'READ_SEQUENCE - Fatal error!'
    write ( *, '(a)' ) '  Maximum number of sequences exceeded.'
    write ( *, '(a,i6)' ) '  Maximum number is ', sequence_max
    stop
  end if

  sequence_length(sequence_number) = 0

  if ( pir_format ) then
    sequence_name(sequence_number) = line(5:20)
  else
    sequence_name(sequence_number) = line(2:17)
  end if
!
!  Skip an extra line if PIR format is used.
!
  if ( pir_format ) then
    read ( input_unit, '(a)' ) line
  end if

  do

    read ( input_unit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      exit
    end if
!
!  Look for the start of the next sequence.
!
    if ( line(1:1) == '>' ) then
      backspace ( input_unit )
      go to 1
    end if

    do i = 1, inlen

      ch = line(i:i)

      if ( ch /= ' ' ) then

        if ( index ( digits, ch ) == 0 ) then
!
!  Look for a * if PIR format is used (= end of sequence)
!
          if ( pir_format .and. ch == '*' ) then

3           continue

            read ( input_unit, '(a)', iostat = ios ) line

            if ( ios /= 0 ) then
              go to 90
            end if
!
!  Look for the start of the next sequence
!
            if ( line(1:1) /= '>' ) then
              go to 3
            end if

            backspace ( input_unit )
            go to 1

          end if
!
!  Convert to upper case.
!
          if ( ichar ( ch ) >= 97 ) then
            ch = char ( ichar ( ch ) - 32 )
          end if

          sequence_length(sequence_number) = &
            sequence_length(sequence_number) + 1

          if ( sequence_length(sequence_number) > sequence_length_max ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'READ_SEQUENCE - Fatal error!'
            write ( *, '(a)' ) '  Maximum sequence length exceeded.'
            write ( *, '(a,i6)' ) '  Maximum length is ', sequence_length_max
            stop
          end if

          sequence(sequence_number,sequence_length(sequence_number)) = ch       

        end if

      end if

    end do
!
!  Go back to read more sequences.
!
  end do

90    continue

  conlen = sequence_length(1)

  do i = 1, sequence_number

    if ( sequence_length(i) /= conlen ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'READ_SEQUENCE - Fatal error!'
      write ( *, '(a)' ) '  Sequences MUST be aligned already!!'
      write ( *, '(a,i6)' ) '  Sequence number: ', i
      write ( *, '(a,i6)' ) '  Length = ', sequence_length(i)
      stop
    end if

    conlen = max ( conlen, sequence_length(i) )

  end do
!
!  Count the number of non-blanks in each sequence.
!
  do i = 1, sequence_number
    sequence_length(i) = 0
    do j = 1, conlen
      if ( sequence(i,j) /= '-' ) then
        sequence_length(i) = sequence_length(i) + 1
      end if
    end do
  end do

  return
end
subroutine set_aamat ( aamat, aacids )

!*****************************************************************************80
!
!! SET_AAMAT sets up the amino acid distance matrix.
!
!  Author:
!
!    Des Higgins, 
!    EMBL Data Library, 
!    April 1992
!
!  Modified:
!
!    06 April 2000
!
!  Parameters:
!
!    Output, real AAMAT(0:20,0:20), the matching score matrix.
!
!    Input, character ( len = 20 ) AACIDS, the one letter codes for
!    the amino acids.
!
  implicit none

  character ( len = 20 ) aacids
  real aamat(0:20,0:20)
  integer ( kind = 4 ) i
  real idmatch
  integer ( kind = 4 ) j
  real match1
  real match2
  real nomatch

  idmatch = 0.0E+00
  match1  = 1.0E+00
  match2  = 2.0E+00
  nomatch = 3.0E+00
!
!  Initialize the matrix to no matching.
!
  aamat(0:20,0:20) = nomatch
!
!  Set groups that are somewhat close.
!
  call set_aascore ( 'DEKRHNQST', match2, aamat, aacids )
  call set_aascore ( 'ILVFWYCM',  match2, aamat, aacids )
!
!  Set subgroups that are very close.
!
  call set_aascore ( 'DE',        match1, aamat, aacids )
  call set_aascore ( 'KRH',       match1, aamat, aacids )
  call set_aascore ( 'NQ',        match1, aamat, aacids )
  call set_aascore ( 'ST',        match1, aamat, aacids )
  call set_aascore ( 'ILV',       match1, aamat, aacids )
  call set_aascore ( 'FWY',       match1, aamat, aacids )
  call set_aascore ( 'AG',        match1, aamat, aacids )
!
!  Set the score for the identity match.
!
  do i = 1, 20
    aamat(i,i) = idmatch
  end do

  return
end
subroutine set_aascore ( group, match, aamat, aacids )

!*****************************************************************************80
!
!! SET_AASCORE updates a match matrix with a group of matching data.
!
!  Author:
!
!    Des Higgins, 
!    EMBL Data Library, 
!    April 1992
!
!  Modified:
!
!    06 April 2000
!
!  Parameters:
!
!    Input, character ( len = * ) GROUP, a string of characters, each of
!    which is a member of a group.
!
!    Input, real MATCH, the matching score to be assigned for matching
!    any two members of the given group.
!
!    Input/output, real AAMAT(0:20,0:20), the matching score matrix.
!
!    Input, character ( len = 20 ) AACIDS, the one letter codes for
!    the amino acids.
!
  implicit none

  character ( len = 20 ) aacids
  real aamat(0:20,0:20)
  integer ( kind = 4 ) col
  integer ( kind = 4 ) grlen
  character ( len = * ) group
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real match
  integer ( kind = 4 ) row

  grlen = len ( group )
 
  do i = 1, grlen
    do j = 1, grlen
      col = index ( aacids, group(i:i) )
      row = index ( aacids, group(j:j) )
      aamat(col,row) = match
      aamat(row,col) = match
    end do
  end do
  
  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    May 31 2001   9:45:54.872 AM
!
!  Modified:
!
!    31 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  character ( len = 8 ) date
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  character ( len = 10 )  time
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y
  character ( len = 5 ) zone

  call date_and_time ( date, time, zone, values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
subroutine toss_gaps ( sequence_max, sequence_length_max, sequence_number, &
  sequence, conlen )

!*****************************************************************************80
!
!! TOSS_GAPS obliterates information in a position where any sequence has a gap.
!
!  Discussion:
!
!    If any sequence has a gap, indicated by the character '-', then
!    the corresponding element of all sequences is set to '-'.
!
!  Author:
!
!    Des Higgins, 
!    EMBL Data Library, 
!    April 1992
!
!  Modified:
!
!    06 April 2000
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) SEQUENCE_MAX, the maximum number of sequences.
!
!    Input, integer ( kind = 4 ) SEQUENCE_LENGTH_MAX, the maximum length of a sequence.
!
!    Input, integer ( kind = 4 ) SEQUENCE_NUMBER, the number of sequences.
!
!    Input/output, character SEQUENCE(SEQUENCE_MAX,SEQUENCE_LENGTH_MAX), ?
!
!    Input, integer ( kind = 4 ) CONLEN, the consensus length.
!
  implicit none

  integer ( kind = 4 ) sequence_length_max
  integer ( kind = 4 ) sequence_max

  integer ( kind = 4 ) conlen
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  character sequence(sequence_max,sequence_length_max)
  integer ( kind = 4 ) sequence_number

  do i = 1, conlen

    do j = 1, sequence_number

      if ( sequence(j,i) == '-' ) then
        sequence(1:sequence_number,i) = '-'
        exit
      end if

    end do

  end do

  return
end
