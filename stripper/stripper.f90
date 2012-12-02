program main

!*****************************************************************************80
!
!! MAIN is the main program for STRIPPER.
!
!  Discussion:
!
!    STRIPPER allows a user to interactively modify a file.
!
!  Discussion:
!
!    This program can do a lot of things, but you'd best go look
!    at the HELP routine to get a sketchy idea of what and how.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 August 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: lunit1 = 2
  integer ( kind = 4 ), parameter :: lunit2 = 3
  integer ( kind = 4 ), parameter :: mmax = 500

  logical back
  character ( len = 6 ) chunk
  character ( len = 2 ) combeg
  character ( len = 2 ) comend
  character ( len = 80 ) command
  logical commentout
  integer ( kind = 4 ) iblank
  integer ( kind = 4 ) icap
  integer ( kind = 4 ) icapf
  integer ( kind = 4 ) icapfc
  integer ( kind = 4 ) icolumn
  integer ( kind = 4 ) icomment
  integer ( kind = 4 ) icon
  character ( len = 40 ) cut
  character ( len = 80 ) input
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) ireph
  integer ( kind = 4 ) irepv
  integer ( kind = 4 ) join
  character ( len = 40 ) keep
  integer ( kind = 4 ) keep_max
  integer ( kind = 4 ) keep_min
  character ( len = 40 ) kill
  character ( len = 10 ) lang
  logical left
  logical lexist
  integer ( kind = 4 ) lwrap
  integer ( kind = 4 ) margel
  integer ( kind = 4 ) marger
  integer ( kind = 4 ) mreci
  integer ( kind = 4 ) mreco
  integer ( kind = 4 ) nchunk
  integer ( kind = 4 ) number
  character ( len = 80 ) output
  logical page
  logical pause
  logical dorot13
  logical s_eqi
  logical shocon
  character ( len = 10 ) showme

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'STRIPPER'
  write ( *, '(a)' ) '  FORTRAN90 version'
!
!  Shut FTNCHEK's whining.
!
  back = .false.
  commentout = .false.
  lang = ' '
!
!  Say hello.
!
  call hello

  command = 'DEFAULTS'

  do
!
!  No command.
!
  if ( command == ' ' ) then
!
!  BACK command
!
  else if ( s_eqi ( command, 'BACK' ) ) then

    back = .not. back

    showme = 'back'
!
!  BREAK command
!
  else if ( s_eqi ( command, 'BREAK' ) ) then

    join = -1

    showme = 'break'
!
!  CHUNK =
!
  else if ( s_eqi ( command(1:5), 'CHUNK' ) ) then

    call s_blank_delete ( command )

    if ( s_eqi ( command(1:6), 'CHUNK=' ) ) then
      chunk = command(7:)
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'Enter chunking (FNAME, LINE, LINEP, or WORD):'
      read ( *, '(a)', iostat = ios ) chunk
      if ( ios /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STRIPPER - Fatal error!'
        write ( *, '(a)' ) '  Error reading user input.'
        stop
      end if
    end if

    showme = 'chunk'
!
!  COMBEG=
!
  else if ( s_eqi ( command(1:6), 'COMBEG' ) ) then

    call s_blank_delete ( command )

    if ( s_eqi ( command(1:7), 'COMBEG=' ) ) then
      combeg = command(8:)
    else
      write ( *, '(a)' ) 'Enter "comment begin" characters:'
      read ( *, '(a)', iostat = ios ) combeg
      if ( ios /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STRIPPER - Fatal error!'
        write ( *, '(a)' ) '  Error reading user input.'
        stop
      end if
    end if

    showme = 'COMBEG'
!
!  COMEND=
!
  else if ( s_eqi ( command(1:6), 'COMEND' ) ) then

    call s_blank_delete ( command )

    if ( s_eqi ( command(1:7), 'COMEND=' ) ) then
      comend = command(8:)
    else
      write ( *, '(a)' ) 'Enter "comment end" characters:'
      read ( *, '(a)', iostat = ios ) comend
      if ( ios /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STRIPPER - Fatal error!'
        write ( *, '(a)' ) '  Error reading user input.'
        stop
      end if
    end if

    showme = 'COMEND'
!
!  COMMENTOUT
!
  else if ( s_eqi ( command, 'COMMENTOUT' ) ) then

    commentout = .not. commentout

    showme = 'COMMENTOUT'
!
!  CUT=
!
  else if ( s_eqi ( command, 'CUT' ) .or. &
            s_eqi ( command(1:4), 'CUT=' ) ) then

    call s_blank_delete ( command )

    if ( s_eqi ( command(1:4), 'CUT=' ) ) then
      cut = command(5:)
    else
      write ( *, '(a)' ) 'Enter the CUT string:'
      read ( *, '(a)', iostat = ios ) cut
      if ( ios /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STRIPPER - Fatal error!'
        write ( *, '(a)' ) '  Error reading user input.'
        stop
      end if
    end if

    showme = 'CUT'
!
!  DEFAULTS
!
  else if ( s_eqi ( command(1:3), 'DEF' ) ) then

    call init ( back, chunk, combeg, comend, commentout, cut, iblank, icap, &
      icapf, icapfc, icolumn, icomment, icon, input, ireph, irepv, join, keep, &
      keep_max, keep_min, kill, lang, left, lwrap, margel, marger, &
      mmax, mreci, mreco, nchunk, number, output, page, pause, dorot13, shocon )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'All options set to default values.'

    showme = 'NONE'
!
!  GO
!
  else if ( s_eqi ( command(1:2), 'GO' ) ) then

    call stripit ( back, chunk, combeg, comend, commentout, cut, iblank, &
      icap, icapf, icapfc, icolumn, icomment, icon, input, ireph, irepv, join, &
      keep, keep_max, keep_min, kill, lang, left, lunit1, lunit2, lwrap, &
      margel, marger, mreci, mreco, nchunk, number, output, page, pause, &
      dorot13, shocon )

    showme = 'NONE'
!
!  HELP
!
  else if ( s_eqi ( command(1:1), 'H' ) ) then

    call help

    showme = 'NONE'
!
!  IBLANK command
!
  else if ( s_eqi ( command(1:6), 'IBLANK' ) ) then

    call s_blank_delete ( command )

    if ( s_eqi ( command(1:7), 'IBLANK=' ) ) then
      read ( command(8:), *, iostat = ios ) iblank
      if ( ios /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STRIPPER - Fatal error!'
        write ( *, '(a)' ) '  Error reading user input.'
        stop
      end if
    else
      write ( *, '(a)' ) 'Enter IBLANK = 0, 1 or 2:'
      read ( *, *, iostat = ios ) iblank
      if ( ios /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STRIPPER - Fatal error!'
        write ( *, '(a)' ) '  Error reading user input.'
        stop
      end if
    end if

    showme = 'IBLANK'
!
!  ICAP=
!
  else if ( s_eqi ( command(1:4), 'ICAP' ) .and. &
      .not. s_eqi ( command(1:5), 'ICAPF' ) ) then

    call s_blank_delete ( command )

    if ( s_eqi ( command(1:5), 'ICAP=' ) ) then
      read ( command(6:), *, iostat = ios ) icap
      if ( ios /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STRIPPER - Fatal error!'
        write ( *, '(a)' ) '  Error reading user input.'
        stop
      end if
    else
      write ( *, '(a)' ) 'Enter ICAP option:'
      write ( *, '(a)' ) '-1: Lowercase each character;'
      write ( *, '(a)' ) ' 0: Leave alone;'
      write ( *, '(a)' ) '+1: Uppercase each character;'
      write ( *, '(a)' ) '+2: Uppercase each word;'
      write ( *, '(a)' ) '+3: Uppercase each line.'
      read ( *, *, iostat = ios ) icap
      if ( ios /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STRIPPER - Fatal error!'
        write ( *, '(a)' ) '  Error reading user input.'
        stop
      end if
    end if

    showme = 'ICAP'
!
!  ICAPF=
!
  else if ( s_eqi ( command(1:5), 'ICAPF' ) .and. &
      .not. s_eqi ( command(1:6), 'ICAPFC' ) ) then

    call s_blank_delete ( command )

    if ( s_eqi ( command(1:6), 'ICAPF=' ) ) then
      read ( command(7:), *, iostat = ios ) icapf
      if ( ios /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STRIPPER - Fatal error!'
        write ( *, '(a)' ) '  Error reading user input.'
        stop
      end if
    else
      write ( *, '(a)' ) 'Enter ICAPF option for ' // trim ( lang ) &
        // ' statements:'
      write ( *, '(a)' ) '-1: Lowercase each character;'
      write ( *, '(a)' ) ' 0: Leave alone;'
      write ( *, '(a)' ) '+1: Uppercase each character;'
      write ( *, '(a)' ) '+2: Uppercase each word.'
      write ( *, '(a)' ) '+3: Uppercase each line.'
      read ( *, *, iostat = ios ) icapf
      if ( ios /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STRIPPER - Fatal error!'
        write ( *, '(a)' ) '  Error reading user input.'
        stop
      end if
    end if

    showme = 'ICAPF'
!
!  ICAPFC=
!
  else if ( s_eqi ( command(1:6), 'ICAPFC' ) ) then

    call s_blank_delete ( command )

    if ( s_eqi ( command(1:7), 'ICAPFC=' ) ) then
      read ( command(8:), *, iostat = ios ) icapfc
      if ( ios /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STRIPPER - Fatal error!'
        write ( *, '(a)' ) '  Error reading user input.'
        stop
      end if
    else
      write ( *, '(a)' ) 'Enter ICAPFC option for ' // trim ( lang ) &
        // ' comments:'
      write ( *, '(a)' ) '-1: Lowercase each character;'
      write ( *, '(a)' ) ' 0: Leave alone;'
      write ( *, '(a)' ) '+1: Uppercase each character;'
      write ( *, '(a)' ) '+2: Uppercase each word.'
      write ( *, '(a)' ) '+3: Uppercase each line.'
      read ( *, *, iostat = ios ) icapfc
      if ( ios /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STRIPPER - Fatal error!'
        write ( *, '(a)' ) '  Error reading user input.'
        stop
      end if
    end if

    showme = 'ICAPFC'
!
!  ICOLUMN =
!
  else if ( s_eqi ( command(1:7), 'ICOLUMN' ) ) then

    call s_blank_delete ( command )

    if ( s_eqi ( command(1:8), 'ICOLUMN=' ) ) then
      read ( command(9:), *, iostat = ios ) icolumn
      if ( ios /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STRIPPER - Fatal error!'
        write ( *, '(a)' ) '  Error reading user input.'
        stop
      end if
    else
      write ( *, '(a)' ) 'Enter column value, or 0:'
      read ( *, *, iostat = ios ) icolumn
      if ( ios /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STRIPPER - Fatal error!'
        write ( *, '(a)' ) '  Error reading user input.'
        stop
      end if
    end if

    showme = 'ICOLUMN'
!
!  ICOMMENT = (0 no action, 1 delete comments, 2 delete noncomments)
!
  else if ( s_eqi ( command(1:8), 'ICOMMENT' ) ) then

    if ( s_eqi ( command(1:9), 'ICOMMENT = ' ) ) then
      read ( command(10:), *, iostat = ios ) icomment
      if ( ios /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STRIPPER - Fatal error!'
        write ( *, '(a)' ) '  Error reading user input.'
        stop
      end if
    else
      write ( *, '(a)' ) 'Enter ICOMMENT = 0/1/2:'
      read ( *, *, iostat = ios ) icomment
      if ( ios /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STRIPPER - Fatal error!'
        write ( *, '(a)' ) '  Error reading user input.'
        stop
      end if
    end if

    showme = 'ICOMMENT'
!
!  ICON=
!
  else if ( s_eqi ( command(1:4), 'ICON' ) ) then

    call s_blank_delete ( command )

    if ( s_eqi ( command(1:5), 'ICON=' ) ) then
      read ( command(6:), *, iostat = ios ) icon
      if ( ios /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STRIPPER - Fatal error!'
        write ( *, '(a)' ) '  Error reading user input.'
        stop
      end if
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'Enter control character option ICON:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '-1: Symbols become control characters.'
      write ( *, '(a)' ) ' 0: Control characters will be preserved.'
      write ( *, '(a)' ) ' 1: Control characters become symbols.'
      write ( *, '(a)' ) ' 2: Control characters replaced by blanks.'
      write ( *, '(a)' ) ' 3: Control characters removed.'

      read ( *, *, iostat = ios ) icon
      if ( ios /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STRIPPER - Fatal error!'
        write ( *, '(a)' ) '  Error reading user input.'
        stop
      end if

    end if

    showme = 'ICON'
!
!  INPUT = filename
!    or
!  < filename
!
  else if ( s_eqi ( command(1:5), 'INPUT' ) .or. command(1:1) == '<' ) then

    call s_blank_delete ( command )

    if ( s_eqi ( command(1:6), 'INPUT=' ) ) then
      input = command(7:)
    else if ( command(1:1) == '<' ) then
      input = command(2:)
    else
      write ( *, '(a)' ) 'Enter the name of the input file:'
      read ( *, '(a)', iostat = ios ) input
      if ( ios /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STRIPPER - Fatal error!'
        write ( *, '(a)' ) '  Error reading user input.'
        stop
      end if
    end if

    showme = 'INPUT'

    if ( input /= '*' ) then

      inquire ( file = input, exist = lexist )

      if ( .not. lexist ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STRIPPER - Fatal error!'
        write ( *, '(a)' ) '  The input file does not exist!'
        showme = 'NONE'
      end if

    end if
!
!  IREPH=
!
  else if ( s_eqi ( command(1:5), 'IREPH' ) ) then

    call s_blank_delete ( command )

    if ( s_eqi ( command(1:6), 'IREPH=' ) ) then
      read ( command(7:), *, iostat = ios ) ireph
      if ( ios /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STRIPPER - Fatal error!'
        write ( *, '(a)' ) '  Error reading user input.'
        stop
      end if
    else
      write ( *, '(a)' ) 'Enter horizontal line repeat factor.'
      read ( *, *, iostat = ios ) ireph
      if ( ios /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STRIPPER - Fatal error!'
        write ( *, '(a)' ) '  Error reading user input.'
        stop
      end if
    end if

    showme = 'IREPH'
!
!  IREPV=
!
  else if ( s_eqi ( command(1:5), 'IREPV' ) ) then

    call s_blank_delete ( command )

    if ( s_eqi ( command(1:6), 'IREPV=' ) ) then
      read ( command(7:), *, iostat = ios ) irepv
      if ( ios /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STRIPPER - Fatal error!'
        write ( *, '(a)' ) '  Error reading user input.'
        stop
      end if
    else
      write ( *, '(a)' ) 'Enter vertical line repeat factor.'
      read ( *, *, iostat = ios ) irepv
      if ( ios /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STRIPPER - Fatal error!'
        write ( *, '(a)' ) '  Error reading user input.'
        stop
      end if
    end if

    showme = 'IREPV'
!
!  JOIN command
!
  else if ( s_eqi ( command, 'JOIN' ) ) then

    join = +1

    showme = 'join'
!
!  KEEP=
!
  else if ( s_eqi ( command, 'KEEP' ) .or. &
            s_eqi ( command(1:5), 'KEEP=' ) ) then

    call s_blank_delete ( command )

    if ( s_eqi ( command(1:5), 'KEEP=' ) ) then
      keep = command(6:)
    else
      write ( *, '(a)' ) 'Enter the KEEP string:'
      write ( *, '(a)' ) '  (Use "<" or ">" to specify the string must'
      write ( *, '(a)' ) '  begin in the first column, or end in the last.)'
      write ( *, '(a)' ) '  (Use ANY_ALPHA to keep lines with any'
      write ( *, '(a)' ) '  alphabetic characters)'
      read ( *, '(a)', iostat = ios ) keep
      if ( ios /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STRIPPER - Fatal error!'
        write ( *, '(a)' ) '  Error reading user input.'
        stop
      end if
    end if

    showme = 'KEEP'
!
!  KEEP_MAX=
!
  else if ( s_eqi ( command(1:8), 'KEEP_MAX' ) ) then

    call s_blank_delete ( command )

    if ( command(9:9) == '=' ) then
      read ( command(10:), *, iostat = ios ) keep_max
      if ( ios /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STRIPPER - Fatal error!'
        write ( *, '(a)' ) '  Error reading user input.'
        stop
      end if
    else
      write ( *, '(a)' ) 'Enter the maximum KEEP length:'
      read ( *, *, iostat = ios ) keep_max
      if ( ios /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STRIPPER - Fatal error!'
        write ( *, '(a)' ) '  Error reading user input.'
        stop
      end if
    end if

    showme = 'KEEP_MAX'
!
!  KEEP_MIN=
!
  else if ( s_eqi ( command(1:8), 'KEEP_MIN' ) ) then

    call s_blank_delete ( command )

    if ( command(9:9) == '=' ) then
      read ( command(10:), *, iostat = ios ) keep_min
      if ( ios /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STRIPPER - Fatal error!'
        write ( *, '(a)' ) '  Error reading user input.'
        stop
      end if
    else
      write ( *, '(a)' ) 'Enter the minimum KEEP length:'
      read ( *, *, iostat = ios ) keep_min
      if ( ios /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STRIPPER - Fatal error!'
        write ( *, '(a)' ) '  Error reading user input.'
        stop
      end if
    end if

    showme = 'KEEP_MIN'
!
!  KILL=
!
  else if ( s_eqi ( command(1:4), 'KILL' ) ) then

    call s_blank_delete ( command )

    if ( s_eqi ( command(1:5), 'KILL=' ) ) then
      kill = command(6:)
    else
      write ( *, '(a)' ) 'Enter the KILL string:'
      write ( *, '(a)' ) '  (Use "<" to begin in the first column, '
      write ( *, '(a)' ) '  or ">" to end in the last.)'
      read ( *, '(a)', iostat = ios ) kill
      if ( ios /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STRIPPER - Fatal error!'
        write ( *, '(a)' ) '  Error reading user input.'
        stop
      end if
    end if

    showme = 'KILL'
!
!  LANG=
!
  else if ( s_eqi ( command(1:4), 'LANG' ) ) then

    call s_blank_delete ( command )

    if ( s_eqi ( command(1:5), 'LANG=' ) ) then
      lang = command(6:)
    else if ( s_eqi ( command(1:9), 'LANGUAGE=' ) ) then
      lang = command(10:)
    else
      write ( *, '(a)' ) 'Enter the language'
      write ( *, '(a)' ) '  (ADA/C/C++/F77/F90/TEXT/UNIX):'
      read ( *, '(a)', iostat = ios ) lang
      if ( ios /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STRIPPER - Fatal error!'
        write ( *, '(a)' ) '  Error reading user input.'
        stop
      end if
    end if

  showme = 'LANG'

  if ( s_eqi ( lang, 'ADA' ) ) then
    combeg = '--'
    comend = ' '
  else if ( s_eqi ( lang, 'C' ) ) then
    combeg = '/*'
    comend = '*/'
  else if ( s_eqi ( lang, 'C++' ) ) then
    combeg = '//'
    comend = ' '
  else if ( s_eqi ( lang, 'F77' ) ) then
    combeg = '!'
    comend = ' '
  else if ( s_eqi ( lang, 'F90' ) ) then
    combeg = '!'
    comend = ' '
  else if ( s_eqi ( lang, 'TEXT' ) ) then
    combeg = ' '
    comend = ' '
  else if ( s_eqi ( lang, 'UNIX' ) ) then
    combeg = '#'
    comend = ' '
  end if
!
!  LEFT
!
  else if ( s_eqi ( command, 'LEFT' ) ) then

    left = .not. left

    showme = 'LEFT'
!
!  LWRAP=
!
  else if ( s_eqi ( command(1:5), 'LWRAP' ) ) then

    call s_blank_delete ( command )

    if ( s_eqi ( command(1:6), 'LWRAP=' ) ) then
      read ( command(7:), *, iostat = ios ) lwrap
      if ( ios /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STRIPPER - Fatal error!'
        write ( *, '(a)' ) '  Error reading user input.'
        stop
      end if
    else
      write ( *, '(a)' ) 'Enter line wrapping length'
      read ( *, *, iostat = ios ) lwrap
      if ( ios /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STRIPPER - Fatal error!'
        write ( *, '(a)' ) '  Error reading user input.'
        stop
      end if
    end if

    showme = 'LWRAP'
!
!  MARGEL =
!
  else if ( s_eqi ( command(1:6), 'MARGEL' ) ) then

    call s_blank_delete ( command )

    if ( s_eqi ( command(1:7), 'MARGEL=' ) ) then
      read ( command(8:), *, iostat = ios ) margel
      if ( ios /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STRIPPER - Fatal error!'
        write ( *, '(a)' ) '  Error reading user input.'
        stop
      end if
    else
      write ( *, '(a)' ) 'Enter left margin:'
      read ( *, *, iostat = ios ) margel
      if ( ios /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STRIPPER - Fatal error!'
        write ( *, '(a)' ) '  Error reading user input.'
        stop
      end if
    end if

    showme = 'MARGEL'
!
!  MARGER =
!
  else if ( s_eqi ( command(1:6), 'MARGER' ) ) then

    call s_blank_delete ( command )

    if ( s_eqi ( command(1:7), 'MARGER=' ) ) then
      read ( command(8:), *, iostat = ios ) marger
      if ( ios /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STRIPPER - Fatal error!'
        write ( *, '(a)' ) '  Error reading user input.'
        stop
      end if
    else
      write ( *, '(a)' ) 'Enter right margin:'
      read ( *, *, iostat = ios ) marger
      if ( ios /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STRIPPER - Fatal error!'
        write ( *, '(a)' ) '  Error reading user input.'
        stop
      end if
    end if

    if ( marger > mmax ) then
      write ( *, '(a)' ) 'STRIPPER - Error!'
      write ( *, '(a,i6)' ) '  The maximum column cannot be greater than', mmax
      marger = mmax
    end if

    showme = 'MARGER'
!
!  MRECI =
!
  else if ( s_eqi ( command(1:5), 'MRECI' ) ) then

    call s_blank_delete ( command )

    if ( s_eqi ( command(1:6), 'MRECI=' ) ) then
      read ( command(7:), *, iostat = ios ) mreci
      if ( ios /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STRIPPER - Fatal error!'
        write ( *, '(a)' ) '  Error reading user input.'
        stop
      end if
    else
      write ( *, '(a)' ) 'Enter maximum number of input records'
      read ( *, *, iostat = ios ) mreci
      if ( ios /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STRIPPER - Fatal error!'
        write ( *, '(a)' ) '  Error reading user input.'
        stop
      end if
    end if

    showme = 'MRECI'
!
!  MRECO =
!
  else if ( s_eqi ( command(1:5), 'MRECO' ) ) then

    call s_blank_delete ( command )

    if ( s_eqi ( command(1:6), 'MRECO=' ) ) then
      read ( command(7:), *, iostat = ios ) mreco
      if ( ios /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STRIPPER - Fatal error!'
        write ( *, '(a)' ) '  Error reading user input.'
        stop
      end if
    else
      write ( *, '(a)' ) 'Enter maximum number of output records'
      read ( *, *, iostat = ios ) mreco
      if ( ios /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STRIPPER - Fatal error!'
        write ( *, '(a)' ) '  Error reading user input.'
        stop
      end if
    end if

    showme = 'MRECO'
!
!  NCHUNK=
!
  else if ( s_eqi ( command(1:6), 'NCHUNK' ) ) then

    call s_blank_delete ( command )

    if ( s_eqi ( command(1:7), 'NCHUNK=' ) ) then
      read ( command(8:), *, iostat = ios ) nchunk
      if ( ios /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STRIPPER - Fatal error!'
        write ( *, '(a)' ) '  Error reading user input.'
        stop
      end if
    else
      write ( *, '(a)' ) 'Enter number of characters to read.'
      read ( *, *, iostat = ios )  nchunk
      if ( ios /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STRIPPER - Fatal error!'
        write ( *, '(a)' ) '  Error reading user input.'
        stop
      end if
    end if

    if ( nchunk < 1 ) then
      nchunk = 80
      write ( *, '(a)' ) 'NCHUNK cannot be less than 1!'
    end if

    showme = 'NCHUNK'
!
!  NUMBER command
!
  else if ( s_eqi ( command(1:6), 'NUMBER' ) ) then

    call s_blank_delete ( command )

    if ( s_eqi ( command(1:7), 'NUMBER=' ) ) then
      read ( command(8:), *, iostat = ios )  number
      if ( ios /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STRIPPER - Fatal error!'
        write ( *, '(a)' ) '  Error reading user input.'
        stop
      end if
    else
      write ( *, '(a)' ) 'Enter number option, -1, 0, +1:'
      read ( *, *, iostat = ios )  number
      if ( ios /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STRIPPER - Fatal error!'
        write ( *, '(a)' ) '  Error reading user input.'
        stop
      end if
    end if

    showme = 'NUMBER'
!
!  OUTPUT = filename
!  or
!  > filename
!
  else if ( s_eqi ( command(1:6), 'OUTPUT' ) .or. command(1:1) == '>' ) then

    call s_blank_delete ( command )

    if ( s_eqi ( command(1:7), 'OUTPUT=' ) ) then
      output = command(8:)
    else if ( command(1:1) == '>' ) then
      output = command(2:)
    else
      write ( *, '(a)' ) 'Enter the name of the output file:'
      read ( *, '(a)', iostat = ios ) output
      if ( ios /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STRIPPER - Fatal error!'
        write ( *, '(a)' ) '  Error reading user input.'
        stop
      end if
    end if

    showme = 'OUTPUT'
!
!  PAGE command
!
  else if ( s_eqi ( command, 'PAGE' ) ) then

    page = .not. page

    showme = 'PAGE'
!
!  PAUSE command
!
  else if ( s_eqi ( command, 'PAUSE' ) ) then

    pause = .not. pause

    showme = 'PAUSE'
!
!  QUIT command
!
  else if ( s_eqi ( command(1:2), 'QY' ) .or. &
            s_eqi ( command, 'EXIT' ) .or. &
            s_eqi ( command, 'QUIT' ) .or. &
            s_eqi ( command, 'STOP' ) ) then

    showme = 'none'
    exit

  else if ( s_eqi ( command(1:1), 'Q' ) ) then

    showme = 'NONE'

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Enter "Y" to confirm you want to stop.'
    read ( *, '(a)', iostat = ios ) command
      if ( ios /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STRIPPER - Fatal error!'
        write ( *, '(a)' ) '  Error reading user input.'
        stop
      end if

    if ( s_eqi ( command(1:1), 'Y' ) ) then
      exit
    end if
!
!  ROT13
!
  else if ( s_eqi ( command, 'ROT13' ) ) then

    dorot13 = .not. dorot13

    showme = 'ROT13'
!
!  SHOCON
!
  else if ( s_eqi ( command, 'SHOCON' ) ) then

    shocon = .not. shocon

    showme = 'SHOCON'
!
!  SHOW command
!
  else if ( s_eqi ( command, 'SHOW' ) ) then

    showme = 'ALL'
!
!  Unrecognized command.
!
  else

    showme = 'NONE'

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'STRIPPER - Warning!'
    write ( *, '(a)' ) '  Unrecognized command: ' // trim ( command )

  end if
!
!  Show results of current command.
!
  call show ( back, chunk, combeg, comend, commentout, cut, iblank, icap, &
    icapf, icapfc, icolumn, icomment, icon, input, ireph, irepv, join, keep, &
    kill, lang, left, lwrap, margel, marger, mmax, mreci, mreco, &
    nchunk, number, output, page, pause, dorot13, shocon, showme )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Enter command (H for help)'
  read ( *, '(a)', iostat = ios ) command

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'STRIPPER - Fatal error!'
    write ( *, '(a)' ) '  Error reading user input.'
    stop
  end if

  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'STRIPPER:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine ch_cap ( c )

!*****************************************************************************80
!
!! CH_CAP capitalizes a single character.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character C, the character to capitalize.
!
  implicit none

  character c
  integer ( kind = 4 ) itemp

  itemp = ichar ( c )

  if ( 97 <= itemp .and. itemp <= 122 ) then
    c = char ( itemp - 32 )
  end if

  return
end
function ch_is_alpha ( c )

!*****************************************************************************80
!
!! CH_IS_ALPHA returns TRUE if C is an alphabetic character.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C, a character to check.
!
!    Output, logical CH_IS_ALPHA is TRUE if C is an alphabetic character.
!
  implicit none

  character c
  logical ch_is_alpha

  if ( ( lle ( 'a', c ) .and. lle ( c, 'z' ) ) .or. &
       ( lle ( 'A', c ) .and. lle ( c, 'Z' ) ) ) then
    ch_is_alpha = .true.
  else
    ch_is_alpha = .false.
  end if

  return
end
function ch_is_control ( c )

!*****************************************************************************80
!
!! CH_IS_CONTROL reports whether a character is a control character or not.
!
!  Discussion:
!
!    A "control character" has ASCII code <= 31 or ASCII code => 127.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C, the character to be tested.
!
!    Output, logical CH_IS_CONTROL, TRUE if C is a control character, and
!    FALSE otherwise.
!
  implicit none

  character c
  logical ch_is_control

  if ( ichar ( c ) <= 31 .or. ichar ( c ) >= 127 ) then
    ch_is_control = .true.
  else
    ch_is_control = .false.
  end if

  return
end
subroutine ch_low ( c )

!*****************************************************************************80
!
!! CH_LOW lowercases a single character.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character C, the character to be lowercased.
!
  implicit none

  character c
  integer ( kind = 4 ) itemp

  itemp = ichar ( c )

  if ( 65 <= itemp .and. itemp <= 90 ) then
    c = char ( itemp + 32 )
  end if

  return
end
function ch_to_rot13 ( c )

!*****************************************************************************80
!
!! CH_TO_ROT13 converts a character to its ROT13 equivalent.
!
!  Discussion:
!
!    Two applications of CH_TO_ROT13 to a character will return the original.
!
!  Example:
!
!    Input:  Output:
!
!    a       n
!    C       P
!    J       W
!    5       5
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C, the character to be converted.
!
!    Output, character CH_TO_ROT13, the ROT13 equivalent of the character.
!
  implicit none

  character c
  character ch_to_rot13
  integer ( kind = 4 ) itemp

  itemp = ichar ( c )

  if ( itemp >= 65 .and. itemp <= 77 ) then
    itemp = itemp + 13
  else if ( itemp >= 78 .and. itemp <= 90 ) then
    itemp = itemp - 13
  else if ( itemp >= 97 .and. itemp <= 109 ) then
    itemp = itemp + 13
  else if ( itemp >= 110 .and. itemp <= 122 ) then
    itemp = itemp - 13
  end if

  ch_to_rot13 = char ( itemp )

  return
end
subroutine ch_to_sym ( c, sym )

!*****************************************************************************80
!
!! CH_TO_SYM returns a printable symbol for any ASCII character.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C, the character to be represented.
!
!    Output, character ( len = 4 ) SYM, is the printable symbol for CHR.
!
  implicit none

  character c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iput
  character ( len = 4 ) sym

  i = ichar ( c )

  sym = ' '

  iput = 0
!
!  Characters 128-255 are symbolized with a ! prefix.
!  Then shift them down by 128.
!  Now all values of I are between 0 and 127.
!
  if ( i >= 128 ) then
    i = mod ( i, 128 )
    iput = iput + 1
    sym(iput:iput) = '!'
  end if
!
!  Characters 0-31 are symbolized with a ^ prefix.
!  Shift them up by 64.  Now all values of I are between 32 and 127.
!
  if ( i <= 31 ) then
    i = i + 64
    iput = iput + 1
    sym(iput:iput) = '^'
  end if
!
!  Characters 32 through 126 are themselves.
!
  if ( i <= 126 ) then
    iput = iput + 1
    sym(iput:iput) = char ( i )
!
!  Character 127 is DEL.
!
  else
    iput = iput + 1
    sym(iput:iput+2) = 'DEL'
  end if

  return
end
subroutine chra_to_s ( s1, s2 )

!*****************************************************************************80
!
!! CHRA_TO_S replaces control characters by printable symbols.
!
!  Table:
!
!    ICHAR(c)    Symbol
!    --------    ------
!      0          ^@
!      1          ^A
!    ...         ...
!     31          ^_
!     32         (space)
!    ...         ...
!    126         ~
!    127         DEL
!    128         !^@
!    ...         ...
!    159         !^_
!    160         !(space)
!    ...         ...
!    254         !~
!    255         !DEL
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S1, the string to be operated on.
!
!    Output, character ( len = * ) S2, a copy of S1, except that each
!    control character has been replaced by a symbol.
!
  implicit none

  logical ch_is_control
  integer ( kind = 4 ) iget
  integer ( kind = 4 ) iput
  integer ( kind = 4 ) lsym
  integer ( kind = 4 ) nchar1
  character ( len = * ) s1
  character ( len = * ) s2
  character ( len = 4 ) sym

  nchar1 = len_trim ( s1 )
  s2 = ' '

  iput = 1

  do iget = 1, nchar1

    if ( ch_is_control ( s1(iget:iget) ) ) then

      call ch_to_sym ( s1(iget:iget), sym )
      lsym = len_trim ( sym )

      s2(iput:iput+lsym-1) = sym(1:lsym)
      iput = iput + lsym

    else

      s2(iput:iput) = s1(iget:iget)
      iput = iput + 1

    end if

  end do

  return
end
subroutine chrdt6 ( line )

!*****************************************************************************80
!
!! CHRDT6 replaces TAB characters by 6 spaces.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) LINE, the line to be modified.  On
!    output, some significant characters at the end of LINE may have
!    been lost.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) iget
  integer ( kind = 4 ) iput
  integer ( kind = 4 ) lenc
  integer ( kind = 4 ) lens
  character ( len = * ) line
  integer ( kind = 4 ) ntab
  character, parameter :: TAB = char(9)
!
!  If no TAB's occur in the line, there is nothing to do.
!
  if ( index ( line, TAB ) == 0 ) then
    return
  end if
!
!  Otherwise, find out how long the line is.
!
  lenc = len_trim ( line )
  lens = len ( line )
!
!  Count the number of TAB's in the line.
!
  ntab = 0
  do i = 1, lenc
    if ( line(i:i) == TAB ) then
      ntab = ntab + 1
    end if
  end do
!
!  Now copy the line onto itself, going backwards.
!  As soon as we've processed the first TAB, we're done.
!
  iput = lenc + 5 * ntab

  do iget = lenc, 1, - 1

    if ( line(iget:iget) /= TAB ) then

      if ( iput <= lens ) then
        line(iput:iput) = line(iget:iget)
      end if

      iput = iput - 1

    else

      do i = iput, iput - 5, -1
        if ( i <= lens ) then
          line(i:i) = ' '
        end if
      end do

      iput = iput - 6
      ntab = ntab - 1

      if ( ntab == 0 ) then
        return
      end if

    end if

  end do

  return
end
subroutine chrs_to_a ( s1, s2 )

!*****************************************************************************80
!
!! CHRS_TO_A replaces all control symbols by control characters.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 August 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S1 is the string to be operated on.
!
!    Output, character ( len = * ) S2 is a copy of S1, except that each
!    control symbol has been replaced by a control character.
!
  implicit none

  character c
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) iput
  integer ( kind = 4 ) nchar1
  integer ( kind = 4 ) nchar2
  character ( len = * ) s1
  character ( len = * ) s2

  nchar1 = len_trim ( s1 )
  nchar2 = len ( s2 )

  ihi = 0
  iput = 0

  do

    if ( ihi >= nchar1 ) then
      return
    end if

    ilo = ihi + 1

    call sym_to_ch ( s1, c(ilo:), ihi )

    iput = iput + 1

    if ( iput > nchar2 ) then
      exit
    end if

    s2(iput:iput) = c

  end do

  return
end
subroutine doml ( line, margel, marger )

!*****************************************************************************80
!
!! DOML handles the left margin of the line of text.
!
!  Discussion:
!
!    DOML takes the current line of text, and "resets" it
!    so that the MARGEL-th character becomes the first character.
!
!    MARGEL may be positive, or nonpositive.  Nonpositive values
!    result in left-padding by blanks.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) LINE, the line to be reset.
!
!    Input, integer ( kind = 4 ) MARGEL, MARGER, the left and right "margins".
!    In effect, the string LINE(MARGEL:MARGER) will overwrite the
!    string LINE.  MARGEL may be zero or negative, in which case
!    a certain number of blank characters will be generated.
!
  implicit none

  character ( len = * ) line
  character ( len = 256 ) line2
  integer ( kind = 4 ) margel
  integer ( kind = 4 ) marger
!
!  Blank out everything to the right of the upper bound.
!
  if ( marger > 0 ) then
    line(marger+1:) = ' '
  end if
!
!  Either delete characters in positions 1 through MARGEL-1, or
!  insert characters in positions 1..1-MARGEL.
!
  if ( margel > 1 ) then
    line2 = line(margel:)
    line = line2
  else if ( margel < 1 ) then
    line2(1:1-margel) = ' '
    line2(2-margel:) = line
    line = line2
  end if

  return
end
subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is an integer between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is an integer between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IUNIT, the free unit number.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  logical lopen

  iunit = 0

  do i = 1, 99

    if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then

      inquire ( unit = i, opened = lopen, iostat = ios )

      if ( ios == 0 ) then
        if ( .not. lopen ) then
          iunit = i
          return
        end if
      end if

    end if

  end do

  return
end
subroutine getinp ( chunk, ierror, input, line, lunit1, mreci, nchunk, nreci )

!*****************************************************************************80
!
!! GETINP reads the next chunk of input from the input file.
!
!  Discussion:
!
!    Depending on the user's request, it reads
!    * a line,
!    * a given number of characters,
!    * a word, or
!    * a FORTRAN name.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) CHUNK, specifies how input is to be read.
!    'FNAME', means return the next FORTRAN token on each call;
!    'LINE', means return a line on each call;
!    'LINEP', means read NCHUNK characters, if possible.
!    'WORD', means return the next blank-separated word on each call;
!    'RETURN', means that the text of LINE should be 'returned' to the
!      input file, until called for later.
!
!    Output, integer ( kind = 4 ) IERROR, records whether an error occurred.
!
!    Input, character ( len = * ) INPUT, the name of the input file, or '*'
!    if input comes directly from a user terminal.
!
!    Output, character ( len = * ) LINE, the next chunk of input.
!
!    Input, integer ( kind = 4 ) LUNIT1, is the FORTRAN unit number of the file from
!    which input is read, unless the input file is '*'.
!
!    Input, integer ( kind = 4 ) MRECI, is the maximum number of input records that
!    may be read.
!
!    Input/output integer NRECI, is the current number of input records
!    read.
!
  implicit none

  integer ( kind = 4 ), parameter :: mmax = 500

  character ( len = 6 )  chunk
  character ( len = 20 ) error
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ), save :: indx = 0
  character ( len = 80 ) input
  integer ( kind = 4 ) iread
  integer ( kind = 4 ) lenmy
  logical s_eqi
  character ( len = * ) line
  integer ( kind = 4 ) lunit1
  integer ( kind = 4 ) mreci
  character ( len = mmax ), save :: myline = ' '
  integer ( kind = 4 ) nchunk
  integer ( kind = 4 ) nreci
  integer ( kind = 4 ) nuread

  ierror = 0
!
!  If CHUNK = 'RETURN', then the user is returning some input that
!  was not needed at this time, but which should be stored for
!  later requests.
!
  if ( s_eqi ( chunk, 'RETURN' ) ) then
!
!  Read a word or FORTRAN name at a time.
!
  else if ( s_eqi ( chunk, 'FNAME' ) .or. s_eqi ( chunk, 'WORD' ) ) then

    if ( indx == 0 ) then

      if ( nreci >= mreci ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'GETINP - Warning!'
        write ( *, '(a)' ) '  Maximum number of input records reached.'
        write ( *, '(a)' ) '  (Use MRECI =  to increase this limit.)'
        ierror = 1
        return
      end if

      call getlin ( 'READ', error, input, myline, lunit1, nreci )

      if ( error /= ' ' ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'GETINP - Warning!'
        write ( *, '(a)' ) '  GETLIN reported an error of type'
        write ( *, '(3x,a)' ) trim ( error )
        write ( *, '(a)' ) '  while trying to read a line.'
        ierror = 1
        return
      end if

      indx = 1

    end if

10      continue

    if ( s_eqi ( chunk, 'WORD' ) ) then
      call word_index ( myline, indx, ilo, ihi )
    else if ( s_eqi ( chunk, 'FNAME' ) ) then
      call token_ndx ( myline, indx, ilo, ihi )
    end if

    if ( ilo == 0 ) then

      if ( nreci >= mreci ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'GETINP - Warning!'
        write ( *, '(a)' ) '  Maximum number of input records reached.'
        write ( *, '(a)' ) '  (Use MRECI =  to increase this limit.)'
        ierror = 1
        return
      end if

      call getlin ( 'READ', error, input, myline, lunit1, nreci )

      if ( error /= ' ' ) then
        ierror = 1
        return
      end if

      indx = 1
      go to 10

    else

      indx = indx + 1
      line = myline(ilo:ihi)

    end if
!
!  Read one line of input.
!
  else if ( s_eqi ( chunk, 'LINE' ) ) then

    if ( nreci >= mreci ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'GETINP - Warning!'
      write ( *, '(a)' ) '  Maximum number of input records reached.'
      write ( *, '(a)' ) '  (Use MRECI =  to increase this limit.)'
      ierror = 1
      return
    end if

    call getlin ( 'READ', error, input, line, lunit1, nreci )

    if ( error /= ' ' ) then
      ierror = 1
      return
    end if
!
!  "Packed line": read NREAD characters from input.
!
  else if ( s_eqi ( chunk, 'LINEP' ) ) then

    if ( nreci >= mreci ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'GETINP - Warning!'
      write ( *, '(a)' ) '  Maximum number of input records reached.'
      write ( *, '(a)' ) '  (Use MRECI =  to increase this limit.)'
      ierror = 1
      return
    end if

    iread = 0
    line = ' '

20  continue

    call getlin ( 'READ', error, input, myline, lunit1, nreci )
!
!  Suppress an error return if we read some information successfully.
!
    if ( error /= ' ' ) then
      if ( iread > 0 ) then
        ierror = 0
      else
        ierror = 1
      end if
      return
    end if

    lenmy = len_trim ( myline )

    if ( lenmy > 0 ) then

      nuread = min ( lenmy, nchunk - iread )

      line(iread+1:iread+nuread) = myline(1:nuread)
      iread = iread + nuread
!
!  Put remainder of line, if any, back into input.
!
      if ( nuread < lenmy ) then
        myline = myline(nuread+1: )
        call getlin ( 'UNREAD', error, input, myline, lunit1, nreci )
      end if

    end if

    if ( iread < nchunk ) then
      go to 20
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GETINP - Fatal error!'
    write ( *, '(a)' ) '  Illegal value of CHUNK = ' // trim ( chunk )
    write ( *, '(a)' ) '  Legal values: '
    write ( *, '(a)' ) '    FNAME, LINE, LINEP, RETURN, WORD.'
    stop

  end if

  return
end
subroutine getlin ( action, error, input, line, lunit, nrec )

!*****************************************************************************80
!
!! GETLIN reads (or "unreads") lines from a file.
!
!  Discussion:
!
!    It can also be asked to flush its buffer.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION, tells what GETLIN is being asked
!    to do:
!    'FLUSH', means discard any current lines in the buffer;
!    'READ', means return a line, reading from the file if necessary;
!    'UNREAD', means a line is being returned to the buffer.
!
!    Output, character ( len = * ) ERROR, tells if an error occurred:
!    ' ', no error.
!    'End-of-File' means that no more input could be read.
!    'Read-error' means an "ERR=" error occurred during a read.
!
!    Input, character ( len = * ) INPUT, is the name of the input file, or
!    '*' if input is being read from the user terminal.
!
!    Input/output, character ( len = * ) LINE.
!
!    If ACTION is 'UNREAD', then LINE is an input quantity, and is the
!    line to be returned to the buffer.
!
!    If ACTION is 'READ', then LINE is an output quantity, and is the
!    line of input read from the file.
!
!    If ACTION is 'FLUSH', then LINE is unused.
!
!    Input, integer ( kind = 4 ) LUNIT, the logical unit from which input is to be
!    read.
!
!    Input/output, integer ( kind = 4 ) NREC, is the actual number of successful
!    FORTRAN READ statements carried out.  The user is responsible for
!    initializing this quantity.  GETLIN simply increments it when a READ
!    is carried out.
!
  implicit none

  integer ( kind = 4 ), parameter :: MAXLIN = 10
  integer ( kind = 4 ), parameter :: MMAX = 1000

  character ( len = * ) action
  character ( len = * ) error
  integer ( kind = 4 ) i
  character ( len = * ) input
  integer ( kind = 4 ) ios
  character ( len = * ) line
  character ( len = MMAX ) lines(MAXLIN)
  integer ( kind = 4 ) lunit
  integer ( kind = 4 ), save :: nline = 0
  integer ( kind = 4 ) nrec
  character ( len = MMAX ) nuline
  logical s_eqi

  error = ' '

  if ( s_eqi ( action, 'FLUSH' ) ) then

    nline = 0

  else if ( s_eqi ( action, 'READ'  ) ) then

    if ( nline <= 0 ) then

      if ( input == '*' ) then

        read ( *, '(a)', iostat = ios ) nuline

        if ( ios /= 0 ) then
          error = 'End-of-File'
          return
        end if

        if ( nuline(1:1) == '.' ) then
          error = 'End-of-File'
          return
        end if

      else

        read ( lunit, '(a)', iostat = ios ) nuline

        if ( ios /= 0 ) then
          error = 'End-of-File'
          return
        end if

      end if

      nrec = nrec + 1

      nline = 1
      lines(nline) = nuline

    end if

    line = lines(nline)
    nline = nline-1

  else if ( s_eqi ( action, 'UNREAD' ) ) then

    if ( nline >= MAXLIN ) then
      do i = 1, MAXLIN-1
        lines(i) = lines(i+nline+1-MAXLIN)
      end do
      nline = MAXLIN - 1
    end if

    nline = nline + 1
    lines(nline) = line

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GETLIN - Error!'
    write ( *, '(a)' ) '  The input argument ACTION is unrecognized:'
    write ( *, '(2x,a)' ) trim ( action )
    stop

  end if

  return
end
subroutine hello

!*****************************************************************************80
!
!! HELLO prints out the program's name and last date of revision.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 October 2002
!
!  Author:
!
!    John Burkardt
!
  implicit none


  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  CUT command chops off lines after a string.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  # lines accepted as C/C++/F77/F90 comments.'
  write ( *, '(a)' ) '  F77/F90 lines now automatically de-TABbed.'
  write ( *, '(a)' ) '  Added "<" and ">" options.'
  write ( *, '(a)' ) '  Corrected ICAP/ICAPF/ICAPFC confusion.'
  write ( *, '(a)' ) '  Implemented F77/F90 join.'
  write ( *, '(a)' ) '  Added CHUNK = LINEP option.'
  write ( *, '(a)' ) '  HELP output pauses after a screenful.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  F77/F90 comment lines begin with "!"'
  write ( *, '(a)' ) '    and are corrected automatically.'
  write ( *, '(a)' ) '  "D" now accepted as F77/F90 comment character.'
  write ( *, '(a)' ) '  Set LWRAP to MMAX ( = 500).'
  write ( *, '(a)' ) '  Backslash replaced by CHAR(92).'
  write ( *, '(a)' ) '  Started implementing "breaks";'
  write ( *, '(a)' ) '  (just a "wrap" right now");'
  write ( *, '(a)' ) '  Input is now "buffered" can be "unread";'
  write ( *, '(a)' ) '  STRIPIT can join UNIX lines;'
  write ( *, '(a)' ) '  CHRPAD changed, to fix MARGEL/MARGER problems.'

  return
end
subroutine help

!*****************************************************************************80
!
!! HELP prints out the command options.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 August 1999
!
!  Author:
!
!    John Burkardt
!
  implicit none

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HELP - Note:'
  write ( *, '(a)' ) '  Legal commands to STRIPPER:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Back       Remove character+BACKSPACE sequences;'
  write ( *, '(a)' ) 'Break      Break long lines;'
  write ( *, '(a)' ) '             (ADA/C/C++/F77/F90/TEXT/UNIX not tested.)'
  write ( *, '(a)' ) 'Chunk =    FNAME, LINE, LINEP, or WORD;'
  write ( *, '(a)' ) 'COMBEG =   Set comment begin characters;'
  write ( *, '(a)' ) 'COMEND =   Set comment end characters;'
  write ( *, '(a)' ) 'CommentOut Comment out the text;'
  write ( *, '(a)' ) 'CUT =      Cut off lines after given string.'
  write ( *, '(a)' ) 'Defaults   Restore defaults;'
  write ( *, '(a)' ) 'Go	 Carry out the stripping;'
  write ( *, '(a)' ) 'Help	 Print this information;'
  write ( *, '(a)' ) 'IBlank =   Remove 0 no, 1 all, or 2 multi blank lines;'
  write ( *, '(a)' ) 'ICap =	 -1 lower, +1 upper, +2 first upper text;'
  write ( *, '(a)' ) 'ICapF =    -1 lower, +1 upper, +2 first upper language text;'

  read ( *, * )

  write ( *, '(a)' ) 'ICapFC =   -1 lower, +1 upper, +2 first upper language comments;'
  write ( *, '(a)' ) 'IColumn =  Read column ICOLUMN of a table;'

  write ( *, '(a)' ) 'IComment = 0 do nothing, 1 = delete comments, '
  write ( *, '(a)' ) '           2 = delete noncomments;'
  write ( *, '(a)' ) 'ICon =     -1: Symbol-->control, 0: preserve '
  write ( *, '(a)' ) '           1: Control-->symbol, 2: blanks 3: nothing;'
  write ( *, '(a)' ) 'Input =	 Specify input file, * = screen;'
  write ( *, '(a)' ) '< filename Specify input file, * for screen;'
  write ( *, '(a)' ) 'IRepH =    Horizontal line repeats;'
  write ( *, '(a)' ) 'IRepV =    Vertical line repeats;'
  write ( *, '(a)' ) 'Join =     Join ADA/C/C++/F77/F90/TEXT/UNIX lines;'
  write ( *, '(a)' ) '  	 (ADA/C/C++/TEXT not implemented).'
  write ( *, '(a)' ) 'Keep =     Define keep string;'
  write ( *, '(a)' ) '  	 (Use "<" or ">" to force the string to'
  write ( *, '(a)' ) '           start in column 1, or end in the last.)'
  write ( *, '(a)' ) 'Keep_Max = Maximum length of keeper lines.'
  write ( *, '(a)' ) 'Keep_Min = Minimum length of keeper lines.'
  write ( *, '(a)' ) 'Kill =     Define kill string;'
  write ( *, '(a)' ) '  	 (Use "<" or ">" to force the string to'
  write ( *, '(a)' ) '           start in column 1, or end in the last.)'
  write ( *, '(a)' ) 'Lang =     Set language (ADA/C/C++/F77/F90/TEXT/UNIX);'
  write ( *, '(a)' ) 'Left       Output file will be left justified;'

  read ( *, * )

  write ( *, '(a)' ) 'LWrap =    Set line wrapping margin;'
  write ( *, '(a)' ) 'MargeL =   Set left margin of input;'
  write ( *, '(a)' ) 'MargeR =   Set right margin of output;'
  write ( *, '(a)' ) 'MRecI =    Maximum number of input records;'
  write ( *, '(a)' ) 'MRecO =    Maximum number of output records;'
  write ( *, '(a)' ) 'Number =   -1 delete, 0 preserve, +1 add line numbers;'
  write ( *, '(a)' ) 'Output =   Specify the output file, * = screen;'
  write ( *, '(a)' ) '> filename Specify output file, * for screen;'
  write ( *, '(a)' ) 'Page	 Remove form feeds (new page controls);'
  write ( *, '(a)' ) 'Pause      Pause when printing output to screen;'
  write ( *, '(a)' ) 'Quit       Stop the program;'
  write ( *, '(a)' ) 'ROT13      Encode/decode the output.'
  write ( *, '(a)' ) 'ShoCon     Show control characters;'
  write ( *, '(a)' ) 'Show       Show current settings.'

  return
end
subroutine init ( back, chunk, combeg, comend, commentout, cut, iblank, &
  icap, icapf, icapfc, icolumn, icomment, icon, input, ireph, irepv, join, &
  keep, keep_max, keep_min, kill, lang, left, lwrap, margel, marger, mmax, &
  mreci, mreco, nchunk, number, output, page, pause, dorot13, shocon )

!*****************************************************************************80
!
!! INIT sets or restores the default values of the user options.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, logical BACK, is TRUE if character+BACKSPACE combinations
!    are to be removed.
!
  implicit none

  logical back
  character ( len = 6 )  chunk
  character ( len = 2 ) combeg
  character ( len = 2 ) comend
  logical commentout
  character ( len = 40 ) cut
  integer ( kind = 4 ) iblank
  integer ( kind = 4 ) icap
  integer ( kind = 4 ) icapf
  integer ( kind = 4 ) icapfc
  integer ( kind = 4 ) icolumn
  integer ( kind = 4 ) icomment
  integer ( kind = 4 ) icon
  character ( len = 80 ) input
  integer ( kind = 4 ) ireph
  integer ( kind = 4 ) irepv
  integer ( kind = 4 ) join
  character ( len = 40 ) keep
  integer ( kind = 4 ) keep_max
  integer ( kind = 4 ) keep_min
  character ( len = 40 ) kill
  character ( len = 10 ) lang
  logical left
  integer ( kind = 4 ) lwrap
  integer ( kind = 4 ) margel
  integer ( kind = 4 ) marger
  integer ( kind = 4 ) mmax
  integer ( kind = 4 ) mreci
  integer ( kind = 4 ) mreco
  integer ( kind = 4 ) nchunk
  integer ( kind = 4 ) number
  character ( len = 80 ) output
  logical page
  logical pause
  logical dorot13
  logical shocon

  back = .false.
  chunk = 'LINE'
  combeg = '!'
  comend = ' '
  commentout = .false.
  cut = ' '
  iblank = 0
  icap = 0
  icapf = 0
  icapfc = 0
  icolumn = 0
  icomment = 0
  icon = 0
  input = ' '
  ireph = 1
  irepv = 1
  join = 0
  keep = ' '
  keep_max = mmax
  keep_min = 0
  kill = ' '
  lang = 'TEXT'
  left = .FALSE.
  lwrap = mmax
  margel = 1
  marger = mmax
  mreci = 1000000
  mreco = 1000000
  nchunk = 80
  number = 0
  output = ' '
  page = .false.
  pause = .false.
  dorot13 = .false.
  shocon = .false.

  return
end
function lalpha ( string )

!*****************************************************************************80
!
!! LALPHA returns .TRUE. if STRING contains only alphabetic characters.
!
!  Discussion:
!
!    Alphabetic characters are 'A' through 'Z' and 'a' through 'z' and
!    blanks.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) STRING, the string to be checked.
!
!    Output, logical LALPHA, .TRUE. if STRING contains only
!    alphabetic characters and blanks, .FALSE. otherwise.
!
  implicit none

  character chrtmp
  integer ( kind = 4 ) i
  integer ( kind = 4 ) itemp
  logical lalpha
  integer ( kind = 4 ) lenc
  character ( len = * ) string

  lenc = len ( string )

  lalpha = .false.

  do i = 1, lenc

    chrtmp = string(i:i)
    itemp = ichar ( string(i:i) )

    if ( chrtmp /= ' ' ) then

      if ( .not. ( itemp >= 65 .and. itemp <= 90 ) ) then
        if ( .not. ( itemp >= 97 .and. itemp <= 122 ) ) then
          return
        end if
      end if

    end if

  end do

  lalpha = .true.

  return
end
function ldigit ( string )

!*****************************************************************************80
!
!! LDIGIT returns .TRUE. if STRING contains only digits or blanks.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) STRING, the string to be checked.
!
!    Output, logical LDIGIT, .TRUE. if STRING contains only digits and
!    blanks, .FALSE. otherwise.
!
  implicit none

  character chrtmp
  integer ( kind = 4 ) i
  logical ldigit
  integer ( kind = 4 ) lenc
  character ( len = * ) string

  lenc = len ( string )

  ldigit = .false.

  do i = 1, lenc

    chrtmp = string(i:i)

    if ( chrtmp /= ' ' ) then
      if ( llt ( chrtmp, '0' ) .or. lgt ( chrtmp, '9' ) ) then
        return
      end if
    end if

  end do

  ldigit = .true.

  return
end
subroutine line_get ( chunk, ierror, input, join, lang, line, lunit1, lwrap, &
  mreci, nchunk, nreci )

!*****************************************************************************80
!
!! LINE_GET gets the next input line.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ), parameter :: mmax = 500

  character ( len = 6 )  chunk
  character ( len = 20 ) error
  integer ( kind = 4 ) ierror
  character ( len = 80 ) input
  integer ( kind = 4 ) join
  character ( len = 10 ) lang
  integer ( kind = 4 ) lenl
  character ( len = * ) line
  character ( len = mmax ) line2
  integer ( kind = 4 ) lunit1
  integer ( kind = 4 ) lwrap
  integer ( kind = 4 ) mreci
  integer ( kind = 4 ) nchunk
  integer ( kind = 4 ) nreci
  character s_c_last
  logical s_eqi

  ierror = 0
!
!  Get the next input line (or word).
!
  call getinp ( chunk, ierror, input, line, lunit1, mreci, nchunk, nreci )

  if ( ierror /= 0 ) then
    return
  end if
!
!  If the JOIN option is enabled (join = +1), then ...
!  (We have to say "CHAR(92)" rather than "\" because idiotic UNIX
!  compilers will otherwise take "\" literally.
!
  if ( join == 1 ) then

    if ( s_eqi ( lang, 'ADA' ) ) then

    else if ( s_eqi ( lang, 'C' ) ) then

    else if ( s_eqi ( lang, 'C++' ) ) then
!
!  How do we "look ahead" here?
!
    else if ( s_eqi ( lang, 'F77' ) ) then

      if ( line(1:1) /= '!' ) then

91      continue

        if ( s_c_last ( line ) == '&' ) then

          lenl = len_trim ( line )
          line(lenl:lenl) = ' '

          call getinp ( chunk, ierror, input, line2, lunit1, mreci, nchunk, &
            nreci )

          if ( ierror == 0 ) then
            call s_cat ( line, line2, line )
            go to 91
          end if

        end if

      end if

    else if ( s_eqi ( lang, 'F90' ) ) then

      if ( line(1:1) /= '!' ) then

        do while ( s_c_last ( line ) == '&' )

          lenl = len_trim ( line )
          line(lenl:lenl) = ' '

          call getinp ( chunk, ierror, input, line2, lunit1, mreci, nchunk, &
            nreci )

          if ( ierror /= 0 ) then
            exit
          end if

          call s_cat ( line, line2, line )

        end do

      end if
!
!  TEXT logic.
!
!    Have a current buffer < max.
!    Have a current new line.
!
!    If length of current buffer + new line <= max,
!      add line to buffer,
!      get another line.
!
!    If length of current buffer + new line up to some word <= max,
!      add partial line to buffer, print partial line, come again.
!
!    If length of current buffer + new line word(1) > max, then
!       if ( current buffer > 0 )
!          print buffer and come again
!       else
!         add line up to max -1, add a continuing dash,
!         print buffer, come again.
!
    else if ( s_eqi ( lang, 'TEXT' ) ) then

    else if ( s_eqi ( lang, 'UNIX' ) ) then

      if ( line(1:1) /= '#' ) then

21      continue

        if ( s_c_last ( line ) == char(92) ) then

          lenl = len_trim ( line )
          line(lenl:lenl) = ' '

          call getinp ( chunk, ierror, input, line2, lunit1, mreci, &
            nchunk, nreci )

          if ( ierror == 0 ) then
            call s_cat ( line, line2, line )
            go to 21
          end if

        end if

      end if

    end if
!
!  ...else if the BREAK operation is enabled (JOIN = -1) then
!
  else if ( join == -1 ) then

    if ( s_eqi ( chunk, 'LINE' ) ) then

      if ( lenl > lwrap ) then

        if ( s_eqi ( lang, 'ADA' ) ) then

          call getlin ( 'UNREAD', error, input, line(lwrap+1:lenl), &
            lunit1, nreci )

          if ( error /= ' ') then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'LINE_GET - Warning!'
            write ( *, '(a)' ) '  GETLIN reported an error of type:'
            write ( *, '(2x,a)' ) error
            write ( *, '(a)' ) '  trying to "UNREAD" a long line.'
          end if

          line = line(1:lwrap)

        else if ( s_eqi ( lang, 'C' ) ) then

          call getlin ( 'UNREAD', error, input, line(lwrap+1:lenl), &
            lunit1, nreci )

          if ( error /= ' ') then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'LINE_GET - Warning!'
            write ( *, '(a)' ) '  GETLIN reported an error of type:'
            write ( *, '(3x,a)' ) error
            write ( *, '(a)' ) '  trying to "UNREAD" a long line.'
          end if

          line = line(1:lwrap)

        else if ( s_eqi ( lang, 'C++') ) then

        else if ( s_eqi ( lang, 'F77' ) ) then

          call getlin ( 'UNREAD', error, input, line(lwrap+1:lenl), &
            lunit1, nreci )

          if ( error /= ' ' ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'LINE_GET - Warning!'
            write ( *, '(a)' ) '  GETLIN reported an error of type:'
            write ( *, '(2x,a)' ) error
            write ( *, '(a)' ) '  trying to "UNREAD" a long line.'
          end if

          line = line(1:lwrap)

        else if ( s_eqi ( lang, 'F90' ) ) then

          call getlin ( 'UNREAD', error, input, line(lwrap+1:lenl), &
            lunit1, nreci )

          if ( error /= ' ' ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'LINE_GET - Warning!'
            write ( *, '(a)' ) '  GETLIN reported an error of type:'
            write ( *, '(3x,a)' ) error
            write ( *, '(a)' ) '  trying to "UNREAD" a long line.'
          end if

          line = line(1:lwrap)

        else if ( s_eqi ( lang, 'TEXT' ) ) then

          call getlin ( 'UNREAD', error, input, line(lwrap+1:lenl), &
            lunit1, nreci )

          if ( error /= ' ' ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'LINE_GET - Warning!'
            write ( *, '(a)' ) '  GETLIN reported an error of type:'
            write ( *, '(3x,a)' ) error
            write ( *, '(a)' ) '  trying to "UNREAD" a long line.'
          end if

          line = line(1:lwrap)

        else if ( s_eqi ( lang, 'UNIX' ) ) then

          call getlin ( 'UNREAD', error, input, line(lwrap+1:lenl), &
            lunit1, nreci )

          if ( error /= ' ' ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'LINE_GET - Warning!'
            write ( *, '(a)' ) '  GETLIN reported an error of type:'
            write ( *, '(3x,a)' ) error
            write ( *, '(a)' ) '  trying to "UNREAD" a long line.'
          end if

          line = line(1:lwrap)

        end if

      end if

    end if

  end if

  return
end
subroutine line_put ( ierror, ireph, irepv, isay, line, lunit2, lwrap, mreco, &
  nreco, output, pause )

!*****************************************************************************80
!
!! LINE_PUT writes the line to the output file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error occurred.
!    nonzero, an error occurred.
!
!    Input, integer ( kind = 4 ) IREPH, the number of horizontal repetitions of each line.
!
!    Input, integer ( kind = 4 ) IREPV, the number of vertical repetitions of each line.
!
!    ???, integer ISAY, ???
!
!    Input, character ( len = * ) LINE, ???
!
!    Input, integer ( kind = 4 ) LUNIT2, the output unit.
!
!    Input, integer ( kind = 4 ) LWRAP, the "wrapping length".
!
!    Input/output, integer ( kind = 4 ) MRECO, ???
!
!    Input/output, integer ( kind = 4 ) NRECO, ???
!
!    Input, character ( len = * ) OUTPUT, the output file name, or '*' if output
!    is directly to the screen.
!
!    Input, logical PAUSE, is .TRUE. if output to the screen is to
!    pause.
!
  implicit none

  integer ( kind = 4 ), parameter :: mmax = 500

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) ireph
  integer ( kind = 4 ) irepv
  character isay
  integer ( kind = 4 ) jrepeat
  integer ( kind = 4 ) lenl
  character ( len = mmax ) line
  integer ( kind = 4 ) lunit2
  integer ( kind = 4 ) lwrap
  integer ( kind = 4 ) mreco
  integer ( kind = 4 ) nreco
  character ( len = 80 ) output
  logical pause

  ierror = 0

  lenl = len_trim ( line )
!
!  IREPH: Repeat the line, horizontally.
!  Logically, this operation comes AFTER most others!
!
  if ( ireph > 1 ) then
    ihi = min ( ireph * lenl, mmax )
    do i = lenl+1, ihi
      line(i:i) = line(i-lenl:i-lenl)
    end do
    lenl = ihi
  end if
!
!  IREPV: Repeat the line vertically.
!  Logically, this operation comes AFTER most others!
!
  if ( output == '*' ) then

    if ( isay == '?' .and. pause ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'When the output pauses, press RETURN'
      write ( *, '(a)' ) 'to continue, or any other character to quit.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'Press RETURN now.'
      write ( *, '(a)' ) ' '
      read ( *, '(a)', iostat = ios ) isay

      if ( ios /= 0 ) then
        ierror = 1
        return
      end if

      if ( isay /= ' ' ) then
        ierror = 1
        return
      end if

    end if

    do jrepeat = 1, irepv

      ihi = max ( lenl, 1 )

      do ilo = 1, max ( lenl, 1 ), lwrap

        ihi = min ( ilo+lwrap-1, lenl )
        if ( ihi >= ilo ) then
          write ( *, '(a)' ) line(ilo:ihi)
        else
          write ( *, '(a)' ) ' '
        end if

        nreco = nreco + 1

        if ( nreco >= mreco ) then
          ierror = 1
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'LINE_PUT - Warning!'
          write ( *, '(a)' ) '  Maximum number of output lines reached.'
          write ( *, '(a)' ) '  (Use MRECO =  to increase this limit.)'
          return
        end if

      end do

    end do

    if ( pause .and. mod ( nreco, 20 ) == 19 ) then

      read ( *, '(a)', iostat = ios ) isay

      if ( ios /= 0 ) then
        ierror = 1
        return
      end if

      if ( isay .ne. ' ' ) then
        ierror = 1
        return
      end if

    end if

  else

    do jrepeat = 1, irepv

      do ilo = 1, max ( lenl, 1 ), lwrap

        ihi = min ( lenl, ilo+lwrap-1 )
        if ( ihi >= ilo ) then
          write ( lunit2, '(a)' ) line(ilo:ihi)
        else
          write ( lunit2, '(a)' )
        end if

        nreco = nreco + 1

        if ( nreco >= mreco ) then
          ierror = 1
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'LINE_PUT: Warning!'
          write ( *, '(a)' ) '  Maximum number of output lines reached.'
          return
        end if

      end do

    end do

  end if

  return
end
function numcon ( line )

!*****************************************************************************80
!
!! NUMCON counts the number of control characters in a line of text.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 August 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) LINE, the string of characters to be examined.
!
!    Output, integer ( kind = 4 ) NUMCON is the number of control characters in LINE.
!
  implicit none

  character c
  logical ch_is_control
  integer ( kind = 4 ) i
  integer ( kind = 4 ) lchar
  character ( len = * ) line
  integer ( kind = 4 ) numcon

  numcon = 0

  lchar = len_trim ( line )

  do i = 1, lchar
    c = line(i:i)
    if ( ch_is_control ( c ) ) then
      numcon = numcon + 1
    end if
  end do

  return
end
subroutine rubout ( s )

!*****************************************************************************80
!
!! RUBOUT deletes the pair "character" + BACKSPACE from a string.
!
!  Discussion:
!
!    RUBOUT will also remove a backspace if it is the first character
!    on the line.  RUBOUT is recursive.  In other words, given the
!    string of 8 characters:
!      'ABCD###E'
!    where we are using "#" to represent a backspace, RUBOUT will
!    return the string 'AE'.
!
!    RUBOUT was written for use in "cleaning up" UNICOS MAN pages.
!    The raw text of these MAN pages is unreadable for two reasons:
!
!      Passages which are to be underlined are written so:
!      "_#T_#e_#x_#t" when what is meant is that "Text" is to be
!      underlined if possible.  Note that the seemingly equivalent
!      "T#_e#_x#_t#_" is NOT used.  This is because, in the olden
!      days, certain screen terminals could backspace, but would only
!      display the new character, obliterating rather than
!      overwriting the old one.  This convention allows us to know
!      that we want to delete "character" + Backspace, rather than
!      Backspace + "character".
!
!      Passages which are meant to be in BOLDFACE are written so:
!      "U#U#U#Ug#g#g#gl#l#l#ly#y#y#y", when what is meant is that
!      "Ugly" is to be printed as boldly as possible.  These boldface
!      passages may also be cleaned up using the same rule of
!      removing all occurrences of "character" + Backspace.
!
!    It is truly a fright to look at the text of one of these MAN
!    pages with all the ugly Backspace's, which display on VMS as ^H.
!    These files print or type properly, but look awful in an editor.
!    Moreoever, the lavish use of boldface means that text that is
!    meant to fit in 80 columns can sometimes require 7 times as much
!    space to describe.  This can cause a VMS editor to abort, or to
!    skip the line, since 255 characters is the maximum for EDT.
!
!    A FORTRAN program that tries to read a long line like that will
!    also fail if not careful, since a formatted sequential file
!    on VMS has a default maximum record length of something like
!    133 characters.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) S.  On input, the line of
!    text to be cleaned.  On output, any leading backspace
!    character has been deleted, and all pairs of
!    "character"+Backspace have been deleted.
!
  implicit none

  character, parameter :: BS = char ( 8 )
  integer ( kind = 4 ) i
  integer ( kind = 4 ) nchar
  character ( len = * ) s

  nchar = len_trim ( s )
  i = 1

  do while ( i <= nchar )

    if ( s(i:i) == BS ) then

      if ( i == 1 ) then
        call s_chop ( s, i, i )
        nchar = nchar - 1
        i = i - 1
      else
        call s_chop ( s, i-1, i )
        nchar = nchar - 2
        i = i - 2
      end if

    end if

    i = i + 1

  end do

  return
end
subroutine s_blank_delete ( s )

!*****************************************************************************80
!
!! S_BLANK_DELETE removes blanks from a string, left justifying the remainder.
!
!  Discussion:
!
!    All TAB characters are also removed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) S, the string to be transformed.
!
  implicit none

  character c
  integer ( kind = 4 ) iget
  integer ( kind = 4 ) iput
  character ( len = * ) s
  character ( len = 1 ), parameter :: TAB = char ( 9 )

  iput = 0

  do iget = 1, len ( s )

    c = s(iget:iget)

    if ( c /= ' ' .and. c /= TAB ) then
      iput = iput + 1
      s(iput:iput) = c
    end if

  end do

  s(iput+1:) = ' '

  return
end
function s_c_last ( s )

!*****************************************************************************80
!
!! S_C_LAST returns the last nonblank character in a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string to be examined.
!
!    Output, character S_C_LAST, the last nonblank character in S,
!    or ' ' if S is all blank.
!
  implicit none

  integer ( kind = 4 ) lenc
  character ( len = * ) s
  character s_c_last

  lenc = len_trim ( s )

  if ( lenc > 0 ) then
    s_c_last = s(lenc:lenc)
  else
    s_c_last = ' '
  end if

  return
end
subroutine s_cap ( s )

!*****************************************************************************80
!
!! S_CAP replaces any lowercase letters by uppercase ones in a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) S, the string to be transformed.
!
  implicit none

  character c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) nchar
  character ( len = * ) s

  nchar = len_trim ( s )

  do i = 1, nchar

    c = s(i:i)
    call ch_cap ( c )
    s(i:i) = c

  end do

  return
end
subroutine s_cat ( s1, s2, s3 )

!*****************************************************************************80
!
!! S_CAT concatenates two strings to make a third string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S1, the "prefix" string.
!
!    Input, character ( len = * ) S2, the "postfix" string.
!
!    Output, character ( len = * ) S3, the string made by
!    concatenating S1 and S2, ignoring any trailing blanks.
!
  implicit none

  character ( len = * ) s1
  character ( len = * ) s2
  character ( len = * ) s3

  s3 = trim ( s1 ) // trim ( s2 )

  return
end
subroutine s_chop ( s, ilo, ihi )

!*****************************************************************************80
!
!! S_CHOP "chops out" a portion of a string, and closes up the hole.
!
!  Example:
!
!    S = 'Fred is not a jerk!'
!
!    call s_chop ( S, 9, 12 )
!
!    S = 'Fred is a jerk!    '
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) S, the string to be transformed.
!
!    Input, integer ( kind = 4 ) ILO, IHI, the locations of the first and last
!    characters to be removed.
!
  implicit none

  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ihi2
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) ilo2
  integer ( kind = 4 ) lens
  character ( len = * ) s

  lens = len ( s )

  ilo2 = max ( ilo, 1 )
  ihi2 = min ( ihi, lens )

  if ( ilo2 > ihi2 ) then
    return
  end if

  s(ilo2:lens+ilo2-ihi2-1) = s(ihi2+1:lens)
  s(lens+ilo2-ihi2:lens) = ' '

  return
end
function s_contains_any_alpha ( s )

!*****************************************************************************80
!
!! S_CONTAINS_ANY_ALPHA is TRUE if the string contains any alphabetic character.
!
!  Example:
!
!    Input         Output
!
!    Riding Hood   TRUE
!    123 + 34      FALSE
!    Seven Eleven  TRUE
!    1.0E+11       TRUE
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, a string to be checked.
!
!    Output, logical S_CONTAINS_ANY_ALPHA is TRUE if any character in string
!    is an alphabetic character.
!
  implicit none

  logical ch_is_alpha
  integer ( kind = 4 ) i
  integer ( kind = 4 ) lens
  logical s_contains_any_alpha
  character ( len = * ) s

  lens = len ( s )

  s_contains_any_alpha = .true.

  do i = 1, lens
    if ( ch_is_alpha ( s(i:i) ) ) then
      return
    end if
  end do

  s_contains_any_alpha = .false.

  return
end
subroutine s_control_blank ( s )

!*****************************************************************************80
!
!! S_CONTROL_BLANK replaces control characters with blanks.
!
!  Discussion:
!
!    A "control character" has ASCII code <= 31 or ASCII code => 127.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) S, the string to be transformed.
!
  implicit none

  logical ch_is_control
  integer ( kind = 4 ) i
  character ( len = * ) s

  do i = 1, len ( s )
    if ( ch_is_control ( s(i:i) ) ) then
      s(i:i) = ' '
    end if
  end do

  return
end
subroutine s_control_delete ( s )

!*****************************************************************************80
!
!! S_CONTROL_DELETE removes all control characters from a string.
!
!  Discussion:
!
!    The string is collapsed to the left, and padded on the right with
!    blanks to replace the removed characters.
!
!    A "control character" has ASCII code <= 31 or ASCII code => 127.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) S, is the string to be transformed.
!
  implicit none

  logical ch_is_control
  integer ( kind = 4 ) iget
  integer ( kind = 4 ) iput
  character ( len = * ) s

  iput = 0

  do iget = 1, len ( s )

    if ( .not. ch_is_control ( s(iget:iget) ) ) then
      iput = iput + 1
      s(iput:iput) = s(iget:iget)
    end if

  end do
!
!  Pad the end of the string with blanks
!
  s(iput+1:) = ' '

  return
end
function s_eqi ( strng1, strng2 )

!*****************************************************************************80
!
!! S_EQI is a case insensitive comparison of two strings for equality.
!
!  Example:
!
!    S_EQI ( 'Anjana', 'ANJANA' ) is .TRUE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) STRNG1, STRNG2, the strings to compare.
!
!    Output, logical S_EQI, the result of the comparison.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) len1
  integer ( kind = 4 ) len2
  integer ( kind = 4 ) lenc
  logical s_eqi
  character s1
  character s2
  character ( len = * ) strng1
  character ( len = * ) strng2

  len1 = len ( strng1 )
  len2 = len ( strng2 )
  lenc = min ( len1, len2 )

  s_eqi = .false.

  do i = 1, lenc

    s1 = strng1(i:i)
    s2 = strng2(i:i)
    call ch_cap ( s1 )
    call ch_cap ( s2 )

    if ( s1 /= s2 ) then
      return
    end if

  end do

  do i = lenc + 1, len1
    if ( strng1(i:i) /= ' ' ) then
      return
    end if
  end do

  do i = lenc + 1, len2
    if ( strng2(i:i) /= ' ' ) then
      return
    end if
  end do

  s_eqi = .true.

  return
end
function s_indexi ( s, sub )

!*****************************************************************************80
!
!! S_INDEXI is a case-insensitive INDEX function.
!
!  Discussion:
!
!    The function returns the location in the string at which the
!    substring SUB is first found, or 0 if the substring does not
!    occur at all.
!
!    The routine is also trailing blank insensitive.  This is very
!    important for those cases where you have stored information in
!    larger variables.  If S is of length 80, and SUB is of
!    length 80, then if S = 'FRED' and SUB = 'RED', a match would
!    not be reported by the standard FORTRAN INDEX, because it treats
!    both variables as being 80 characters long!  This routine assumes that
!    trailing blanks represent garbage!
!
!    Because of the suppression of trailing blanks, this routine cannot be
!    used to find, say, the first occurrence of the two-character
!    string 'A '.  However, this routine treats as a special case the
!    occurrence where S or SUB is entirely blank.  Thus you can
!    use this routine to search for occurrences of double or triple blanks
!    in a string, for example, although INDEX itself would be just as
!    suitable for that problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string to be searched.
!
!    Input, character ( len = * ) SUB, the substring to search for.
!
!    Output, integer ( kind = 4 ) S_INDEXI.  0 if SUB does not occur in
!    the string.  Otherwise S(S_INDEXI:S_INDEXI+LENS-1) = SUB,
!    where LENS is the length of SUB, and is the first place
!    this happens.  However, note that this routine ignores case,
!    unlike the standard FORTRAN INDEX function.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) llen1
  integer ( kind = 4 ) llen2
  character ( len = * ) s
  logical s_eqi
  integer ( kind = 4 ) s_indexi
  character ( len = * ) sub

  s_indexi = 0

  llen1 = len_trim ( s )
  llen2 = len_trim ( sub )
!
!  In case S or SUB is blanks, use LEN.
!
  if ( llen1 == 0 ) then
    llen1 = len ( s )
  end if

  if ( llen2 == 0 ) then
    llen2 = len ( sub )
  end if

  if ( llen2 > llen1 ) then
    return
  end if

  do i = 1, llen1 + 1 - llen2

    if ( s_eqi ( s(i:i+llen2-1), sub ) ) then
      s_indexi = i
      return
    end if

  end do

  return
end
subroutine s_low ( string )

!*****************************************************************************80
!
!! S_LOW replaces all uppercase letters by lowercase ones.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) STRING, the string to be
!    transformed.  On output, the string is all lowercase.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) nchar
  character ( len = * ) string

  nchar = len ( string )

  do i = 1, nchar
    call ch_low ( string(i:i) )
  end do

  return
end
function s_only_alphab ( s )

!*****************************************************************************80
!
!! S_ONLY_ALPHAB checks if a string is only alphabetic and blanks.
!
!  Discussion:
!
!    Acceptable characters are 'A' through 'Z' and 'a' through 'z' and blanks.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string to be checked.
!
!    Output, logical S_ONLY_ALPHAB, .TRUE. if the string contains only
!    alphabetic characters and blanks, .FALSE. otherwise.
!
  implicit none

  character c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) itemp
  character ( len = * ) s
  logical s_only_alphab

  s_only_alphab = .false.

  do i = 1, len ( s )

    c = s(i:i)

    if ( c /= ' ' ) then

      itemp = ichar ( c )

      if ( .not. ( itemp >= 65 .and. itemp <= 90 ) ) then
        if ( .not. ( itemp >= 97 .and. itemp <= 122 ) ) then
          return
        end if
      end if

    end if

  end do

  s_only_alphab = .true.

  return
end
function s_only_digitb ( s )

!*****************************************************************************80
!
!! S_ONLY_DIGITB returns .TRUE. if the string contains only digits or blanks.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string to be checked.
!
!    Output, logical S_ONLY_DIGITB, .TRUE. if the string contains only digits
!    and blanks, .FALSE. otherwise.
!
  implicit none

  character c
  integer ( kind = 4 ) i
  character ( len = * ) s
  logical s_only_digitb

  s_only_digitb = .false.

  do i = 1, len ( s )

    c = s(i:i)

    if ( c /= ' ' ) then
      if ( llt ( c, '0' ) .or. lgt ( c, '9' ) ) then
        return
      end if
    end if

  end do

  s_only_digitb = .true.

  return
end
subroutine s_to_rot13 ( s )

!*****************************************************************************80
!
!! S_TO_ROT13 "rotates" the alphabetical characters in a string by 13 positions.
!
!  Discussion:
!
!    Two applications of the routine will return the original string.
!
!  Example:
!
!    Input:                      Output:
!
!    abcdefghijklmnopqrstuvwxyz  nopqrstuvwxyzabcdefghijklm
!    Cher                        Pure
!    James Thurston Howell       Wnzrf Guhefgba Ubjryy
!    12345                       12345
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!   07 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) S, a string to be "rotated".
!
  implicit none

  character ch_to_rot13
  integer ( kind = 4 ) i
  integer ( kind = 4 ) lens
  character ( len = * ) s

  lens = len_trim ( s )

  do i = 1, lens
    s(i:i) = ch_to_rot13 ( s(i:i) )
  end do

  return
end
subroutine show ( back, chunk, combeg, comend, commentout, cut, iblank, &
  icap, icapf, icapfc, icolumn, icomment, icon, input, ireph, irepv, join, &
  keep, kill, lang, left, lwrap, margel, marger, mmax, mreci, mreco, nchunk, &
  number, output, page, pause, dorot13, shocon, showme )

!*****************************************************************************80
!
!! SHOW displays the values of the user variables.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 October 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, logical BACK, is TRUE if character+BACKSPACE combinations
!    are to be removed.
!
!    Input, character ( len = 40 ) CUT, is nonblank if lines containing
!    this string are to be cut at that string.
!
  implicit none

  logical back
  character ( len = 6 )  chunk
  character ( len = 2 ) combeg
  character ( len = 2 ) comend
  logical commentout
  character ( len = 40 ) cut
  integer ( kind = 4 ) iblank
  integer ( kind = 4 ) icap
  integer ( kind = 4 ) icapf
  integer ( kind = 4 ) icapfc
  integer ( kind = 4 ) icolumn
  integer ( kind = 4 ) icomment
  integer ( kind = 4 ) icon
  character ( len = 80 ) input
  integer ( kind = 4 ) ireph
  integer ( kind = 4 ) irepv
  integer ( kind = 4 ) join
  character ( len = 40 ) keep
  character ( len = 40 ) kill
  character ( len = 10 ) lang
  integer ( kind = 4 ) lchar
  logical left
  logical s_eqi
  integer ( kind = 4 ) lwrap
  integer ( kind = 4 ) margel
  integer ( kind = 4 ) marger
  integer ( kind = 4 ) mmax
  integer ( kind = 4 ) mreci
  integer ( kind = 4 ) mreco
  integer ( kind = 4 ) nchunk
  integer ( kind = 4 ) number
  character ( len = 80 ) output
  logical page
  logical pause
  logical dorot13
  logical shocon
  character ( len = 10 ) showme

  if ( s_eqi ( showme, 'NONE' ) ) then
    return
  end if
!
!  BACK
!
  if ( s_eqi ( showme, 'ALL' ) .or. s_eqi ( showme, 'BACK' ) ) then

    if ( back ) then
      write ( *, '(a)' ) '  All character+BACKSPACE sequences will be deleted.'
    else
      write ( *, '(a)' ) '  Character+BACKSPACE sequences will NOT be deleted.'
    end if

  end if
!
!  BREAK
!
  if ( s_eqi ( showme, 'ALL' ) .or. s_eqi ( showme, 'BREAK' ) ) then

    if ( join == -1 ) then
      write ( *, '(a)' ) '  Long ' // trim ( lang ) // ' lines will be broken.'
    else
      write ( *, '(a)' ) '  Long ' // trim ( lang ) // &
        ' lines will NOT be broken.'
    end if

  end if
!
!  CHUNK / NCHUNK
!
  if ( s_eqi ( showme, 'ALL' ) .or. s_eqi ( showme, 'CHUNK' ) .or. &
       s_eqi ( showme, 'NCHUNK' ) ) then
    write ( *, '(a)' ) '  Text is digested in chunks of ' // trim ( chunk )
    write ( *, '(a,i6,a)' ) '  The chunk size is NCHUNK = ', nchunk, &
      ' characters.'
  end if
!
!  COMBEG
!
  if ( s_eqi ( showme, 'ALL' ) .or. s_eqi ( showme, 'COMBEG' ) ) then
    lchar = len_trim ( combeg )
    if ( lchar > 0 ) then
      write ( *, '(a)' ) '  Comments will begin with "' // combeg(1:lchar) &
        // '".'
    end if
  end if
!
!  COMEND
!
  if ( s_eqi ( showme, 'ALL' ) .or. s_eqi ( showme, 'COMEND' ) ) then
    lchar = len_trim ( comend )
    if ( lchar > 0 ) then
      write ( *, '(a)' ) '  Comments will end with "' // comend(1:lchar) &
        // '".'
    end if
  end if
!
!  COMMENTOUT
!
  if ( s_eqi ( showme, 'ALL' ) .or. s_eqi ( showme, 'COMMENTOUT' ) ) then
    if ( commentout ) then
      write ( *, '(a)' ) '  Text will commented out.'
    else
      write ( *, '(a)' ) '  Text will NOT be commented out.'
    end if
  end if
!
!  CUT
!
  if ( s_eqi ( showme, 'ALL') .or. s_eqi ( showme, 'CUT' ) ) then
    write ( *, '(a)' ) '  The cut string is "' // trim ( cut ) // '".'
  end if
!
!  IBLANK
!
  if ( s_eqi ( showme, 'ALL' ) .or. s_eqi ( showme, 'IBLANK' ) ) then
    if ( iblank == 1 ) then
      write ( *, '(a)' ) '  All blank lines will be removed.'
    else if ( iblank == 2 ) then
      write ( *, '(a)' ) '  Double blank lines will be removed.'
    else
      write ( *, '(a)' ) '  Blank lines will NOT be removed.'
    end if
  end if
!
!  ICAP
!
  if ( s_eqi ( showme, 'ALL') .or. s_eqi ( showme, 'icap') ) then
    if ( icap == -1 ) then
      write ( *, '(a)' ) '  Lowercase text characters.'
    else if ( icap == 0 ) then
      write ( *, '(a)' ) '  No case changes made to text.'
    else if ( icap == 1 ) then
      write ( *, '(a)' ) '  Capitalize text characters.'
    else if ( icap == 2 ) then
      write ( *, '(a)' ) '  Initial capitalize every word.'
    else if ( icap == 3 ) then
      write ( *, '(a)' ) '  Initial capitalize every line.'
    end if
  end if
!
!  ICAPF
!
  if ( s_eqi ( showme, 'ALL') .or. s_eqi ( showme, 'ICAPF') ) then
    if ( icapf == -1 ) then
      write ( *, '(a)' ) '  Lowercase '// trim ( lang ) // ' characters.'
    else if ( icapf == 0 ) then
      write ( *, '(a)' ) '  No case changes will be made for ' // &
        trim ( lang ) // '.'
    else if ( icapf == 1 ) then
      write ( *, '(a)' ) '  Capitalize '// trim ( lang ) //' characters.'
    else if ( icapf == 2 ) then
      write ( *, '(a)' ) '  Initial capitalize '// trim ( lang ) // ' words.'
    else if ( icapf == 3 ) then
      write ( *, '(a)' ) '  Initial capitalize '// trim ( lang ) // ' lines.'
    end if
  end if
!
!  ICAPFC
!
  if ( s_eqi ( showme, 'ALL') .or. s_eqi ( showme, 'icapfc') ) then
    if ( icapfc == -1 ) then
      write ( *, '(a)' ) '  Lowercase '// trim ( lang ) // ' comments.'
    else if ( icapfc == 0 ) then
      write ( *, '(a)' ) '  No case changes will be made for ' &
        // trim ( lang ) // ' comments.'
    else if ( icapfc == 1 ) then
      write ( *, '(a)' ) '  Capitalize '// trim ( lang ) // ' comments.'
    else if ( icapfc == 2 ) then
      write ( *, '(a)' ) &
        '  Initial capitalize '// trim ( lang ) // ' comment words.'
    else if ( icapfc == 3 ) then
      write ( *, '(a)' ) '  Initial capitalize ' // trim ( lang ) // &
        ' comment lines.'
    end if
  end if
!
!  ICOLUMN
!
  if ( s_eqi ( showme, 'ALL') .or. s_eqi ( showme, 'icolumn') ) then
    if ( icolumn == 0 ) then
      write ( *, '(a)' ) '  No special column of the input is chosen.'
    else
      write ( *, '(a,i6,a)' ) '  User picks column ', icolumn, ' from table.'
    end if
  end if
!
!  ICON
!
  if ( s_eqi ( showme, 'ALL') .or. s_eqi ( showme, 'ICON') ) then
    if ( icon == -1 ) then
      write ( *, '(a)' ) '  Symbols become control characters.'
    else if ( icon == 0 ) then
      write ( *, '(a)' ) '  Control characters will be preserved.'
    else if ( icon == 1 ) then
      write ( *, '(a)' ) '  Control characters become symbols.'
    else if ( icon == 2 ) then
      write ( *, '(a)' ) '  Control characters replaced by blanks.'
    else if ( icon == 3 ) then
      write ( *, '(a)' ) '  Control characters will be removed.'
    end if
  end if
!
!  ICOMMENT
!
  if ( s_eqi ( showme, 'ALL') .or. s_eqi ( showme, 'ICOMMENT') ) then
    if ( icomment == 0 ) then
      write ( *, '(a)' ) &
        '  No special delete/save of ' // trim ( lang ) //' comments.'
    else if ( icomment == 1 ) then
      write ( *, '(a)' ) '  All '// trim ( lang ) //' comments will be removed.'
    else
      write ( *, '(a)' ) &
        '  ALL '// trim ( lang ) //' NONcomments will be removed.'
    end if
  end if
!
!  INPUT
!
  if ( s_eqi ( showme, 'ALL') .or. s_eqi ( showme, 'INPUT') ) then
    if ( len_trim ( input ) > 0 ) then
      if ( input == '*' ) then
        write ( *, '(a)' ) '  The input file is the screen.'
        write ( *, '(a)' ) '  Terminate input with a period in column 1.'
      else
        write ( *, '(a)' ) '  The input file is "'// trim ( input ) // '".'
      end if
    end if
  end if
!
!  IREPH
!
  if ( s_eqi ( showme, 'ALL') .or. s_eqi ( showme, 'IREPH') ) then
    if ( ireph > 1 ) then
      write ( *, '(a,i6,a)' ) &
        '  Each line is repeated ', ireph, ' times horizontally.'
    else
      write ( *, '(a)' ) '  Each line is NOT repeated horizontally.'
    end if
  end if
!
!  IREPV
!
  if ( s_eqi ( showme, 'ALL' ) .or. s_eqi ( showme, 'IREPV' ) ) then
    if ( irepv > 1 ) then
      write ( *, '(a,i6,a)' ) &
        '  Each line will is repeated ', irepv, ' times vertically.'
    else
      write ( *, '(a)' ) '  Each line is NOT repeated vertically.'
    end if
  end if
!
!  JOIN
!
  if ( s_eqi ( showme, 'ALL') .or. s_eqi ( showme, 'JOIN' ) ) then

    if ( join == 1 ) then
      write ( *, '(a)' ) '  Broken ' // trim ( lang ) // &
        ' lines will be joined.'
    else
      write ( *, '(a)') &
        '  Broken ' // trim ( lang ) // ' lines will NOT be joined.'
    end if

  end if
!
!  KEEP
!
  if ( s_eqi ( showme, 'ALL') .or. s_eqi ( showme, 'KEEP' ) ) then
    write ( *, '(a)' ) '  The keep string is "' // trim ( keep ) // '".'
  end if
!
!  KILL
!
  if ( s_eqi ( showme, 'ALL') .or. s_eqi ( showme, 'KILL' ) ) then
    write ( *, '(a)' ) '  The kill string is "' // trim ( kill ) // '".'
  end if
!
!  LANG
!
  if ( s_eqi ( showme, 'ALL' ) .or. s_eqi ( showme, 'lang' ) ) then
    write ( *, '(a)' ) '  The text format is '// trim ( lang )
  end if
!
!  LEFT
!
  if ( s_eqi ( showme, 'ALL' ) .or. s_eqi ( showme, 'left' ) ) then
    if ( left ) then
      write ( *, '(a)' ) '  Output will be left justified.'
    else
      write ( *, '(a)' ) '  Output will NOT be left justified.'
    end if
  end if
!
!  LWRAP
!
  if ( s_eqi ( showme, 'ALL' ) .or. s_eqi ( showme, 'lwrap' ) ) then
    write ( *, '(a,i6,a)' ) '  Lines will wrap after ', lwrap, ' characters.'
  end if
!
!  MARGEL
!
  if ( s_eqi ( showme, 'ALL') .or. s_eqi ( showme, 'margel') ) then
    write ( *, '(a,i6)' ) '  First column read will be ', margel
  end if
!
!  MARGER
!
  if ( s_eqi ( showme, 'ALL') .or. s_eqi ( showme, 'marger') ) then
    write ( *, '(a,i6)' ) '  Last column read will be ', marger
  end if
!
!  MMAX
!
  if ( s_eqi ( showme, 'ALL') .or. s_eqi ( showme, 'MMAX') ) then
    write ( *, '(a,i6,a)' ) '  Maximum input line length = ', mmax, &
      ' characters.'
  end if
!
!  MRECI
!
  if ( s_eqi ( showme, 'ALL') .or. s_eqi ( showme, 'MRECI') ) then
    write ( *, '(a,i12)' ) '  Maximum number of input records is ', mreci
  end if
!
!  MRECO
!
  if ( s_eqi ( showme, 'ALL') .or. s_eqi ( showme, 'MRECO') ) then
    write ( *, '(a,i12)' ) '  Maximum number of output records is ', mreco
  end if
!
!  NUMBER command
!
  if ( s_eqi ( showme, 'ALL' ) .or. s_eqi ( showme, 'NUMBER' ) ) then
    if ( number < 0 ) then
      write ( *, '(a)' ) '  Initial line numbers will be stripped.'
    else if ( number == 0 ) then
      write ( *, '(a)' ) '  No special numbering option.'
    else if ( number > 0 ) then
      write ( *, '(a)' ) '  Initial line numbers will be inserted.'
    end if
  end if
!
!  OUTPUT
!
  if ( s_eqi ( showme, 'ALL') .or. s_eqi ( showme, 'output' ) ) then
    if ( lchar > 0 ) then
      if ( output == '*' ) then
        write ( *, '(a)' ) '  The output file is the screen.'
      else
        write ( *, '(a)' ) '  The output file is "' // trim ( output ) // '".'
      end if
    end if
  end if
!
!  PAGE
!
  if ( s_eqi ( showme, 'ALL') .or. s_eqi ( showme, 'PAGE' ) ) then
    if (page) then
      write ( *, '(a)' ) '  New Page (FF) controls will be removed.'
    else
      write ( *, '(a)' ) '  New Page (FF) controls will NOT be removed.'
    end if
  end if
!
!  PAUSE
!
  if ( s_eqi ( showme, 'ALL') .or. s_eqi ( showme, 'PAUSE' ) ) then
    if (pause) then
      write ( *, '(a)' ) '  Output to the screen will pause.'
    else
      write ( *, '(a)' ) '  Output to the screen will NOT pause.'
    end if
  end if
!
!  ROT13
!
  if ( s_eqi ( showme, 'ALL' ) .or. s_eqi ( showme, 'ROT13' ) ) then

    if ( dorot13 ) then
      write ( *, '(a)' ) '  Output will be ROT13 encoded/decoded.'
    else
      write ( *, '(a)' ) '  Output will NOT be ROT13 encoded/decoded.'
    end if

  end if
!
!  SHOCON
!
  if ( s_eqi ( showme, 'ALL' ) .or. s_eqi ( showme, 'SHOCON' ) ) then

    if ( shocon ) then
      write ( *, '(a)' ) '  Control characters will be identified.'
    else
      write ( *, '(a)' ) '  Control characters will NOT be identified.'
    end if

  end if

  return
end
subroutine stripit ( back, chunk, combeg, comend, commentout, cut, iblank, &
  icap, icapf, icapfc, icolumn, icomment, icon, input, ireph, irepv, join, &
  keep, keep_max, keep_min, kill, lang, left, lunit1, lunit2, lwrap, margel, &
  marger, mreci, mreco, nchunk, number, output, page, pause, dorot13, shocon )

!*****************************************************************************80
!
!! STRIPIT processes the file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 October 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, logical BACK, is TRUE if character+BACKSPACE combinations
!    are to be removed.
!
  implicit none

  integer ( kind = 4 ), parameter :: mmax = 500

  logical back
  character c
  logical ch_is_control
  character ( len = 6 )  chunk
  character ( len = 2 ) combeg
  character ( len = 2 ) comend
  logical comment
  logical commentout
  character ( len = 40 ) cut
  character ffchar
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iblank
  integer ( kind = 4 ) icap
  integer ( kind = 4 ) icapf
  integer ( kind = 4 ) icapfc
  integer ( kind = 4 ) icolumn
  integer ( kind = 4 ) icomment
  integer ( kind = 4 ) icon
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) imhere
  character ( len = 80 ) input
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) ireph
  integer ( kind = 4 ) irepv
  integer ( kind = 4 ) isblnk
  character isay
  integer ( kind = 4 ) itwo
  integer ( kind = 4 ) join
  integer ( kind = 4 ) jsblnk
  character ( len = 40 ) keep
  integer ( kind = 4 ) keep_max
  integer ( kind = 4 ) keep_min
  character ( len = 40 ) kill
  character ( len = 10 ) lang
  integer ( kind = 4 ) lchar
  logical ldigit
  logical left
  integer ( kind = 4 ) lenl
  logical s_eqi
  character ( len = mmax ) line
  character ( len = mmax ) line2
  logical lnexcom
  logical l_temp
  integer ( kind = 4 ) lunit1
  integer ( kind = 4 ) lunit2
  integer ( kind = 4 ) lwrap
  integer ( kind = 4 ) margel
  integer ( kind = 4 ) marger
  integer ( kind = 4 ) mleni
  integer ( kind = 4 ) mleno
  integer ( kind = 4 ) mnumi
  integer ( kind = 4 ) mnumo
  integer ( kind = 4 ) mreci
  integer ( kind = 4 ) mreco
  integer ( kind = 4 ) nblank
  integer ( kind = 4 ) nchri
  integer ( kind = 4 ) nchro
  integer ( kind = 4 ) nchunk
  integer ( kind = 4 ) nconi
  integer ( kind = 4 ) ncut
  integer ( kind = 4 ) nff
  integer ( kind = 4 ) nkeep
  integer ( kind = 4 ) nkill
  integer ( kind = 4 ) nreci
  integer ( kind = 4 ) nreco
  integer ( kind = 4 ) number
  integer ( kind = 4 ) numcon
  character ( len = 80 ) output
  logical page
  logical pause
  logical dorot13
  logical shocon
  logical s_contains_any_alpha
  integer ( kind = 4 ) s_indexi

  ierror = 0
  isay = '?'
  isblnk = 1
  itwo = 1
  jsblnk = 1
  line2 = ' '
  mleni = 0
  mleno = 0
  mnumi = 0
  mnumo = 0
  nblank = 0
  nchri = 0
  nchro = 0
  nconi = 0
  ncut = 0
  nff = 0
  nkeep = 0
  nkill = 0
  nreci = 0
  nreco = 0
!
!  The user must have specified an input file name.
!
  if ( len_trim ( input ) <= 0 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'STRIPIT: Error!'
    write ( *, '(a)' ) '  Please specify an input file!'
    return

  end if
!
!  The user must have specified an output file name.
!
  if ( len_trim ( output ) <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'STRIPIT: Not enough information!'
    write ( *, '(a)' ) '  Please specify an output file first!'
    return
  end if
!
!  Open the input file.
!
  if ( input /= '*' ) then

    open ( unit = lunit1, file = input, status = 'old', iostat = ios, &
      form = 'formatted', access = 'sequential' )

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'STRIPIT - Warning!'
      write ( *, '(a)' ) '  Error while trying to open the input file:'
      write ( *, '(a)' ) trim ( input )
      return
    end if

  end if
!
!  Open the output file.
!
  if ( output /= '*' ) then

    open ( unit = lunit2, file = output, status = 'new', iostat = ios, &
      form = 'formatted', access = 'sequential' )

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'STRIPIT - Warning!'
      write ( *, '(a)' ) '  Error while trying to open the output file:'
      write ( *, '(a)' ) trim ( output )
      return
    end if

  end if
!
!*******************************************************************************
!
!  READ A LOGICAL LINE.
!
!*******************************************************************************
!
!  A logical line may be the result of joining several physical lines,
!  or dividing a physical input line into smaller pieces.
!
  do

    call line_get ( chunk, ierror, input, join, lang, line, lunit1, &
      lwrap, mreci, nchunk, nreci )

    if ( ierror /= 0 ) then
      exit
    end if
!
!*******************************************************************************
!
!  UPDATE INPUT LINE STATISTICS.
!
!*******************************************************************************
!
    lenl = len_trim ( line )
    nchri = nchri + lenl

    if ( lenl > mleni ) then
      mleni = lenl
      mnumi = nreci
    end if

    nconi = nconi + numcon ( line )
!
!*******************************************************************************
!
!  ADJUST THE INPUT LINE.
!
!*******************************************************************************
!
!  MARGIN: Handle the left and right margins.
!
    call doml ( line, margel, marger )
!
!  If we're only reading one column, extract column ICOLUMN.
!
    if ( icolumn > 0 ) then
      call word_index ( line, icolumn, ilo, ihi )
      line = line(ilo:ihi)
    end if
!
!  DETAB the input if it is F77/F90.
!
    if ( s_eqi ( lang, 'F77' ) .or. s_eqi ( lang, 'F90' ) ) then
      call chrdt6 ( line )
    end if
!
!*******************************************************************************
!
!  PROCESS THE INPUT LINE.
!
!*******************************************************************************
!
!  Determine if this line begins, continues, or ends a
!  series of comment lines.
!
!  Let's pretend there are only two kinds of comments:
!    INLINE comments, which trail a statement, and
!    MULTILINE comments, comprising one or more full line comments.
!
!  In that case, I've handled the MULTILINE case, and can't
!  handle the INLINE case yet.
!
!    LNEXCOM  =
!      TRUE:  the next line MUST be a comment line.
!      FALSE: the next line may be a comment or executable line.
!
    if ( s_eqi ( lang, 'ADA' ) ) then

      lnexcom = .false.

      if ( line(1:2) == '--' ) then
        comment = .true.
      else
        comment = .false.
      end if

    else if ( s_eqi ( lang, 'C' ) ) then

      if ( line(1:2) == '/*' ) then

        comment = .true.
        if ( index(line(3:), '*/' ) == 0 ) then
          lnexcom = .true.
        else
          lnexcom = .false.
        end if

      else if ( line(1:1) == '#' ) then

        comment = .true.

      else if ( lnexcom ) then

        comment = .true.

        if ( index ( line, '*/' ) == 0 ) then
          lnexcom = .true.
        else
          lnexcom = .false.
        end if

      else

        comment = .false.

      end if

    else if ( s_eqi ( lang, 'C++' ) ) then

      lnexcom = .false.

      if ( line(1:2) == '//' ) then
        comment = .true.
      else if ( line(1:1) == '#' ) then
        comment = .true.
      else
        comment = .false.
      end if

    else if ( s_eqi ( lang, 'F77' ) .or. s_eqi ( lang, 'F90' ) ) then

      lnexcom = .false.

      if ( line(1:1) == 'c' .or. line(1:1) == 'C' .or. &
           line(1:1) == 'd' .or. line(1:1) == 'D' .or. &
           line(1:1) == '*' .or. line(1:1) == '#' .or. &
           line(1:1) == '!' ) then

        comment = .true.
        line(1:1) = '!'

      else

        comment = .false.

      end if

    else if ( s_eqi ( lang, 'UNIX' ) ) then

      lnexcom = .false.

      if ( line(1:1) == '#' ) then
        comment = .true.
      else
        comment = .false.
      end if

    else

      comment = .false.
      lnexcom = .false.

    end if
!
!  If KEEP is nonblank, see if LINE contains the token, and either
!  save or discard those lines.  Sneakily pre- and postpend "<" and
!  ">" to allow the specification of first and last columns in a
!  token.
!
    if ( keep /= ' ' ) then

      if ( s_eqi ( keep, 'ANY_ALPHA' ) ) then

        l_temp = s_contains_any_alpha ( line )

        if ( l_temp ) then
          nkeep = nkeep + 1
        else
          go to 60
        end if

      else

        if ( len_trim ( line ) <= 0 ) then

          go to 60

        else

          line2 = '<' // trim ( line ) // '>'
          lchar = len_trim ( keep )
          imhere = s_indexi ( line2, keep(1:lchar) )

          if ( imhere == 0 ) then
            go to 60
          else
            nkeep = nkeep + 1
          end if

        end if

      end if

    end if
!
!  Check that the line length is between KEEP_MIN and KEEP_MAX
!
    lenl = len_trim ( line )
    if ( lenl < keep_min .or. lenl > keep_max ) then
      go to 60
    end if
!
!  If KILL is nonblank, see if LINE contains the token, and either
!  save or discard those lines.  Sneakily pre- and postpend "<" and
!  ">" to allow the specification of first and last columns in a
!  token.
!
    if ( kill /= ' ' ) then

      if ( len_trim ( line ) > 0 ) then

        line2 = '<' // trim ( line ) // '>'
        imhere = s_indexi ( line2, trim ( kill ) )

        if ( imhere /= 0 ) then
          nkill = nkill + 1
          go to 60
        end if

      end if
    end if
!
!  If CUT is nonblank, see if LINE contains the token and cut the line.
!
    if ( cut /= ' ' ) then

      if ( len_trim ( line ) <= 0 ) then
        go to 60
      end if

      lchar = len_trim ( cut )
      imhere = s_indexi ( line, cut(1:lchar) )

      if ( imhere /= 0 ) then

        ncut = ncut + 1

        if ( imhere == 1 ) then
          line = ' '
          lenl = 0
        else
          line = line(1:imhere-1)
          lenl = imhere - 1
        end if

      end if

    end if
!
!  BACK:  character+BACKSPACE deletion
!
    if ( back ) then
      call rubout ( line )
      lenl = len_trim ( line )
    end if
!
!  COMMENTOUT
!
    if ( commentout ) then

      lchar = len_trim ( combeg )

      if ( len_trim ( combeg ) > 0 ) then
        line(lchar+3:) = trim ( line )
        line(1:lchar) = combeg(1:lchar)
        line(lchar+1:lchar+2) = '  '
      end if

      lchar = len_trim ( comend )

      if ( len_trim ( comend ) > 0 ) then
        lenl = len_trim ( line )
        line(lenl+1:lenl+1) = ' '
        line(lenl+2:lenl+lchar+1) = comend(1:lchar)
      end if

    end if
!
!  If LANG = F77, force the continuation character to be '&'.
!
    if ( s_eqi ( lang, 'F77' ) ) then
      if ( .not. comment .and. line(6:6) /= ' ' ) then
        line(6:6) = '&'
      end if
    end if
!
!  IBLANK: Handle blank or double blank lines.
!
    jsblnk = isblnk
    isblnk = 0

    if ( lenl == 0 ) then
      isblnk = 1
      if ( iblank == 1 ) then
        nblank = nblank + 1
        go to 60
      else if ( iblank == 2 .and. jsblnk == 1 ) then
        nblank = nblank + 1
        go to 60
      end if
    end if
!
!  ICAP
!
    if ( icap == -1 ) then
      call s_low ( line )
    else if ( icap == 1 ) then
      call s_cap ( line )
    else if ( icap == 2 ) then
      call word_cap ( line )
    else if ( icap == 3 ) then
      call ch_cap ( line(1:1) )
    end if
!
!  ICAPF
!
    if ( s_eqi ( lang, 'ADA' ) .or. &
         s_eqi ( lang, 'C' ) .or. &
         s_eqi ( lang, 'C++' ) .or. &
         s_eqi ( lang, 'F77' ) .or. &
         s_eqi ( lang, 'F90' ) .or. &
         s_eqi ( lang, 'UNIX' ) ) then

      if ( .not. comment ) then

        if ( icapf == -1 ) then
          call s_low ( line )
        else if ( icapf == 1 ) then
          call s_cap ( line )
        else if ( icapf == 2 ) then
          call word_cap ( line )
        else if ( icapf == 3 ) then
          call ch_cap ( line(1:1) )
        end if

      end if
    end if
!
!  ICAPFC
!
    if ( s_eqi ( lang, 'ADA' ) .or. s_eqi ( lang, 'C' ) .or. &
         s_eqi ( lang, 'C++' ) .or. s_eqi ( lang, 'F77' ) .or. &
         s_eqi ( lang, 'F90' ) .or. s_eqi ( lang, 'UNIX' ) ) then

      if ( comment ) then

        if ( icapfc == -1 ) then
          call s_low ( line )
        else if ( icapf == 1 ) then
          call s_cap ( line )
        else if ( icapf == 2 ) then
          call word_cap ( line )
        else if ( icapf == 3 ) then
          call ch_cap ( line(1:1) )
        end if

      end if
    end if
!
!  ICON:
!    -1 Symbols become controls;
!     0 No action;
!    +1 Controls go to symbols;
!    +2 Controls go to blanks;
!    +3 Controls go to nothing.
!
    if ( icon == -1 ) then
      call chrs_to_a ( line, line2 )
      line = line2
    else if ( icon == 1 ) then
      call chra_to_s ( line, line2 )
      line = line2
    else if ( icon == 2 ) then
      call s_control_blank ( line )
    else if ( icon == 3 ) then
      call s_control_delete ( line )
    end if
!
!  ICOMMENT: Do nothing, or remove comments or remove noncomments.
!
    if ( s_eqi ( lang, 'ADA' ) .or. s_eqi ( lang, 'C' ) .or. &
         s_eqi ( lang, 'C++' ) .or. s_eqi ( lang, 'F77' ) .or. &
         s_eqi ( lang, 'F90' ) .or. s_eqi ( lang, 'UNIX' ) ) then

      if ( icomment == 1 ) then
        if ( comment ) go to 60
      else if ( icomment == 2 ) then
        if ( .not. comment ) go to 60
      end if

    end if
!
!  LEFT: Left justify the text.
!
    if ( left ) then
      line = adjustl ( line )
      lenl = len_trim ( line )
    end if
!
!  NUMBER: Remove or insert initial digits.
!
    if ( number < 0 ) then

      lenl = len_trim ( line )
      do i = 1, lenl
        if ( ldigit ( line(i:i) ) ) then
          line(i:i) = ' '
        else
          exit
        end if
      end do

      line = adjustl ( line )
      lenl = len_trim ( line )

    else if ( number > 0 ) then

      lenl = len_trim ( line )
      line(1+5:lenl+5) = line(1:lenl)
      write ( line(1:4), '(i4)' ) nreci
      line(5:5) = ' '

    end if
!
!  PAGE: remove FF (new page) controls.
!
    if ( page ) then

      ffchar = char ( 12 )
      lenl = len_trim ( line )

      do i = 1, lenl
        if ( line(i:i) == ffchar ) then
          line(i:i) = ' '
          nff = nff + 1
        end if
      end do

    end if
!
!  SHOCON: Identify control characters
!
    if ( shocon ) then
      do i = 1, lenl
        c = line(i:i)
        if ( ch_is_control ( c ) ) then
          write ( *, '(a,i6,a,i6,a,i6)' ) &
            'Line ', nreci, ' Column = ', i, ' ASCII # = ', ichar ( c )
        end if
      end do
    end if
!
!  Now, just before output, apply ROT13, if requested.
!
    if ( dorot13 ) then
      call s_to_rot13 ( line )
    end if
!
!*******************************************************************************
!
!  UPDATE THE OUTPUT LINE STATISTICS.
!
!*******************************************************************************
!
    lenl = len_trim ( line )
    lenl = max ( 1, lenl )
!
!  Compute the output statistics.
!
    nchro = nchro + lenl

    if ( lenl > mleno ) then
      mleno = lenl
      mnumo = nreco + 1
    end if
!
!*******************************************************************************
!
!  WRITE THE OUTPUT LINE.
!
!*******************************************************************************
!
    call line_put ( ierror, ireph, irepv, isay, line, lunit2, &
      lwrap, mreco, nreco, output, pause )

    if ( ierror /= 0 ) then
      exit
    end if

60  continue
!
!  Print running indicator of number of lines processed so far.
!
    if ( nreci == itwo ) then
      if ( output /= '*' ) then
        if ( itwo == 1 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '   Input  Output'
          write ( *, '(a)' ) ' '
        end if
        write ( *, '(2i8)' ) nreci, nreco
      end if
      itwo = 2 * itwo
    end if

  end do
!
!*******************************************************************************
!
!  PRINT THE STATISTICS.
!
!*******************************************************************************
!
  if ( nreci /= 2 * itwo ) then
    if ( output /= '*' ) then
      write ( *, '(2i8)' ) nreci, nreco
    end if
  end if

  if ( input /= '*' ) then
    close ( unit = lunit1 )
  end if

  if ( output /= '*' ) then
    close ( unit = lunit2 )
  end if

  write ( *, '(a)'          ) ' '
  write ( *, '(a)'          ) '                     Input    Output'
  write ( *, '(a)'          ) ' '
  write ( *, '(a,i8,2x,i8)' ) '  Lines:         ', nreci, nreco
  write ( *, '(a,i8,2x,i8)' ) '  Characters:    ', nchri, nchro
  write ( *, '(a,i8      )' ) '  Controls:      ', nconi
  write ( *, '(a,i8,2x,i8)' ) '  Longest line:  ', mleni, mleno
  write ( *, '(a,i8,2x,i8)' ) '  (line number:) ', mnumi, mnumo
  write ( *, '(a)' ) ' '

  if ( nff > 0 ) then
    write ( *, '(a,i6,a)' ) '  Removed ', nff, ' form feeds.'
  end if

  if ( iblank /= 0 ) then
    write ( *, '(a,i6,a)' ) '  Removed ', nblank, ' blank lines.'
  end if

  if ( keep /= ' ' ) then
    write ( *, '(a,i6,a)' ) '  Kept ', nkeep, ' lines with KEEP string.'
  end if

  if ( kill /= ' ' ) then
    write ( *, '(a,i6,a)' ) '  Killed ', nkill, ' lines with KILL string.'
  end if

  if ( cut /= ' ' ) then
    write ( *, '(a,i6,a)' ) '  Cut ', ncut, ' lines with CUT string.'
  end if


  return
end
subroutine sym_to_ch ( sym, c, ihi )

!*****************************************************************************80
!
!! SYM_TO_CH returns the character represented by a symbol.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) SYM is a string containing printable symbols.
!
!    Output, character C, is the ASCII character represented by the
!    first symbol in SYM.
!
!    Output, integer ( kind = 4 ) IHI, C is represented by SYM(1:IHI).
!    IHI = 0 if there was a problem.
!
  implicit none

  character c
  integer ( kind = 4 ) ialt
  integer ( kind = 4 ) ichr
  integer ( kind = 4 ) ictl
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) nchar
  logical s_eqi
  character ( len = * ) sym

  c = ' '
  nchar = len_trim ( sym )

  if ( nchar <= 0 ) then
    c = ' '
    ihi = 0
    return
  end if

  ialt = 0
  ictl = 0
  ihi = 1
!
!  Could it be an ALT character?
!
  if ( sym(ihi:ihi) == '!' .and. ihi < nchar ) then
    ialt = 1
    ihi = ihi + 1
  end if
!
!  Could it be a control character?
!
  if ( sym(ihi:ihi) == '^' .and. ihi < nchar ) then
    ictl = 1
    ihi = ihi + 1
  end if
!
!  Could it be a DEL character?
!
  ichr = ichar ( sym(ihi:ihi) )

  if ( ihi+2 <= nchar ) then
    if ( s_eqi ( sym(ihi:ihi+2), 'DEL' ) ) then
      ichr = 127
      ihi = ihi + 2
    end if
  end if
!
!  Could it be an SP character?
!
  if ( ihi + 1 <= nchar ) then
    if ( s_eqi ( sym(ihi:ihi+1), 'SP' ) ) then
      ichr = 32
      ihi = ihi + 1
    end if
  end if
!
!  Interpret the character.
!
  if ( ialt == 1 ) then
    ichr = ichr + 128
  end if

  if ( ictl == 1 ) then
    ichr = ichr - 64
  end if

  c = char ( ichr )

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
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
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
subroutine token_ndx ( string, indx, ilo, ihi )

!*****************************************************************************80
!
!! TOKEN_NDX finds the N-th FORTRAN variable name in a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) STRING is the string of words to be analyzed.
!
!    Input, integer ( kind = 4 ) INDX is the index of the desired token.
!
!    Output, integer ( kind = 4 ) ILO is the index of the first character of the
!    INDX-th token, or 0 if there was no INDX-th token.
!
!    Output, integer ( kind = 4 ) IHI is the index of the last character of the
!    INDX-th token, or 0 if there was no INDX-th token.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) indx
  character ( len = * ) string

  ihi = 0
  ilo = 0

  do i = 1, indx

    call token_next ( string, ilo, ihi)

    if ( ilo == 0 ) then
      return
    end if

  end do

  return
end
subroutine token_next ( s, ilo, ihi )

!*****************************************************************************80
!
!! TOKEN_NEXT finds the next FORTRAN variable name in a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S is the string of words to be analyzed.
!
!    Output, integer ( kind = 4 ) ILO is the location of the first character of the
!    next word, or 0 if there was no next word.
!
!    Input/output, integer ( kind = 4 ) IHI.
!    On input, IHI is taken to be the LAST character of the
!    PREVIOUS word, or 0 if the first word is sought.
!
!    On output, IHI is the index of the last character of
!    the next word, or 0 if there was no next word.
!
  implicit none

  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) lchar
  character ( len = * ) s
  logical s_only_alphab
  logical s_only_digitb

  lchar = len_trim ( s )

  ilo = ihi

  if ( ilo < 0 ) then
    ilo = 0
  end if
!
!  Find ILO, the index of the next alphabetic character.
!
  do

    ilo = ilo + 1

    if ( ilo > lchar ) then
      ilo = 0
      ihi = 0
      return
    end if

    if ( s_only_alphab ( s(ilo:ilo) ) ) then
      exit
    end if

  end do
!
!  Find the index of the next character which is neither
!  alphabetic nor numeric.
!
  ihi = ilo

  do

    ihi = ihi + 1

    if ( ihi > lchar ) then
      ihi = lchar
      return
    end if

    if ( .not. ( s_only_alphab ( s(ihi:ihi) ) ) .and. &
         .not. ( s_only_digitb ( s(ihi:ihi) ) ) ) then
      exit
    end if

  end do

  ihi = ihi - 1

  return
end
subroutine word_cap ( s )

!*****************************************************************************80
!
!! WORD_CAP capitalizes the first character of each word in a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) S, the string to be transformed.
!
  implicit none

  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  character ( len = * ) s

  ilo = 0
  ihi = 0

  do

    call word_next ( s, ilo, ihi )

    if ( ilo <= 0 ) then
      exit
    end if

    call ch_cap ( s(ilo:ilo) )

  end do

  return
end
subroutine word_index ( s, indx, ilo, ihi )

!*****************************************************************************80
!
!! WORD_INDEX finds the word of a given index in a string.
!
!  Discussion:
!
!    The routine returns in ILO and IHI the beginning and end of the INDX-th
!    word, or 0 and 0 if there is no INDX-th word.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S is the string of words to be analyzed.
!
!    Input, integer ( kind = 4 ) INDX is the index of the desired token.
!
!    Output, integer ( kind = 4 ) ILO is the index of the first character of the
!    INDX-th word, or 0 if there was no INDX-th word.
!
!    Output, integer ( kind = 4 ) IHI is the index of the last character of the INDX-th
!    word, or 0 if there was no INDX-th word.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) indx
  character ( len = * ) s

  ihi = 0
  ilo = 0

  do i = 1, indx

    call word_next ( s, ilo, ihi )

    if ( ilo == 0 ) then
      return
    end if

  end do

  return
end
subroutine word_next ( s, ilo, ihi )

!*****************************************************************************80
!
!! WORD_NEXT finds the next (blank separated) word in a string.
!
!  Discussion:
!
!    This routine is usually used repetitively on a fixed string.  On each
!    call, it accepts IHI, the index of the last character of the
!    previous word extracted from the string.
!
!    It then computes ILO and IHI, the first and last characters of
!    the next word in the string.
!
!    It is assumed that words are separated by one or more spaces.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) S, the string of words to be analyzed.
!
!    Output, integer ( kind = 4 ) ILO is the location of the first character of the
!    next word, or 0 if there was no next word.
!
!    Input/output, integer ( kind = 4 ) IHI.
!
!    On input, IHI is taken to be the LAST character of the
!    PREVIOUS word, or 0 if the first word is sought.
!
!    On output, IHI is the index of the last character of
!    the next word, or 0 if there was no next word.
!
  implicit none

  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) lchar
  character ( len = * ) s

  lchar = len_trim ( s )
!
!  Find ILO, the index of the first nonblank character after
!  (the old value of) IHI.
!
  if ( ihi < 0 ) then
    ilo = 0
  else
    ilo = ihi
  end if

  do

    ilo = ilo + 1

    if ( ilo > lchar ) then
      ilo = 0
      ihi = 0
      return
    end if

    if ( s(ilo:ilo) /= ' ') then
      exit
    end if

  end do
!
!  Find IHI, the index of the next blank character, or end of line.
!
  ihi = ilo

  do

    ihi = ihi + 1

    if ( ihi >= lchar ) then
      ihi = lchar
      return
    end if

    if ( s(ihi:ihi) == ' ' ) then
      exit
    end if

  end do
!
!  Decrement IHI to point to the previous, nonblank, character.
!
  ihi = ihi - 1

  return
end
subroutine word_next2 ( s, first, last )

!*****************************************************************************80
!
!! WORD_NEXT2 returns the first word in a string.
!
!  Discussion:
!
!    "Words" are any string of characters, separated by commas or blanks.
!
!    The routine returns:
!    * FIRST, the first string of nonblank, noncomma characters;
!    * LAST, the characters of STRING that occur after FIRST and
!      the commas and blanks.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) S is the string to be analyzed.
!
!    Output, character ( len = * ) FIRST, the next word in the string.
!
!    Output, character ( len = * ) LAST, the text of the remaining words in the
!    string.
!
  implicit none

  character c
  character ( len = * ) first
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ido
  integer ( kind = 4 ) ifirst
  integer ( kind = 4 ) ilast
  character ( len = * ) last
  integer ( kind = 4 ) lenf
  integer ( kind = 4 ) lenl
  integer ( kind = 4 ) lens
  character ( len = * ) s

  first = ' '
  last = ' '

  ifirst = 0
  ilast = 0

  lens = len_trim ( s )
  lenf = len ( first )
  lenl = len ( last )

  ido = 0

  do i = 1, lens

    c = s(i:i)

    if ( ido == 0 ) then
      if ( c /= ' ' .and. c /= ',' ) then
        ido = 1
      end if
    end if

    if ( ido == 1 ) then
      if ( c /= ' ' .and. c /= ',' ) then
        ifirst = ifirst + 1
        if ( ifirst <= lenf ) then
          first(ifirst:ifirst) = c
        end if
      else
        ido = 2
      end if
    end if

    if ( ido == 2 ) then
      if ( c /= ' ' .and. c /= ',' ) then
        ido = 3
      end if
    end if

    if ( ido == 3 ) then
      ilast = ilast + 1
      if ( ilast <= lenl ) then
        last(ilast:ilast) = c
      end if
    end if

  end do

  return
end
