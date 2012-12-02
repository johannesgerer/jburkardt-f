%% MATRIX_MULTIPLY is a script to read A and B and compute C = A * C.
%
%  Discussion:
%
%    This program is designed to be used in batch mode.  In particular,
%    once it is done, it issues an EXIT statement, which terminates
%    the execution of MATLAB.
%
%  Modified:
%
%    02 November 2007
%
%  Author:
%
%    John Burkardt
%

%
%  We need to move MATLAB's default directory to the place where the
%  FORTRAN program has place the files, and where the M-files
%  R8MAT_READ.M and R8MAT_WRITE.M are available.
%
  cd /Users/burkardt/public_html/f_src/f90_matlab

  fprintf ( 1, '\n' );
  fprintf ( 1, 'MATRIX_MULTIPLY\n' );
  fprintf ( 1, '  A MATLAB script which reads matrices A and B from files,\n' );
  fprintf ( 1, '  multiplies them, and writes the result C to a file.\n' );

  a_file_name = 'a.txt';
  [ n1, n2, a ] = r8mat_read ( a_file_name );

  fprintf ( 1, '\n' );
  fprintf ( 1, '  Read A matrix of order (%d,%d).\n', n1, n2 );

  b_file_name = 'b.txt';
  [ n3, n4, b ] = r8mat_read ( b_file_name );

  fprintf ( 1, '\n' );
  fprintf ( 1, '  Read B matrix of order (%d,%d).\n', n3, n4 );
%
%  Make sure matrices are conformable.
%
  if ( n2 ~= n3 )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'MATRIX_MULTIPLY - Fatal error!\n' );
    fprintf ( 1, '  Matrices A and B are not conformable.\n' );
    fprintf ( 1, '  A = (%d,%d)\n', n1, n2 );
    fprintf ( 1, '  B = (%d,%d)\n', n3, n4 );
    exit;
  end

  c = a * b;

  c_file_name = 'c.txt';
  r8mat_write ( c_file_name, n1, n4, c );

  fprintf ( 1, '\n' );
  fprintf ( 1, '  Wrote C matrix of order (%d,%d).\n', n1, n4 );

  fprintf ( 1, '\n' );
  fprintf ( 1, 'MATRIX_MULTIPLY:\n' );
  fprintf ( 1, '  Normal end of execution.\n' );
%
%  Terminate the MATLAB session.
%
  exit;
