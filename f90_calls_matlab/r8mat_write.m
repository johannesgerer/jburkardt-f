function r8mat_write ( file_name, m, n, a )

%% R8MAT_WRITE writes an R8MAT to a file.
%
%  Modified:
%
%    02 November 2007
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, string FILE_NAME is the name of the file to read.
%
%    Input, integer M, N, the order of the matrix.
%
%    Input, real A(M,N), the matrix.
%
  file_unit = fopen ( file_name, 'w' );

  for i = 1 : m
    for j = 1 : n
      fprintf ( file_unit, '  %f', a(i,j) );
    end
    fprintf ( file_unit, '\n' );
  end

  fclose ( file_unit );
