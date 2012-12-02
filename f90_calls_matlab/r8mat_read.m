function [ m, n, a ] = r8mat_read ( file_name )

%% R8MAT_READ reads an R8MAT from a file.
%
%  Discussion:
%
%    The file is assumed to contain the information about an M by N pointset,
%    with each row of data on a separate line of the file.
%
%    Comment lines may be included, usually at the beginning of the file.
%    A comment line begins with the '#' character.
%
%  Example:
%
%    # cvt_02_00010.txt
%    # created by CVT_DATASET
%    # at April 11 2003  12:04:56.303 PM
%    #
%    #  Spatial dimension M =   2
%    #  Number of points N =  10
%    #
%    #  Initial SEED =    123456789
%    #  Initialization by UNIFORM.
%    #  Sampling by UNIFORM.
%    #  Number of sample points =       500000
%    #  Number of sampling iterations =    100
%    #  L2 norm of dataset change on last step =   0.001501
%    #
%       0.168259       0.878328    
%       0.834417       0.833004    
%       0.521361       0.499896    
%       0.506248       0.165244    
%       0.180542       0.627410    
%       0.179467       0.372410    
%       0.505360       0.833925    
%       0.834464       0.166314    
%       0.841834       0.499935    
%       0.169745       0.122347    
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
%    Output, integer M, N, the order of the matrix.
%
%    Output, real A(M,N) is the data read from the file.
%
  a = [];

  file_unit = fopen ( file_name, 'r' );

  m = 0;
  n = 0;
  nline = 0;

  while ( 1 )

    line = fgets ( file_unit );

    if ( line == -1 )
      break;
    end 

    nline = nline + 1;
%
%  If the line begins with a "#" character, it's a comment, so skip it.
%
    if ( line(1) == '#' ) 
      continue;
    end

    m = m + 1;

    [ array, count ] = sscanf ( line, '%f' );

    if ( m == 1 ) 
      n = count;
    elseif ( n ~= count )
      fprintf ( 1, '\n' );
      fprintf ( 1, 'TABLE_READ - Fatal error!\n' );
      fprintf ( 1, '  Input lines should have %d values.\n', count );
      fprintf ( 1, '  But input line %d has %d.\n', nline, count );
      return;
    end

    a(m,1:n) = array(1:n)';
  
  end

  fclose ( file_unit );
