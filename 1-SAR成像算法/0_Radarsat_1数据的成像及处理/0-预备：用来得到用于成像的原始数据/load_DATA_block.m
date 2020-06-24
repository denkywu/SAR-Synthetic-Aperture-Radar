function [data] = load_DATA_block( file_pre, output_path, ...
            Nrg_lines_blk, Nrg_cells, AGC_atten_dB, block, UseMATfiles )
%
% This program - reads /loads data for a block 
%              - converts to floating point
%              - compansates for the receiver attenuation
% -------------------------------------------------------------------------
% file_pre      : File prefix (path + prefix + block number)
% output_path   : Path where the SAR data MAT files are stored
% Nrg_lines_blk : Number of range lines in a block
% Nrg_cells     : Number of range cells
% AGC_atten_dB  : Nominal attenuation (Rx_AGC_setting) for each line (dB)
% block         : block number
% UseMATfiles   : Data is stored in MAT files
% -------------------------------------------------------------------------
% Created:  Nov 01, 2004  by Kaan Ersahin
% Modified: Nov 29, 2004  by Kaan Ersahin
% Modified: Dec  4, 2004  by Ian Cumming
%           - added single precision option for MATLAB 7
%           - added UseMATfiles option
%           - moved AGC correction to this function
% -------------------------------------------------------------------------

vv        = version;   vers = str2num(vv(1));  % Find MATLAB version

%  Read 4-bit unsigned data as single (MATLAB 7) or double
%  The data is stored as [I Q I Q I Q...] 4-bit words

if UseMATfiles
    % fprintf('\nload_DATA_block:  from MAT file  CDdata%1.0f.mat\n', block )
    disp ' '
    mat_file = strcat( output_path, 'CDdata', num2str(block) )
    load( mat_file )   % The data is already decoded in this case
    
else   %  Use fread to get all the data and parameters
    % fprintf('\nload_DATA_block:  from binary files\n')
    
    file_in  = strcat( file_pre, '_data.dat' )
    fid2     = fopen( file_in, 'r' );
    if vers < 7
      data = fread( fid2, [2*Nrg_cells, Nrg_lines_blk], 'ubit4=>double' );
    else
      data = fread( fid2, [2*Nrg_cells, Nrg_lines_blk], 'ubit4=>single' );
    end
    data = data';  %  Arrange the range lines as rows of the matrix
    
    %  Compensate for packed data format --> convert to signed numbers
    data = 2*( data - 16*(data > 7)) + 1;

    %  Separate the I and Q channels and make a complex array
    data = complex( data(:, 1:2:2*Nrg_cells), data(:, 2:2:2*Nrg_cells) );
    fclose(fid2);
end

% Apply gain correction to the data.  The attenuation varies from 2 to 17
% dB for this CD, so the linear gain factor varies from 1.26 to 7.08.
% If you want to store the decoded data in one byte arrays, use an additional 
%  factor of 1.5, so that the maximum abs value is less than 127.
fact = 1.5; 
if vers < 7      %  For MATLAB version 6
    linear_gain_factor = double( fact * 10.^(AGC_atten_dB/20) );
    ones_array         = ones(1,Nrg_cells);
    data = (linear_gain_factor * ones_array) .* double(data);
else             %  For MATLAB version 7
    linear_gain_factor = single( fact * 10.^(AGC_atten_dB/20) );
    ones_array         = single(ones(1,Nrg_cells));
    data = (linear_gain_factor * ones_array) .* single(data);
end
