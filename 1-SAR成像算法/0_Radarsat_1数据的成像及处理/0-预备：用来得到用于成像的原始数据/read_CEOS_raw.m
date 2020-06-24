function read_ceos_raw ( output_file_pre, start_line_blk, blk )
%
% -------------------------------------------------------------------------
% This function reads the CEOS format RSAT-1 RAW data from the CD
% This function is called by extract_data.m
% -------------------------------------------------------------------------
% output_file_pre : output file prefix ( path + prefix + block number)
% start_line_blk  : extract data starting from row number 'start_line_blk'
% blk             : the block number
%
% input_path     : Input file path
% length_replica : Length (I & Q) of replica (for Fine Beam -> 2880 bytes)
% tot_Nrg_cells  : Total number of columns available in the dataset
% Nrg_lines_blk  : Number of samples in azimuth
% first_rg_cell  : Extract data starting from column number 'first_rg_cell'
% Nrg_cells      : Number of samples in ranged
% UseMATfiles    : Do you want to write MAT files for the data rather than
%                  binary files for all the data and parameters?
% -------------------------------------------------------------------------

% date created  : May 05, 2004  by  Kaan Ersahin & Millie Sikdar
% date modified : Nov 30, 2004  by  Kaan Ersahin
% Modified:  Dec 17 by Ian Cumming
%             - Added UseMATfiles option
% -------------------------------------------------------------------------

% 'RSI - D4' refers to the document RSI-GS-026, Rev 3.0 - May 8, 2000
%  that describes the detailed format of the RADARSAT-1 raw data CDs.
tic
num_aux            = 50;                        % RSI - D4 Pg.32
num_header         = 192;                       % RSI - D4 Pg.32
file_header_length = 16252;                     % RSI - D4 Pg.103
load CD_run_params

%  Specify input data file path and open to read
file_in  = strcat( input_path, CD_data_file_name );
fid      = fopen( file_in, 'r' );

%  Create output file names and open to write if MAT files not used
if ~UseMATfiles
    file_replica 	 = '_replica.dat';     % one replica per 8-line
    file_data_header = '_data_header.dat'; % header of each line
    file_aux     	 = '_aux.dat';         % auxiliary data from each line
    file_data    	 = '_data.dat';        % radar signal data
    file_file_header = '_file_header.dat'; % header of file 'dat_01.001'
    
    file_replica     = strcat(output_file_pre, file_replica);
    file_data_header = strcat(output_file_pre, file_data_header);
    file_aux         = strcat(output_file_pre, file_aux);
    file_data        = strcat(output_file_pre, file_data);
    file_file_header = strcat(output_file_pre, file_file_header);
    
    fod1             = fopen( file_replica    , 'w' );
    fod2             = fopen( file_data_header, 'w' );
    fod3             = fopen( file_aux        , 'w' );
    fod4             = fopen( file_data       , 'w' );
    fod5             = fopen( file_file_header, 'w' );

    %  Read file header
    file_header      = fread( fid, [file_header_length,1], 'uint8' );
    %  Write file header in a binary file as unsigned integers
    fwrite( fod5, file_header, 'uint8');
    fclose( fod5 );
end  % of if ~UseMATfiles

% some useful calculations
num_pixel_data   = 2*tot_Nrg_cells;                               %  18576
h                = num_header;                                    %    192
ha               = num_header + num_aux;                          %    242
hpa              = num_header + num_aux + num_pixel_data;         %  18818

% Define 'rep_block' as the 8-line block of data including one replica
% 'start_rep' and 'end_rep' locates the replica in the 'rep_block'
% 'start_col' and 'end_col' locates the area of interest in the columns
rep_block_length = length_replica + 8*hpa ;

start_rep        = (first_replica-1)*hpa + num_header + num_aux + 1; 
end_rep          = start_rep + length_replica - 1 ;

start_col        = ha + 2*(first_rg_cell - 1) + 1 ;
end_col          = start_col + 2*Nrg_cells - 1;

% The number of bytes to skip from the beginning of the file (Nbytes_skip)
Nrep_blks_row1   = ceil(start_line_blk/8) - 1;
Nbytes_blocks    = Nrep_blks_row1 * rep_block_length;
Nbytes_skip      = file_header_length + Nbytes_blocks;

fseek( fid, Nbytes_skip , 'bof' );  % go to start of the 1st data block
data = int8( zeros(Nrg_lines_blk, 2*Nrg_cells) );   % Allocate data array

N_8line_blocks = ceil(Nrg_lines_blk/8);
fprintf('\nReading %1.0f small 8-line blocks from range line %5.0f\n',...
    N_8line_blocks, start_line_blk )


for kb = 1 : N_8line_blocks   %  Read and decode 8 lines at a time
    
    % read one 8-line block of radar data and replica
    temp = uint8( fread( fid, [rep_block_length,1], 'uint8' ) );
    
    % separate the replica and 8 lines of [header, aux and data]
    replica = temp(start_rep : end_rep)';
    temp = [ temp(1 : start_rep-1)' temp(end_rep+1 : rep_block_length)' ];
    temp = reshape(temp, hpa, 8)';   % Shape array as 8 rows by all columns

    % Extract the desired set of range cells
    data(1+8*(kb-1):8*kb, :) = temp( : , start_col : end_col );

    if ~UseMATfiles  % Write replica, header, aux and data in a binary file
        header (:,1:8) = temp( : ,         1 : h       )';
        aux    (:,1:8) = temp( : ,       h+1 : ha      )';
        count1 = fwrite( fod1, replica, 'ubit4' ); % 4-bit replica
        count2 = fwrite( fod2, header , 'uint8' ); % 8-bit header
        count3 = fwrite( fod3, aux    , 'uint8' ); % 8-bit aux data
        count4 = fwrite( fod4, data(1+8*(kb-1):8*kb, :)'   , 'ubit4' );
    end  % of if UseMATfiles
end

fclose(fid);           % Close input file
if ~UseMATfiles        % Close output files
    fclose(fod1);   fclose(fod2);   fclose(fod3);   fclose(fod4);
    fprintf('This block written to binary file:\n%s\n', file_data )
end

%--------------------------------------------------------------------------
%  Decode the data, save as MAT file , if UseMATfiles == 1
%--------------------------------------------------------------------------

if UseMATfiles == 1
    disp 'Decode data by subtracting 16 from upper register'
    vv    = version;   vers = str2num(vv(1));  % Find MATLAB version
    seven = int8( 7 );
        
    if vers == 7
        one     = int8( 1 );
        two     = int8( 2 );
        sixteen = int8( 16 );
        data    = two*(data - sixteen*int8(data > seven)) + one;
    else
        data    = int8( 2*(double(data) - 16*double(data > seven) ) + 1 );
    end
    
    %  Separate the I and Q channels and form a complex array
    data = complex( data(:, 1:2:2*Nrg_cells), data(:, 2:2:2*Nrg_cells) );
    
    if SaveV6
       eval(['save -v6 ' output_path 'CDdata' num2str(blk) ' data'])
       disp  'Save int8 data in MAT file in V6 format'
    else
       eval(['save ' output_path 'CDdata' num2str(blk) ' data'])
       fprintf('Save int8 data in file:  CDdata%1.0f.mat\n', blk)
    end
    % somedata = data(1,1:4)  % Check a few data samples
end   % of if UseMATfiles
toc