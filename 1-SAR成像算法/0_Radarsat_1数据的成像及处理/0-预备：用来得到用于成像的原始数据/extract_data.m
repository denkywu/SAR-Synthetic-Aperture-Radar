% extract_data.m
% --------------
% This program extracts the raw signal data from the CD included in the 
% Cumming/Wong book.  The data on CD are in CEOS format.  
% It is assumed that the run parameters are stored in CD_run_params.mat 
% in the current directory.
% Run "specify_run_parameters.m" first to create this file.
% -------------------------------------------------------------------------
% Created :   Nov 01, 2004  by Kaan Ersahin
% Modified:   Nov 22, 2004  by Kaan Ersahin
% Modified:   Dec 3,  2004  by Ian Cumming
%              - changed function to m-file
%              - fixed the parameter file name
%              - added more radar parameters
% -------------------------------------------------------------------------

clear,   home,   format compact

%  Load the input parameters from a matlab data file 
load CD_run_params

disp ' '
disp '------------------------------------------------'
disp ' UBC RRSG - CEOS Reader for RADARSAT-1 RAW DATA'
disp '------------------------------------------------'
disp ' '

% -------------------------------------------------------------------------
% Quantize the range line limits and block size, if necessary
% -------------------------------------------------------------------------

% Move the first range line to the beginning of an 8-line block 
first_rg_line = 8 * ( ceil(first_rg_line / 8) - 1 ) + 1;

% Make 'Nrg_lines_blk' a multiple of 8, to get complete the 8-line blocks
Nrg_lines_blk = 8 * ceil(Nrg_lines_blk / 8); 

% Find the number of complete blocks required to cover the area of interest
Nblocks = ceil(Nrg_lines / Nrg_lines_blk);

% Make 'Nrg_lines' a multiple of 'Nblocks', to get complete blocks
Nrg_lines = Nrg_lines_blk * Nblocks; 

% =========================================================================
% These values are specific to the data set, DO NOT CHANGE for this CD
% =========================================================================

length_replica  =  2880;         % Total length (I&Q) of replica record
tot_Nrg_cells   =  9288;         % Total number of range cells per line
tot_Nrg_lines   = 19432;         % Total number of range lines (records)
first_replica   =     7;         % First record that contains the replica 
PRF             = 1256.98;       % Pulse Reputation Frequency (Hz)
Fr              = 32.317e+6;     % Radar sampling rate (Hz)
f0              = 5.300e+9;      % Radar center frequency (Hz)
c               = 2.9979e+8;     % Speed of light (m/s)
R0              = 0.0065956*c/2; % Slant range of first radar sample (m)
Nrepl           = 1349;          % No. of valid samples in the replica
Kr              = 0.72135e+12;   % FM rate of radar pulse (Hz/s)
Tr              = 41.75e-6;      % Chirp duration (s)

% -------------------------------------------------------------------------
% Save parameters in a MAT file that can be used by subsequent programs

if     SaveV6,   save -v6 CD_run_params
else,            save     CD_run_params,      end
% -------------------------------------------------------------------------

fprintf('Total number of range lines      : %5d \n', tot_Nrg_lines )
fprintf('Total number of range cells      : %5d \n', tot_Nrg_cells )
fprintf('First range cell to be extracted : %5d \n', first_rg_cell )
fprintf('First range line to be extracted : %5d \n', first_rg_line )
fprintf('Number of range cells            : %5d \n', Nrg_cells )
fprintf('Number of range lines            : %5d \n', Nrg_lines )
fprintf('Number of range lines per block  : %5d \n', Nrg_lines_blk )
fprintf('Number of blocks                 : %5d \n', Nblocks )
disp ' '
disp '------------------------------------------------'

% -------------------------------------------------------------------------
%  Check the dimensions of the selected area to be processed.
% -------------------------------------------------------------------------

if (first_rg_line <= 0) | ((first_rg_line + Nrg_lines - 1) > tot_Nrg_lines)
    disp ' ',  disp '*****************************************************'
    disp ' ERROR: Check the limits of the range lines !', beep, return
end
if (first_rg_cell <= 0) | ((first_rg_cell + Nrg_cells - 1) > tot_Nrg_cells)
    disp ' ',  disp '*****************************************************'
    disp ' ERROR: Check the limits of the range cells !', beep, return
end
    
% -------------------------------------------------------------------------
% EXTRACT DATA from the area of interest and write data files  
% -------------------------------------------------------------------------

for blk = 1 : Nblocks
    % find the first range line of block 'blk' 
    start_line_blk = Nrg_lines_blk * (blk-1) + first_rg_line;
    fprintf('\nExtracting block number%3.0f,  RC1 =%5.0f,  RL1 =%6.0f\n',...
        blk, first_rg_cell, start_line_blk )
    
    % create the output file name for block 'blk'
    output_file_pre = strcat(output_path,output_prefix,'_',num2str(blk));

    % Call 'read_ceos_raw' function to extract the data for block 'blk'
    read_CEOS_raw( output_file_pre, start_line_blk, blk );

end  % of the 'blk' for loop

beep,   pause(0.3),   beep
