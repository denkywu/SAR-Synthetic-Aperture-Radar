% specify_parameters.m
%
% Run this program first to set up parameters for each data run.  The run
% parameters are stored in the file "CD_run_params.mat" in the current
% directory.  The program "extract_data.m" reads this file and adds some
% fixed radar parameters to the file.  Subsequent programs can read this
% file to obtain the parameters needed for the SAR processing.
%
% The radar data are stored on the CD supplied with the book, but they 
% can be read from a disk file copied from the CD to your hard drive.
%
% The only file needed from the CD is "dat_01.001", which contains the
% radar data as well as the chirp replica and auxiliary data for each range
% line.  As none of the parameters like PRF and range gate delay change in
% this data set, and the replica is well defined by a linear FM chirp, the
% only data that are really needed from the CD is the raw radar signal data
% contained in "dat_01.001".  The radar parameters are stored in the
% program "extract_data.m".
%
% The radar data can be written in one block, or broken up into a number of
% blocks, with a set of files for each block.
%
% The first range cell and the number of range cells can be set to any
% value within the boundaries of the data array.  The first range line
% should be set to 1 plus a multiple of 8, and the number of lines per
% block set to a factor of 8.  Also, the number of lines should be set to a
% multiple of the number of lines per block.  If these rules are not
% obeyed, the "extract_data.m" program will quantize the values accordingly.

% Parameters to specify for each run:
% -----------------------------------
% input_path    : location of the radar data file 'dat_01.001' 
% output_path   : location for the extracted data files
% output_prefix : prefix to the extracted data file names
%
% first_rg_cell : first range cell to extract
% Nrg_cells     : number of range cells to extract
% first_rg_line : first range line (azimuth sample) to extract
% Nrg_lines     : number of range lines to extract
% Nrg_lines_blk : number of range lines per block
% UseMATfiles   : Set to 0 if you want binary header, data, aux and replica 
%                 files (.DAT) rather than just the MATLAB data files (.MAT)
% -------------------------------------------------------------------------
% Created :   Nov 29, 2004 by Kaan Ersahin 
% Modified:   Nov 30, 2004 by Kaan Ersahin
% Modified:   Dec 3,  2004 by Ian Cumming 
%             - moved fixed parameters to "extract_data.m"
%             - now use a fixed parameter file name for all programs
%             - added UseMATfiles option (=1 for MAT files)
% -------------------------------------------------------------------------

clear,   home,   format compact

% Define location of data files
% input_path      = 'Y:\RSAT\scene01\'; % These values are for runs at UBC
% output_path     = 'Y:\RSAT\EXTRACTED_DATA\';
%input_path        = 'C:\a_RADARSAT_Processing\Read_the_CD\RAW_CD\scene01\';
input_path        = 'F:\00_Radarsat_1数据的成像及处理\scene01\';
%output_path       = 'C:\a_RADARSAT_Processing\Read_the_CD\EXTRACTED_DATA\';
output_path       = 'F:\00_Radarsat_1数据的成像及处理\EXTRACTED_DATA\';
CD_data_file_name = 'dat_01.001';
output_prefix     = 'raw';

%==========================================================================
%  Define area of interest of the radar data to be extracted
%  Note that the azimuth parameters are quantized in "extract_data.m"
%
%  Here are some possible values for first_rg_cell and first_rg_line:
%  For coal & ferry terminals, use range =  970, azimuth =  1801
%  For Vancouver airport,      use range = 1060, azimuth =  5561
%  For UBC and one ship,       use range =  760, azimuth =  7097
%  For Stanley Park & city,    use range = 1850, azimuth =  7657
%  For English Bay ships,      use range = 1050, azimuth =  7769
%  For Squamish & Garibaldi,   use range = 2640, azimuth = 16169
%  For Brackendale,            use range = 2800, azimuth = 17897 (max az)

first_rg_cell   =  1;%1050;    % Define the range limits
first_rg_line   =  1;%7769;    % Define the azimuth limits  (19432 max)

Nrg_cells       =  9288 %2048     % Suggest 2048 cells
Nrg_lines_blk   =  19424/2;%6*256;   % Suggest 1536, it should be larger than the
                            %        size of the azimuth match filter (705) 
                            %              so that an image can be created.                    
Nrg_lines   =  2*Nrg_lines_blk;

UseMATfiles =  1; % (1) Use   MAT  files to save the extracted data
                  % (0) Use binary files to save the extracted data
SaveV6      =  0; % (1) Save data in MATLAB v6 format when using v7
                  % (0) Save data in MATLAB v7 format
                  % This variable is only used if you are running MATLAB 7
%==========================================================================

%  Save the data and run parameters in a matlab data file (*.mat)

vv = version;   vers = str2num(vv(1));  % Find MATLAB version
if vers == 6,   SaveV6 = 0;   end

if SaveV6
  save -v6 CD_run_params
  fprintf('\nData parameters written to:  CD_run_params.mat  in V6 format')
else
  save CD_run_params
  fprintf('\nData parameters written to:  CD_run_params.mat')
end

fprintf('\nNow run "extract_data.m"\n\n')
beep,   pause(0.4),   beep

%% This setup program is now finished, but you can add calls to subsequent 
%%  data extraction and processing programs below, if desired.

% extract_data           % Extract the SAR data from the CD

% compute_azim_spectra   % Find the baseband Doppler centroid

% compress_data          % Perform SAR processing
