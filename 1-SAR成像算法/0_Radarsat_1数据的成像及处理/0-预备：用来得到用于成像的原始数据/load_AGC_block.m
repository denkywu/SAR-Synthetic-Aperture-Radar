function [AGC_atten_dB] = load_AGC_block( file_pre, first_rg_line, ...
                                        Nrg_lines_blk, block, UseMATfiles )
%
% -------------------------------------------------------------------------
% This program returns the RADARSAT-1 receiver attenuation values for 
% one block.  They are integers from 2 to 17, in 1 dB steps.
% -------------------------------------------------------------------------
% file_pre        : File prefix ( path + prefix + block number)
% first_rg_line   : First range line for the extracted data  
% Nrg_lines_blk   : Number of range lines in a block
% block           : the current block number
% UseMATfiles     : (0) using extract_AGC and the '_aux.dat' file
%                   (1) load 'AGC_attenuation_values.mat',
%                          which has the AGC values for every range line.
% -------------------------------------------------------------------------
% Created:  Nov 01, 2004  by Kaan Ersahin
% Modified: Nov 29, 2004  by Kaan Ersahin
% Modified: Dec 13, 2004  by Ian Cumming
%            - added UseMATfiles option
% -------------------------------------------------------------------------

if     UseMATfiles == 0  % extract AGC values from binary aux data file
    
    file_in = strcat(file_pre, '_aux.dat');
    tic    
    AGC_atten_dB = extract_AGC( file_in, Nrg_lines_blk );  
    toc
    
elseif UseMATfiles == 1   % load all AGC values from .mat file
    tic    
    file_in = 'AGC_attenuation_values.mat'
    load( file_in )
    
    first = first_rg_line + (Nrg_lines_blk * (block-1));
    last  = first + Nrg_lines_blk - 1;
    AGC_atten_dB = nom_attenuation( first : last );
    toc
end

% Use check sum to verify read operation
% lenAGC = length(AGC_atten_dB);
% multfac = ([1:lenAGC] - lenAGC/2)';
% checksumAGC = sum( multfac .* AGC_atten_dB );
% fprintf('\nCheck sum of AGC values:%10.0f\n', checksumAGC )

%{
figure(201),   clf      %  Plot AGC values
plot( AGC_atten_dB ),   grid
axis([-9 length(AGC_atten_dB)+10  0.5 17.5])
xlabel('Azimuth (sample no. within this block)  \rightarrow', 'FontS', 12 )
ylabel('Magnitude  (dB)  \rightarrow', 'FontS', 12 )
title('AGC attenuation values for this block of data', 'FontS', 13 )
% text( lenAGC/2, 1.6, sprintf(...
%     '\nCheck sum of AGC values:%10.0f\n', checksumAGC ), 'Hor', 'c' )
%}
pause(0.1)