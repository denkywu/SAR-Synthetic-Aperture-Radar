function [AGC_atten_dB] = extract_AGC( file_in, Nlines )
%
% -------------------------------------------------------------------------
% This program reads the RADARSAT aux data and calculates the AGC settings
% -------------------------------------------------------------------------
%  file_in      : input file name including it's path
%  Nlines       : number of azimuth lines
%  AGC_atten_dB : Receiver attenuation in dB for each range line
% -------------------------------------------------------------------------
% Created  : May 05, 2004 by Kaan Ersahin & Millie Sikdar
% Modified : Nov 22, 2004 by Kaan Ersahin
% -------------------------------------------------------------------------
 

% The parameters encoded in the auxilary data bits are defined in RSI-D6
% also known as RSCSA-IC0009 (X-band ICD) 
% -------------------------------------------------------------------------
% PARAMETER NAME                          LOCATION         LENGTH   ID 
% -------------------------------------------------------------------------
% aux_sync_marker         = aux_bits (:,   1:  32);     % 32 bit -  1 
% image_ref_id            = aux_bits (:,  33:  64);     % 32 bit -  2 
% payload_status          = aux_bits (:,  65:  80);     % 16 bit -  3 
% replica_AGC             = aux_bits (:,  81:  86);     %  6 bit -  4 
% CALN_atten_LPT_pow_set  = aux_bits (:,  89:  96);     %  8 bit -  6 
% pulse_waveform_number   = aux_bits (:,  97: 100);     %  4 bit -  7 
% temperature             = aux_bits (:, 113: 144);     % 32 bit -  9 
% beam_sequence           = aux_bits (:, 145: 160);     % 16 bit - 10 
% ephemeris               = aux_bits (:, 161: 176);     % 16 bit - 11
% number_of_beams         = aux_bits (:, 177: 178);     %  2 bit - 12 
% ADC_rate                = aux_bits (:, 179: 180);     %  2 bit - 13 
% pulse_count_1           = aux_bits (:, 185: 192);     %  8 bit - 15
% pulse_count_2           = aux_bits (:, 193: 200);     %  8 bit - 16
% PRF_beam                = aux_bits (:, 201: 213);     % 13 bit - 17 
% beam_select             = aux_bits (:, 214: 215);     %  2 bit - 18  
% Rx_window_start_time    = aux_bits (:, 217: 228);     % 12 bit - 20 
% Rx_window_duration      = aux_bits (:, 233: 244);     % 12 bit - 22 
% altitude                = aux_bits (:, 249: 344);     % 96 bit - 24 
% time                    = aux_bits (:, 345: 392);     % 48 bit - 25 
% SC_T02_defaults         = aux_bits (:, 393: 393);     %  1 bit - 26 
% first_replica           = aux_bits (:, 394: 394);     %  1 bit - 27 
% Rx_AGC_setting          = aux_bits (:, 395: 400);     %  6 bit - 28 
% -------------------------------------------------------------------------
%                                    Total  => 50 bytes (400 bits)              
% -------------------------------------------------------------------------



% -------------------------------------------------------------------------
% read the binary aux data file
% -------------------------------------------------------------------------
file_in
fid_aux  = fopen(file_in,'r');
aux_bytes = fread(fid_aux,[50, Nlines],'uint8'); 
aux_bytes = aux_bytes';
fclose(fid_aux);

aux_bits = char(zeros(Nlines,400));
% -------------------------------------------------------------------------
% convert 50 bytes to 400 bits in general (it is only done for 50th byte)
% -------------------------------------------------------------------------
for byte = 50:50 %1:50 %% in general
   aux_bits(:,8*(byte-1)+1:8*byte) = num2str(dec2bin(aux_bytes(:,byte),8));
end
% -------------------------------------------------------------------------
%  Convert last 6 bits from binary to decimal
% -------------------------------------------------------------------------
d_Rx_AGC = bin2dec( aux_bits (:, 395: 400) );

% -------------------------------------------------------------------------
% For values greater than 31, the binary representation defined in the 
% documentation is different and there is the need to substract 24 from 
% decimal values to get the AGC setting in dB
% -------------------------------------------------------------------------
AGC_atten_dB = d_Rx_AGC - ( 24 * (d_Rx_AGC > 31));  
