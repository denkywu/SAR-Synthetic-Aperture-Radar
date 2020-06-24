% compute_azim_spectra.m
% ----------------------
%
% Divides the range swath (range cells) of this block into Nspect sections 
% and finds the average power spectrum for each segment.
%
% Fits a parabola to the baseband Doppler centroid vs. range
% This fit is most useful if the whole range swath is selected, but the
% computing times get large if more than 1024 range lines (azimuth cells)
% are selected in one run.
%
% Run "specify_run_parameters.m" and "extract_data.m" first to extract
% the data from the CD. The data can be stored in MAT or binary files.
% -------------------------------------------------------------------------

%  Load the parameters for this run
clear,    format compact
set( 0, 'DefaultTextFontSize',   12 )  % Plotting defaults
set( 0, 'DefaultLineLineWidth', 1.5 )
set( 0, 'DefaultAxesFontSize',    8 )
tfs = 13;    lfs = 11;
load CD_run_params

disp ' '
disp '---------------------------------------------------------'
fprintf(' UBC RRSG - Plot the azimuth spectrum of each data block')
disp ' '
disp '---------------------------------------------------------'

Nrowsg = 3;     % Number of subplots in row direction of the figure 
Ncolsg = 3;     % Number of subplots in column direction of the figure 
Nspect = Nrowsg*Ncolsg;         % Total number of spectra to calculate
Nrglpb = floor( Nrg_cells/Nspect );  %  No. of range cells per spectra
wd = 0.81/Ncolsg;   dx = wd + 0.045;   x0 = 0.07;    
ht = 0.39/Nrowsg;   dy = 0.28;         y0 = 1-dy;  

for b = 1 : Nblocks
    
    file_pre = strcat( output_path, output_prefix, '_', num2str(b) );
    
    disp ' '
    disp (['Load or Extract AGC setting and Data for block ' num2str(b) ])
    %  Load a block of 'AGC_values'
    AGC_values = load_AGC_block( file_pre, first_rg_line, ...
                                          Nrg_lines_blk, b , UseMATfiles );
                                      
    %  Load a block of raw SAR data
    data = load_DATA_block( file_pre, output_path, Nrg_lines_blk, ...
                             Nrg_cells, AGC_values, b, UseMATfiles );
        
    disp 'Compute and plot azimuth power spectra'
    tic
    figure(202), clf
    freq = [0:Nrg_lines_blk-1]*PRF/Nrg_lines_blk;
    b1 = first_rg_line + (b-1) * Nrg_lines_blk;
    b2 = first_rg_line + b * Nrg_lines_blk - 1;
    
    %  Find azimuth power spectra and average
    for krg = 1 : Nspect
        r1 = 1 + (krg-1)*Nrglpb;   r2 = r1 + Nrglpb - 1;
        DATA = fft( data(:,r1:r2) );
        DATA_aver(krg,:) = mean( abs(DATA.').^2 )/1000000;
        ysc(krg) = 1.05*max(DATA_aver(krg,:));
    end  % of for krg = 1 : Nspect
    ysc0 = max( ysc );      %  Common vertical scaling for all the spectra
    
    %  Plot the azimuth spectrum
    
    for krg = 1 : Nspect
        subplot(Nrowsg, Ncolsg, krg)
        plot( freq, DATA_aver(krg,:) ),   grid,   hold on
        set( gca, 'Pos',...
           [x0+dx*mod((krg-1),Ncolsg)  y0-dy*floor((krg-1)/Ncolsg)  wd ht])
        axis([0 PRF  0 ysc0]);
        
        azim_spec = fft( DATA_aver(krg,:) )/ Nrg_lines_blk;
        angle_first_harmonic = -angle( azim_spec(2) );
        Ffrac(krg) = PRF * angle_first_harmonic / (2*pi);
        if Ffrac(krg) < 0,   Ffrac(krg) = Ffrac(krg) + PRF;   end
        sine_fit = real(azim_spec(2)) * cos(2*pi*freq/PRF) - ...
                imag(azim_spec(2)) * sin(2*pi*freq/PRF) + 0.5*azim_spec(1);
        plot( freq, 2*sine_fit, 'r--' )
       
        r1 = 1 + (krg-1)*Nrglpb;   r2 = r1 + Nrglpb - 1;
        title(sprintf('RC %4d - %4d   Fdc =%6.0f',...
            r1+first_rg_cell-1, r2+first_rg_cell-1, Ffrac(krg) ), ...
            'FontS', lfs );
        if krg > Nspect - Ncolsg
            xlabel('Azimuth frequency  (Hz)  \rightarrow', 'FontS', lfs )
        end
        if mod(krg,Ncolsg) == 1
            ylabel('Power  \rightarrow', 'FontS', lfs )
        end
        if krg == 1
            text( 1.55*PRF, 1.7*double(ysc0), sprintf(...
             'Azimuth power spectra of range lines%6.0f  to%6.0f  and sine fit',...
              first_rg_line+(b-1)*Nrg_lines_blk, first_rg_line+b*Nrg_lines_blk-1 ),...
              'Hor', 'center', 'FontS', tfs )
        end
        pause(0.1)
    end  % of for krg = 1 : Nspect
    toc
    
    %  Plot Ffrac vs. range
    
    figure(203),   clf
    plot( Ffrac, 'bx-', 'MarkerS', 9 ),   grid,   hold on
    eval(['save Ffrac_' num2str(first_rg_line) ' Ffrac'])
    axis([0.5 Nspect+0.5  200 800])
    xlabel('Range swath number', 'FontS', lfs+1 )
    ylabel('Frequency  (Hz)', 'FontS', lfs+1 )
    title(sprintf(...
        'Baseband Doppler centroid over%5.0f lines starting at%6.0f',...
        Nrg_lines_blk, first_rg_line+(b-1)*Nrg_lines_blk ), 'FontS', tfs )
    
    coeff = polyfit( [1:Nspect], Ffrac, 2 );
    Fit = polyval( coeff, [1:Nspect] );
    plot( Fit, 'r-' )
    text( 0.79*Nspect, 750, sprintf('C_2 =%6.2f', coeff(1) ), 'FontS',13 )
    text( 0.79*Nspect, 650, sprintf('C_1 =%7.1f', coeff(2) ), 'FontS',13 )
    text( 0.79*Nspect, 550, sprintf('C_0 =%5.0f', coeff(3) ), 'FontS',13 )
    text( 0.13*Nspect, 250, sprintf(...
        'Range cells %3.0f to %4.0f', first_rg_cell, ...
         first_rg_cell+Nrg_cells-1 ), 'FontS',13 )
    pause(0.1)
  % Save the plot as .eps file 
  % file_eps=strcat(output_path,output_prefix,'azimuth_',num2str(b),'.eps')
  % saveas(gcf,file_eps,eps)
end  % of for b = 1 : Nblocks
beep,   pause(0.3),   beep,   pause(0.3),   beep