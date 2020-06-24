%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                   用来读出光盘中的数据
%               生成可以直接用于成像的原始数据
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 主程序：  read_data_Radarsat_1.m
% 该程序的目的： 
%   ——用来从光盘中读出数据；
%   ——经过下变频得到基带复信号，且转换成单精度浮点数；
%   ——进行AGC增益补偿；
%   ——最后，再转换为双精度浮点数
%   ——得到我们可以直接用于成像的原始数据，即为data。   
% 
% 关于得到的数据的说明：
%   其中， b = 1，是从分区1读出数据，我们记为 data_1
%          b = 2，是从分区2读出数据，我们记为 data_2
% 	利用 data_1 和 data_2 ，我们可以分别得到这两个分区的成像结果。
% 	此外，还可以将 data_1 和 data_2 合并为一步更大的数据，得到它们整体的成像结果。
%   
% 下面的程序就是得到这两个数据 data_1 和 data_2。
%
% read_data_Radarsat_1.m
% ----------------------------------------------------------
% 使用现成的程序‘compute.azim.spectra.m’中读出数据的方法；
% 利用函数 'load_DATA_block.m'，实现
%                - reads /loads data for a block 
%                - converts to floating point
%                - compansates for the receiver attenuation
% 变量 b -- 需要设置数据取自哪个分区
%                - b = 1 , from CDdata1
%                - b = 2 , from CDdata2
% 得到所需要的数据，也即可以直接进行后续 processing 的数据 data。
% ----------------------------------------------------------
%
% 该程序修改，截止到： 2014.10.10.  14:51 p.m.

%%
% ----------------------------------------------------------
% 得到可以进行后续信号处理的原始数据data（s_echo）
% ----------------------------------------------------------
%  Load the parameters for this run
clear;
clc;
close all;
format compact
set( 0, 'DefaultTextFontSize',   12 )  % Plotting defaults
set( 0, 'DefaultLineLineWidth', 1.5 )
set( 0, 'DefaultAxesFontSize',    8 )
tfs = 13;    lfs = 11;
load CD_run_params

% 设置b的值（ 1或者2 ），我们可以分别得到第一分区或者第二分区的数据。
b = 1;                      % 第一分区，即取自 CDdata1 数据。

file_pre = strcat( output_path, output_prefix, '_', num2str(b) );   
disp ' '
disp (['Load or Extract AGC setting and Data for block ' num2str(b) ])
%  Load a block of 'AGC_values'
AGC_values = load_AGC_block( file_pre, first_rg_line, ...
                                      Nrg_lines_blk, b , UseMATfiles );                                     
%  Load a block of raw SAR data
data = load_DATA_block( file_pre, output_path, Nrg_lines_blk, ...
                           Nrg_cells, AGC_values, b, UseMATfiles );
% 此时得到的 data 是经过下变频的基带复信号，且已经转换成单精度浮点数，也进行了AGC增益补偿。
data = double(data);    % 再将data转换成双精度浮点数，用于进行后续信号处理。
% 此时，我们就可得到了可以直接用于后续成像的原始数据：data

% 作图显示
figure;
imagesc(abs(data));
title('原始数据');     	% 原始回波数据（未处理）的幅度图像
% colormap(gray);




