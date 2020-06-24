%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                       Radarsat_1 光盘中数据
%                             CSA 成像
%
%
%                               WD
%                       2014.10.19. 13:53 p.m.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 程序说明：
% 主程序是：  Radarsat_1_CSA.m
%
% （1）原始数据说明：
% 文件夹中的 data_1 和 data_2 是已经经过下列方法得到的原始数据，
% 可以直接进行后续成像
% ----------------------------------------------------------
% 使用现成的程序‘compute.azim.spectra.m’中读出数据的方法；
% 利用函数 'laod_DATA_block.m'，实现
%                - reads /loads data for a block 
%                - converts to floating point
%                - compansates for the receiver attenuation
% 变量 b -- 需要设置数据取自哪个分区
%                - b = 1 , from CDdata1
%                - b = 2 , from CDdata2
% 得到所需要的数据，也即可以直接进行后续 processing 的数据 data。
% ----------------------------------------------------------
% 因此，文件夹中的 data_1和data_2 分别是分区1和分区2的数据，经过了下变频，
%       转换为了浮点数，进行了AGC增益补偿，最后转换为了double双精度浮点数。
%       因此，直接载入这两个数据就可以进行后续成像。
%
% （2） 本文件夹中还有一个文件：CD_run_params
%           ——这里面是仿真中需要用的许多参数，直接载入即可。
%
% （3）成像程序说明：
%       由CSA的点目标程序修改而来；
% （4）成像流程：
%   ——原始数据
%   ——经过方位向FFT，变换到距离多普勒域，进行“补余RCMC”
%   ——经过距离向FFT，变换到二维频域，进行“距离压缩”、“SRC”、“一致RCMC”
%   ——经过距离向IFFT，变换到距离多普勒域，进行“方位压缩”和“附加相位校正”
%   ——经过方位向IFFT，回到图像域，成像结束。
%
% 本程序修改截止到： 2014.10.19. 13:53 p.m.
%
% 注：修改后的程序中，主要是附加了一步：对原始数据进行补零，再进行后续处理。

%%
clear;
clc;
close all;
% ----------------------------------------------------------
% 得到可以进行后续信号处理的原始数据data（s_echo）
% ----------------------------------------------------------
% 载入参数
load CD_run_params;

% 载入数据
b = 1;              % 选择对于哪一部分成像
% b = 1，则对分区1成像
% b = 2，则对分区2成像
% b = 3，则对整个数据（分区1和分区2）成像

if b == 1
    load data_1;                % 分区1的数据
    s_echo = data_1;            % 原始数据记为s_echo，用于后续成像。
end
clear data_1;                   % 清除data_1，以腾出内存

if b == 2
    load data_2;                % 分区2的数据
    s_echo = data_2;            % 原始数据记为s_echo，用于后续成像。
end
clear data_2;                   % 清除data_2，以腾出内存

if b == 3
    load data_1;                % 分区1的数据    
    s_echo1 = data_1;
    load data_2;                % 分区2的数据
    s_echo2 = data_2;
    s_echo = [s_echo1;s_echo2]; % 将分区1和分区2的数据合成整个数据块，用于成像
end
clear data_1;clear data_2;clear s_echo1;clear s_echo2;

%{
% 作图显示
figure;
imagesc(abs(s_echo));
title('原始数据');              % 原始回波数据（未处理）的幅度图像
% colormap(gray);
%}

%%
% --------------------------------------------------------------------
% 定义一些参数
% --------------------------------------------------------------------
Kr = -Kr;                       % 将调频率Kr改成负值
BW_range = 30.111e+06;          % 脉冲宽度
Vr = 7062;                      % 有效雷达速率
Ka = 1733;                      % 方位调频率
fnc = -6900;                    % 多普勒中心频率
Fa = PRF;                       % 方位向采样率
lamda = c/f0;                   % 波长
T_start = 6.5959e-03;           % 数据窗开始时间

Nr = round(Tr*Fr);              % 线性调频信号采样点数
Nrg = Nrg_cells;                % 距离线采样点数
if b == 1 || b == 2
    Naz = Nrg_lines_blk;     	% 每一个数据块的距离线数
else
    Naz = Nrg_lines;          	% 两个数据块，总共的距离线数
end
NFFT_r = Nrg;                   % 距离向FFT长度
NFFT_a = Naz;                   % 方位向FFT长度

R_ref = R0;                     % 参考目标选在场景中心，其最近斜距为 R_ref  
fn_ref = fnc;                   % 参考目标的多普勒中心频率

%%
%
% --------------------------------------------------------------------
% 对原始数据进行补零
% --------------------------------------------------------------------
if b == 1 || b == 2 
    data = zeros(1*2048,3000);
else
    data = zeros(2*2048,3000);
end
data(1:Naz,1:Nrg) = s_echo;
clear s_echo;
s_echo = data;
clear data;
[Naz,Nrg] = size(s_echo);

NFFT_r = Nrg;               	% 距离向FFT长度
NFFT_a = Naz;                   % 方位向FFT长度

% 作图显示
figure;
imagesc(abs(s_echo));
title('补零后的原始数据');       % 补零后的原始回波数据（未处理）的幅度图像
%}

%%
% --------------------------------------------------------------------
% 距离（方位）向时间，频率相关定义
% --------------------------------------------------------------------
% 距离
tr = 2*R0/c + ( -Nrg/2 : (Nrg/2-1) )/Fr;                % 距离时间轴
fr = ( -NFFT_r/2 : NFFT_r/2-1 )*( Fr/NFFT_r );          % 距离频率轴
% 方位
ta = ( -Naz/2: Naz/2-1 )/Fa;                            % 方位时间轴
fa = fnc + fftshift( -NFFT_a/2 : NFFT_a/2-1 )*( Fa/NFFT_a );	% 方位频率轴

% 生成距离（方位）时间（频率）矩阵
tr_mtx = ones(Naz,1)*tr;    % 距离时间轴矩阵，大小：Naz*Nrg
ta_mtx = ta.'*ones(1,Nrg);  % 方位时间轴矩阵，大小：Naz*Nrg
fr_mtx = ones(Naz,1)*fr;    % 距离频率轴矩阵，大小：Naz*Nrg
fa_mtx = fa.'*ones(1,Nrg);  % 方位频率轴矩阵，大小：Naz*Nrg

%%
% --------------------------------------------------------------------
% 变换到距离多普勒域，进行“补余RCMC”
% --------------------------------------------------------------------
s_rd = s_echo.*exp(-1j*2*pi*fnc.*(ta.'*ones(1,Nrg))); 	% 数据搬移
S_RD = fft(s_rd,NFFT_a,1);  % 进行方位向傅里叶变换，得到距离多普勒域频谱

D_fn_Vr = sqrt(1-lamda^2.*(fa.').^2./(4*Vr^2));     % 大斜视角下的徙动因子，列向量
D_fn_Vr_mtx = D_fn_Vr*ones(1,Nrg);  % 形成矩阵，大小：Nrg*Naz

D_fn_ref_Vr = sqrt(1-lamda^2*fn_ref^2/(4*Vr^2));    % 参考频率fn_ref处的徙动因子，是常数。

K_src = 2*Vr^2*f0^3.*D_fn_Vr.^3./(c*R_ref*(fa.').^2);   % 列向量，使用R_ref处的值 
K_src_mtx = K_src*ones(1,Nrg);  % 形成矩阵
Km = Kr./(1-Kr./K_src_mtx);     % 矩阵，这是变换到距离多普勒域的距离调频率。
                                % 使用 R_ref 处的值

% 下面生成 变标方程 s_sc
s_sc = exp(1j*pi.*Km.*(D_fn_ref_Vr./D_fn_Vr_mtx-1).*(tr_mtx-2*R_ref./(c.*D_fn_Vr_mtx)).^2);

% 下面将距离多普勒域的信号与变标方程相乘，实现“补余RCMC”
S_RD_1 = S_RD.*s_sc;            % 相位相乘，实现“补余RCMC”

disp(' 距离多普勒域，完成“补余RCMC” ');
%{
% 作图
figure;
imagesc(abs(S_RD));
title('原始数据变换到距离多普勒域，幅度');
figure;
imagesc(abs(S_RD_1));
title('距离多普勒域，补余RCMC后，幅度');
%}
clear S_RD;

%% 
% --------------------------------------------------------------------
% 变换到二维频域，进行“距离压缩，SRC，一致RCMC”
% --------------------------------------------------------------------
S_2df_1 = fft(S_RD_1,NFFT_r,2);         % 进行距离向FFT，变换到二维频域。距离零频在两端

% 完成距离压缩，SRC，一致RCMC这三者相位补偿的滤波器为：
H1 = exp(1j*pi.*D_fn_Vr_mtx./(D_fn_ref_Vr.*Km).*fr_mtx.^2)...
    .*exp(1j*4*pi/c.*(1./D_fn_Vr_mtx-1/D_fn_ref_Vr).*R_ref.*fr_mtx);
% 上面的H1距离零频在中心
W_ref = ones(Naz,1)*(kaiser(Nrg,3).');	% 距离向，构建Kaiser窗，此为矩阵形式，距离零频在中心。
% H1 = W_ref.*H1;             % 加入距离平滑窗，以抑制旁瓣，距离零频在中心。
% 下面通过fftshift将H1的距离零频调整到两端
H1 = fftshift(H1,2);        % 左右半边互换，距离零频在两端。

S_2df_2 = S_2df_1.*H1;    	% 在二维频域，相位相乘，实现距离压缩，SRC，一致RCMC

S_RD_2 = ifft(S_2df_2,NFFT_r,2);    % 进行距离IFFT，回到距离多普勒域，完成所有距离处理。

disp(' 在二维频域进行相位相乘，完成距离压缩，SRC，一致RCMC后，回到距离多普勒域 ');
%{
% 作图
figure;
imagesc(abs(S_2df_1));
title('变换到二维频域');
figure;
imagesc(abs(S_2df_2));
title('相位相乘，实现距离压缩，SRC，一致RCMC后，二维频域');
%
figure;
imagesc(abs(S_RD_2));
title('完成距离压缩，SRC，一致RCMC后，距离多普勒域');
%}
clear S_RD_1;
clear S_2df_1;
clear H1;
clear S_2df_2;

%%
% --------------------------------------------------------------------
% 距离多普勒域，完成“方位压缩”和“附加相位校正”
% --------------------------------------------------------------------
R0_RCMC = (c/2).*tr;   % 随距离线变化的R0，记为R0_RCMC，用来计算方位MF。

% 生成方位向匹配滤波器
Haz = exp(1j*4*pi.*(D_fn_Vr*R0_RCMC).*f0./c);       % 方位MF

% 附加相位校正项
H2 = exp(-1j*4*pi.*Km./(c^2).*(1-D_fn_Vr_mtx./D_fn_ref_Vr)...
    .*((1./D_fn_Vr)*R0_RCMC-R_ref./D_fn_Vr_mtx).^2); 	% 附加相位校正项

% 下面进行相位相乘，在距离多普勒域，同时完成方位MF和附加相位校正
S_RD_3 = S_RD_2.*Haz.*H2;           % 距离多普勒域，相位相乘

% 最后通过IFFT回到图像域，完成方未处理
s_image = ifft(S_RD_3,NFFT_a,1); 	% 完成成像过程，得到成像结果为：s_image

disp(' 完成“方位压缩”和“附加相位校正” ');
disp(' 成像结束 ');
%{
% 作图
figure;
imagesc(abs(S_RD_3));
title('距离多普勒域，进行了相位相乘后（方位MF和附加相位校正）');
%}
clear S_RD_2;
clear Haz;
clear H2;
clear S_RD_3;

%% 
% 下面对亮度进行非线性变换，减小对比度
sout = abs(s_image)/max(max(abs(s_image)));
G = 20*log10(sout+eps);             % dB显示
clim = [-55 0];                     % 动态显示范围
%{
figure;
imagesc(((0:Nrg-1)+first_rg_cell)/Fr*c/2+R0,((0:Naz-1)+first_rg_line)/Fa*Vr,G,clim);
axis xy;
title('RADARSAT-1数据，使用CS算法，成像结果')
xlabel('Range(m)')
ylabel('Azimuth(m)')
% colormap(gray);
%}

% 将图像向左移位：
%   基于CSA算法的成像位置是压至参考频率对应的距离单元，而非压至零多普勒处
%   得到的图像结果相比于压至零多普勒，是向右偏移的
% 因此进行以下向左移位
% 此外，还要进行上下半边互换
% 经过以上操作后，得到结果：
tmp = round(2*(R0/D_fn_ref_Vr-R0)/c*Fr);
s_tmp(:,1:Nrg-tmp+1) = G(:,tmp:end);
s_tmp(:,Nrg-tmp+1+1:Nrg) = G(:,1:tmp-1);
figure;
imagesc(((0:Nrg-1)+first_rg_cell)/Fr*c/2+R0,((0:Naz-1)+first_rg_line)/Fa*Vr,fftshift(s_tmp,1),clim);
axis xy;
title('RADARSAT-1数据，使用CS算法，成像结果')
xlabel('Range(m)')
ylabel('Azimuth(m)')

if b ==3
    % 对两个分区一起成像时，使用这部分来成像。
    % 作用是：将上下部分进行一定的移位
    %       （ 原来的图像的第2900行到最后一行应该在新图像的最开头 ）
    ss_tmp(1:Naz-2900+1,:) = s_tmp(2900:Naz,:);
    ss_tmp(Naz-2900+1+1:Naz,:) = s_tmp(1:2900-1,:);
    figure;
    imagesc(((0:Nrg-1)+first_rg_cell)/Fr*c/2+R0,((0:Naz-1)+first_rg_line)/Fa*Vr,ss_tmp,clim);
    axis xy;
    title('RADARSAT-1数据，使用CS算法，成像结果')
    xlabel('Range(m)')
    ylabel('Azimuth(m)')
end


