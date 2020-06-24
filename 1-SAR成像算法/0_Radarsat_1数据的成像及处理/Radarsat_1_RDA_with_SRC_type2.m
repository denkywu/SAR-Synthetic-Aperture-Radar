%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                       Radarsat_1 光盘中数据
%                             RDA 成像
%                     
%                           采用方式2实现SRC
%                    采用大斜视角下的程序进行成像
%
%                               WD
%                       2014.10.10. 13:49 p.m.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 程序说明：
% 主程序是：  Radarsat_1_RDA_with_SRC_type2.m
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
% 由采用 方式2 实现SRC 的大斜视角下的点目标程序修改而来；
% （4）成像流程：
%   ——原始数据
%   ——变换到距离频域、方位时域，进行距离压缩
%   ——回到二维时域，完成距离压缩
%   ——变换到二维频域，采用方式2，实现SRC
%           （在二维频域，相位相乘来实现SRC）
%   ——变换到距离多普勒域，实现RCMC
%   ——方位压缩
%   ——回到图像域，成像结束。
%
% 本程序修改截止到： 2014.10.10. 13:49 p.m.

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
b = 3;              % 选择对于哪一部分成像
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
    s_echo = [s_echo2;s_echo1]; % 将分区1和分区2的数据合成整个数据块，用于成像
end
clear data_1;clear data_2;clear s_echo1;clear s_echo2;

% 作图显示
figure;
imagesc(abs(s_echo));
title('原始数据');                          % 原始回波数据（未处理）的幅度图像
% colormap(gray);

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

%%
% --------------------------------------------------------------------
% 距离（方位）向时间，频率相关定义
% --------------------------------------------------------------------
% 距离
tr = T_start + ( -Nrg/2 : (Nrg/2-1) )/Fr;                       % 距离时间轴
fr = ( -NFFT_r/2 : NFFT_r/2-1 )*( Fr/NFFT_r );                  % 距离频率轴
% 方位
ta = ( -Naz/2: Naz/2-1 )/Fa;                                    % 方位时间轴
fa = fnc + fftshift( -NFFT_a/2 : NFFT_a/2-1 )*( Fa/NFFT_a );	% 方位频率轴

%%  
% --------------------------------------------------------------------
% 距离压缩
% --------------------------------------------------------------------
S_range = fft(s_echo,NFFT_r,2);     % 进行距离向傅里叶变换，零频在两端。

%　生成距离向匹配滤波器
% ====================================================
% 采用方式2
% 时域复制脉冲，末端补零，fft，再取复共轭。
t_ref = ( -Nr/2 : (Nr/2-1) )/Fr;    % 用来生成距离MF的距离时间轴
t_ref_mtx = ones(Naz,1)*t_ref;      % 矩阵形式
w_ref = kaiser(Nr,2.5);             % 距离向，构建Kaiser窗，此为列向量。
w_ref = ones(Naz,1)*(w_ref.');      % 构成矩阵形式，每一行都相同的加窗。

s_ref = exp((1j*pi*Kr).*((t_ref_mtx).^2)); % 复制（发射）脉冲，未加窗。
% s_ref = w_ref.*exp((1j*pi*Kr).*((t_ref_mtx).^2)); % 复制（发射）脉冲，加了窗。

s_ref = [s_ref,zeros(Naz,Nrg-Nr)];      % 对复制脉冲，后端补零。
 
S_ref = fft(s_ref,NFFT_r,2);            % 复制脉冲的距离傅里叶变换，零频在两端。
H_range = conj(S_ref);                  % 距离向匹配滤波器，零频在两端。
S_range_c = S_range.*H_range;           % 乘以匹配滤波器，零频在两端。    
s_rc = ifft(S_range_c,[],2);            % 完成距离压缩，回到二维时域。
% s_rc的长度为：Naz*Nrg。未去除弃置区。

disp('完成距离向匹配滤波，得到处理后的信号，是 s_rc')
% 作图显示
figure;
imagesc(abs(s_rc));
title('距离压缩后，时域');              % 距离压缩后的数据，幅度图像 
% colormap(gray);

%%
% --------------------------------------------------------------------
% 变换到二维频域，进行SRC
% --------------------------------------------------------------------
s_rc = s_rc.*exp(-1j*2*pi*fnc.*(ta.'*ones(1,Nrg)));    % 数据搬移
S_2df = fft(s_rc,NFFT_a,1);         	% 方位向傅里叶变换，到距离多普勒域
S_2df = fft(S_2df,Nrg,2);               % 距离向傅里叶变换，到二维频域
% ！！！注意：距离向零频在两端。
% ====================================================================
% 设置方位频率轴——这是关键点
fa = fnc + fftshift(-NFFT_a/2:NFFT_a/2-1)/NFFT_a*Fa; 	% 方位频率轴如此设置。
% =====================================================================
D_fn_Vr = sqrt(1-lamda^2.*(fa.').^2./(4*Vr^2));         % 大斜视角下的徙动因子
K_src = 2*Vr^2*f0^3.*D_fn_Vr.^3./(c*R0*(fa.').^2);      % 列向量
K_src_1 = 1./K_src;             % 列向量。为了后面能使用矩阵乘法，这里先求倒数
H_src = exp(-1j*pi.*K_src_1*(fr.^2)); % 二次距离压缩滤波器。距离向，零频在中间。
% 这是矩阵，大小Naz*Nrg
H_src = fftshift(H_src,2);      % （左右半边互换）距离向，零频在两端。 ！！！这很关键！！！

S_2df_src = S_2df.*H_src;       % 这一步点乘时，要注意两者的距离向频率轴应该对应上，不然会出错！！
% 这就是为什么上面的 H_src 要 fftshift 的原因！！

S_rd = ifft(S_2df_src,[],2);   	% 完成二次距离压缩（SRC），回到距离多普勒域。

disp('在二维频域，采取方式2完成SRC，并变换到距离多普勒域');
% 作图
figure;
imagesc(abs(S_rd));
title('完成距离压缩，完成SRC后，距离多普勒域（未RCMC）');

%%
% --------------------------------------------------------------------
% 变换到距离多普勒域，进行距离徙动校正
% --------------------------------------------------------------------
% 每一个最近斜距（R0）都随着距离门的不同而改变。
tr_RCMC = T_start + ( -Nrg/2 : (Nrg/2-1) )/Fr;   % 用于RCMC和方位MF的距离时间轴。
R0_RCMC = (c/2).*tr_RCMC;   % 随距离线变化的R0，记为R0_RCMC，用来计算RCM和Ka。
delta_Rrd_fn = ((1-D_fn_Vr)./D_fn_Vr)*R0_RCMC;      % 大斜视角下的RCM

num_range = c/(2*Fr);   % 一个距离采样单元，对应的长度
delta_Rrd_fn_num = delta_Rrd_fn./num_range; % 每一个方位向频率，其RCM对应的距离采样单元数

R = 8;  % sinc插值核长度
S_rd_rcmc = zeros(NFFT_a,Nrg); % 用来存放RCMC后的值
for p = 1 : NFFT_a
    for q = 1 : Nrg   % 此时距离向的长度是 (Nrg-Nr+1)=Nrg        
        delta_Rrd_fn_p = delta_Rrd_fn_num(p,q);
        Rrd_fn_p = q + delta_Rrd_fn_p;
        
        Rrd_fn_p = rem(Rrd_fn_p,Nrg);  % 由于RCM的长度会超过Nrg，所以这样处理一下。
        
        Rrd_fn_p_zheng = ceil(Rrd_fn_p);        % ceil，向上取整。
        ii = ( Rrd_fn_p-(Rrd_fn_p_zheng-R/2):-1:Rrd_fn_p-(Rrd_fn_p_zheng+R/2-1)  );        
        rcmc_sinc = sinc(ii);
        rcmc_sinc = rcmc_sinc/sum(rcmc_sinc);   % 插值核的归一化
        % ii 是sinc插值过程的变量;
        % g(x)=sum(h(ii)*g_d(x-ii)) = sum(h(ii)*g_d(ll));
               
        % 由于S_rd只有整数点取值，且范围有限。因此插值中要考虑它的取值溢出边界问题。
        % 这里我采取循环移位的思想，用来解决取值溢出问题。
        if (Rrd_fn_p_zheng-R/2) > Nrg    % 全右溢
            ll = (Rrd_fn_p_zheng-R/2-Nrg:1:Rrd_fn_p_zheng+R/2-1-Nrg);
        else
            if (Rrd_fn_p_zheng+R/2-1) > Nrg    % 部分右溢
                ll_1 = (Rrd_fn_p_zheng-R/2:1:Nrg);
                ll_2 = (1:1:Rrd_fn_p_zheng+R/2-1-Nrg);
                ll = [ll_1,ll_2];
            else
                if (Rrd_fn_p_zheng+R/2-1) < 1    % 全左溢（不可能发生，但还是要考虑）
                    ll = (Rrd_fn_p_zheng-R/2+Nrg:1:Rrd_fn_p_zheng+R/2-1+Nrg);
                else
                    if (Rrd_fn_p_zheng-R/2) < 1       % 部分左溢
                        ll_1 = (Rrd_fn_p_zheng-R/2+Nrg:1:Nrg);
                        ll_2 = (1:1:Rrd_fn_p_zheng+R/2-1);
                        ll = [ll_1,ll_2];
                    else
                        ll = (Rrd_fn_p_zheng-R/2:1:Rrd_fn_p_zheng+R/2-1);
                    end                    
                end
            end
        end   
        rcmc_S_rd = S_rd(p,ll);
        S_rd_rcmc(p,q) = sum( rcmc_sinc.*rcmc_S_rd );
    end
end
% S_rd_rcmc 就是RCMC后的距离多普勒域频谱。

disp('完成RCMC后，距离多普勒域');
% 作图，显示，完成RCMC后，在距离多普勒域的图像。
figure;
imagesc(abs(S_rd_rcmc));
title('完成RCMC后，距离多普勒域');
% colormap(gray);

%%
% --------------------------------------------------------------------
% 方位压缩
% --------------------------------------------------------------------
fa_azimuth_MF = fa;         % 方位频率轴，采用和RCMC中所用的频率轴相同。
Haz = exp(1j*4*pi.*(D_fn_Vr*R0_RCMC).*f0./c);   % 改进的方位向MF

S_rd_c = S_rd_rcmc.*Haz;            % 乘以匹配滤波器
s_ac = ifft(S_rd_c,[],1);       	% 完成方位压缩，变到图像域。结束。

disp('完成方位压缩后，回到时域');
% 至此，全部的信号处理过程已经完成。剩下的就是进行成像。

%% 
% --------------------------------------------------------------------
% 绘制SAR图像
% --------------------------------------------------------------------
% 直接显示图像
% figure;
% imagesc(abs(s_ac));
% title('图像，时域');
% colormap(gray);

% 设置动态显示范围后，再显示图像。
sout = (s_ac);
sout = abs(sout)/max(max(abs(sout)));
G = 20*log10(sout+eps);                         % db显示
clim = [-55,0];                                 % 动态显示范围
figure;
imagesc(((0:Nrg-1)+first_rg_cell)/Fr*c/2+R0,((0:Naz-1)+first_rg_line)/Fa*Vr,G,clim);
axis xy;
title('RADARSAT-1数据，使用RD算法，成像结果')
xlabel('Range(m)')
ylabel('Azimuth(m)')
% colormap(gray);






