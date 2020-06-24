%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       RDA成像——点目标仿真
%               大斜视角
%               SRC方式3
%  在距离频域将SRC与距离压缩滤波器合并
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 截止到 2014.10.10. 17:36 p.m.

% 2014.10.09. 修改：计算IRW指标时，采用修改后的函数：zhibiao_2( )
% 2014.10.09. 修改：对取出的成像切片进行二维升采样时，原来是先逐行判断最小值处，
%                   然后补零；再逐列判断最小值处，并补零。
%                       ——这会带来一些问题
%                   于是采取的方法是：修改上面的先后顺序，先对每列补零，再对
%                   每行补零。
%                       ——这样修改后的二维升采样，效果便很理想。
%                       ——对函数 target_analysis_2( ) 进行修改
%
% 2014.10.10. 修改：对天线双程方向图的修改。我认为，其方向图的长度不应该任由
%                   它随着Naz的大小而变，也应该用一个矩形窗（比如合成孔径长度）
%                   加以限制。我原来试过1个合成孔径长度，但是这样的结果，其PSRL
%                   稍大概在17,18dB。因此，我现在用书上图5.16（137页）进行计算，
%                   取 1.135 个合成孔径长度。
%                       ——也就是说，生成天线双程方向图的方法不变。但是加上一个
%                           1.135个合成孔径长度的矩形窗，对数据进行限制。
%
% 还存在的问题：
%   1. 计算IRW之前，由于使用了函数 imrotate 对切片进行旋转，使得距离向或方位向
%      转到水平轴或者垂直轴。这过程当中可能会涉及到插值，所以直接利用点数计算
%      IRW 不正确。
%           —— 由于不知道函数 imrotate 具体怎么操作的，所以在这种情况下计算
%                IRW 不好办。最好的方法是自己写一个转角的程序，但是工作量比较大
%                我现在暂时不考虑这个问题了。
%   2. 绝位置的计算
%           ——始终没有解决

%%
clear;
close all;
clc;
% --------------------------------------------------------------------
% 定义参数
% --------------------------------------------------------------------
R_nc = 20e3;            % 景中心斜距
Vr = 150;               % 雷达有效速度
Tr = 2.5e-6;            % 发射脉冲时宽
Kr = 20e12;             % 距离调频率
f0 = 5.3e9;             % 雷达工作频率
BW_dop = 80;            % 多普勒带宽
Fr = 60e6;              % 距离采样率
Fa = 200;               % 方位采样率
Naz = 1024;          	% 距离线数（即数据矩阵，行数）——这里修改为1024。
Nrg = 5*320;               	% 距离线采样点数（即数据矩阵，列数）
%    ——这里的 Nrg 设计的足够大，使得原始数据能够被完整包括。
sita_r_c = (21.9*pi)/180;	% 波束斜视角，21.9 度，这里转换为弧度
c = 3e8;                % 光速

R0 = R_nc*cos(sita_r_c);	% 与R_nc相对应的最近斜距，记为R0
Nr = Tr*Fr;             % 线性调频信号采样点数
BW_range = Kr*Tr;       % 距离向带宽
lamda = c/f0;           % 波长
fnc = 2*Vr*sin(sita_r_c)/lamda;     % 多普勒中心频率，根据公式（4.33）计算。
La_real = 0.886*2*Vr*cos(sita_r_c)/BW_dop; % 方位向天线长度，根据公式（4.36）
beta_bw = 0.886*lamda/La_real;           % 雷达3dB波束
La = 0.886*R_nc*lamda/La_real;              % 合成孔径长度
a_sr = Fr / BW_range;   % 距离向过采样因子
a_sa = Fa / BW_dop;     % 方位向过采样因子

Mamb = round(fnc/Fa);   % 多普勒模糊

NFFT_r = Nrg;           % 距离向FFT长度
NFFT_a = Naz;           % 方位向FFT长度

% --------------------------------------------------------------------
% 设定仿真点目标的位置
% 以距离向作为x轴正方向
% 以方位向作为y轴正方向
% -------------------------------------------------------------------- 
delta_R0 = 0;       % 将目标1的波束中心穿越时刻，定义为方位向时间零点。
delta_R1 = 120; 	% 目标1和目标2的方位向距离差，120m
delta_R2 = 50;      % 目标2和目标3的距离向距离差，50m

% 目标1
x1 = R0;            % 目标1的距离向距离
y1 = delta_R0 + x1*tan(sita_r_c);	% 目标1的方位向距离

% 目标2
x2 = x1;            % 目标2和目标1的距离向距离相同
y2 = y1 + delta_R1; % 目标2的方位向距离
% 目标3
x3 = x2 + delta_R2;                 % 目标3和目标2有距离向的距离差，为delta_R2
y3 = y2 + delta_R2*tan(sita_r_c);  	% 目标3的方位向距离
% 定义以下数组，便于处理
x_range = [x1,x2,x3];
y_azimuth = [y1,y2,y3];

% 计算三个目标各自的波束中心穿越时刻
nc_1 = (y1-x1*tan(sita_r_c))/Vr;    % 目标1的波束中心穿越时刻。
nc_2 = (y2-x2*tan(sita_r_c))/Vr;    % 目标2的波束中心穿越时刻。
nc_3 = (y3-x3*tan(sita_r_c))/Vr;    % 目标3的波束中心穿越时刻。
nc_target = [nc_1,nc_2,nc_3];       % 定义该数组，便于处理。

%% 
% --------------------------------------------------------------------
% 距离（方位）向时间，频率相关定义
% --------------------------------------------------------------------
% 距离
tr = 2*R0/c + ( -Nrg/2 : (Nrg/2-1) )/Fr;                % 距离时间轴
fr = ( -NFFT_r/2 : NFFT_r/2-1 )*( Fr/NFFT_r );          % 距离频率轴
% 方位
ta = ( -Naz/2: Naz/2-1 )/Fa;                            % 方位时间轴
fa = fnc + ( -NFFT_a/2 : NFFT_a/2-1 )*( Fa/NFFT_a );	% 方位频率轴

% 生成距离（方位）时间（频率）矩阵
tr_mtx = ones(Naz,1)*tr;    % 距离时间轴矩阵，大小：Naz*Nrg
ta_mtx = ta.'*ones(1,Nrg);  % 方位时间轴矩阵，大小：Naz*Nrg

%% 
% --------------------------------------------------------------------
% 生成点目标原始数据
% --------------------------------------------------------------------
s_echo = zeros(Naz,Nrg);    % 用来存放生成的回波数据

A0 = 1;                     % 目标回波幅度，都设置为1.
for k = 1:1                 % 生成k个目标的原始回波数据
    R_n = sqrt( (x_range(k).*ones(Naz,Nrg)).^2 + (Vr.*ta_mtx-y_azimuth(k).*ones(Naz,Nrg)).^2 );% 目标k的瞬时斜距
    w_range = ((abs(tr_mtx-2.*R_n./c)) <= ((Tr/2).*ones(Naz,Nrg)));     % 距离向包络，即距离窗
    % =====================================================================    
    % 方位向包络，也就是 天线的双程方向图作用因子。
    %{
    % 方式1
    % sinc平方型函数，根据公式（4.31）计算    
    % 用每个目标对应的 波束中心穿越时刻 。
    sita = atan( Vr.*(ta_mtx-nc_target(k).*ones(Naz,Nrg))/x_range(k) );
    w_azimuth1 = (sinc(0.886.*sita./beta_bw)).^2; 
    % w_azimuth1是天线双程方向图。和原来一样，这里没有修改。
    
    % 下面的 w_azimuth2 是和方式2的矩形窗相同的构造方法，目的是：对天线双程
    % 方向图进行数据限制：限制为 1.135 个合成孔径长度。 
    w_azimuth2 = (abs(ta - nc_target(k)) <= (1.135*La/2)/Vr);    
    w_azimuth2 = w_azimuth2.'*ones(1,Nrg);	% 用来对 w_azimuth1 的天线双程方向图作数据限制。
    
    % 下面将两者相乘，得到仿真中所用的天线加权
    w_azimuth = w_azimuth1.*w_azimuth2;     % 两者相乘，得到仿真中所用的天线加权
    clear w_azimuth1;
    clear w_azimuth2;
    %}
    %
    % 方式2
    % 利用合成孔径长度，直接构造矩形窗（其实这里只是限制数据范围，没有真正加窗）
    w_azimuth = (abs(ta - nc_target(k)) <= (La/2)/Vr);    % 行向量
    w_azimuth = w_azimuth.'*ones(1,Nrg);    % 生成Naz*Nrg的矩阵
    %}
    % =====================================================================   
    s_k = A0.*w_range.*w_azimuth.*exp(-(1j*4*pi*f0).*R_n./c).*exp((1j*pi*Kr).*(tr_mtx-2.*R_n./c).^2);
    % 上式就是生成的某一个点目标（目标k）的回波信号。
    % 经过几次循环，生成几个点目标的回波信号，相加即可。
    if k == 1
        s_1 = s_k;          % 目标1的回波信号
    end
    if k == 2   
        s_2 = s_k;          % 目标2的回波信号
    end
    if k == 3
        s_3 = s_k;          % 目标3的回波信号
    end
    s_echo = s_echo + s_k;  % 所有点目标回波信号之和   
end
% s_echo 就是我们需要的原始数据，点目标回波信号。

% 作图
% 图1——原始数据
figure;
subplot(2,2,1);
imagesc(real(s_echo));
title('（a）实部');
xlabel('距离时域（采样点）');
ylabel('方位时域（采样点）');
text(1500,-60,'图1，原始数据');       % 给图1进行文字说明 
% text 函数：在图像的指定坐标位置，添加文本框

subplot(2,2,2);
imagesc(imag(s_echo));
title('（b）虚部');
xlabel('距离时域（采样点）');
ylabel('方位时域（采样点）');

subplot(2,2,3);
imagesc(abs(s_echo));
title('（c）幅度');
xlabel('距离时域（采样点）');
ylabel('方位时域（采样点）');

subplot(2,2,4);
imagesc(angle(s_echo));
title('（d）相位');
xlabel('距离时域（采样点）');
ylabel('方位时域（采样点）');
% colormap(gray);

figure;
subplot(2,2,1);
imagesc(abs(fft(s_echo,[],1)));
title('RD 频谱幅度');
subplot(2,2,2);
imagesc(angle(fft(s_echo,[],1)));
title('RD 频谱相位');
subplot(2,2,3);
imagesc(abs(fft2(s_echo)));
title('二维频谱幅度');
subplot(2,2,4);
imagesc(angle(fft2(s_echo)));
title('二维频谱相位');
% colormap(gray);

%%
% --------------------------------------------------------------------
% 距离压缩，同时进行SRC
% --------------------------------------------------------------------
S_range = fft(s_echo,NFFT_r,2);     % 进行距离向傅里叶变换，零频在两端。

%
% 作图
% 图2——距离频域，方位时域，频谱（未距离压缩）
figure;
subplot(1,2,1);
imagesc(real(S_range));
title('（a）实部');
xlabel('距离频域（采样点）');
ylabel('方位时域（采样点）');
text(1500,-60,'图2，距离频域');       % 给图2进行文字说明
text(1700,-20,'未压缩');       

subplot(1,2,2);
imagesc(abs(S_range));
title('（b）幅度');
xlabel('距离频域（采样点）');
ylabel('方位时域（采样点）');
%}

%　生成距离向匹配滤波器
% ====================================================
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

% ------------------------------------------------------------------------
% 在距离频域将 SRC 与 距离压缩 合并处理
% 构造SRC滤波器，如下：
% 在方式3中，SRC滤波器其中的 R0 和 fa 都要取为定值
% R0 取为目标A和B的最近斜距
% fa 取为多普勒中心频率 fnc
D_fn_Vr = sqrt(1-lamda^2*fnc^2/(4*Vr^2));   % 大斜视角下的徙动因子。
% SRC方式3工作在距离频域。此时，忽略方位频率的相关性，将fa取为定值，即采用fnc。
% 要注意的是：后面进行RCMC时，我会重新计算这个徙动因子，那里会随着方位频率而变化。

K_src = 2*Vr^2*f0^3*D_fn_Vr^3/(c*R0*fnc^2);  % 定值
H_src = exp(-1j*pi*(fr.^2)/K_src);   % 二次距离压缩滤波器。距离向，零频在中间！！
H_src = ones(Naz,1)*H_src;
% 这是矩阵，大小 Naz*Nrg
H_src = fftshift(H_src,2);  % （左右半边互换）距离向，零频在两端。 ！！！这很关键！！！
% ------------------------------------------------------------------------
S_range_c = S_range.*H_range.*H_src;   	% 乘以匹配滤波器和SRC滤波器，零频在两端。  

s_rc = ifft(S_range_c,[],2);            % 完成距离压缩，回到二维时域。
% s_rc的长度为：Naz*Nrg。未去除弃置区。

% 对s_rc进行去除弃置区的操作
% 弃置区长度为：2*（Nr-1）
% 完全卷积的结果（我们需要的）长度：（Nrg-Nr+1），记为 N_rg。
N_rg = Nrg-Nr+1;                        % 完全卷积的长度
s_rc_c = zeros(Naz,N_rg);               % 用来存放去除弃置区后的数据
s_rc_c = s_rc(:,1:N_rg);                % 取前 N_rg列。
% ====================================================

%
% 作图
% 图3——距离频域，方位时域，频谱（已距离压缩，并同时完成了SRC（采用方式3））
figure;
subplot(1,2,1);
imagesc(real(S_range_c));
title('（a）实部');
xlabel('距离频域（采样点）');
ylabel('方位时域（采样点）');
text(1500,-60,'图3，距离频域');       % 给图3进行文字说明
text(1400,-20,'已压缩，并完成了SRC');       

subplot(1,2,2);
imagesc(abs(S_range_c));
title('（b）幅度');
xlabel('距离频域（采样点）');
ylabel('方位时域（采样点）');
%}
%
% 作图
% 图4——二维时域（完成距离压缩，完成SRC，变回到二维时域）
figure;
subplot(1,2,1);
imagesc(real(s_rc_c));  %　这及其以下，都直接使用去除弃置区后的结果
title('（a）实部');
xlabel('距离时域（采样点）');
ylabel('方位时域（采样点）');
text(1350,-60,'图4，二维时域');       % 给图4进行文字说明
text(1350,-20,'完成压缩，完成SRC');       

subplot(1,2,2);
imagesc(abs(s_rc_c));
title('（b）幅度');
xlabel('距离时域（采样点）');
ylabel('方位时域（采样点）');
%}

%%
% --------------------------------------------------------------------
% 距离多普勒域，进行距离徙动校正
% --------------------------------------------------------------------
s_rc_c = s_rc_c.*exp(-1j*2*pi*fnc.*(ta.'*ones(1,N_rg)));    % 数据搬移
S_rd = fft(s_rc_c,NFFT_a,1);        % 方位向傅里叶变换，到距离多普勒域。

% 作图
% 图5——距离多普勒域（未RCMC）
figure;
subplot(1,2,1);
imagesc(real(S_rd));
title('（a）实部');
xlabel('距离时域（采样点）');
ylabel('方位频域（采样点）');
text(1350,-60,'图5，距离多普勒域');       % 给图5进行文字说明
text(1550,-20,'未RCMC');       

subplot(1,2,2);
imagesc(abs(S_rd));
title('（b）幅度');
xlabel('距离时域（采样点）');
ylabel('方位频域（采样点）');

% ====================================================================
% 设置方位频率轴——关键点！！！
fa = fnc + fftshift(-NFFT_a/2:NFFT_a/2-1)/NFFT_a*Fa;    % 方位频率轴如此设置。
% =====================================================================

% 每一个最近斜距（R0）都随着距离门的不同而改变。
tr_RCMC = 2*R0/c + ( -N_rg/2 : (N_rg/2-1) )/Fr;	% 在新的距离线长度下的时间轴。
R0_RCMC = (c/2).*tr_RCMC;   % 随距离线变化的R0，记为R0_RCMC，用来计算RCM。

D_fn_Vr = sqrt(1-lamda^2.*(fa.').^2./(4*Vr^2)); % 大斜视角下的徙动因子，列向量。
delta_Rrd_fn = ((1-D_fn_Vr)./D_fn_Vr)*R0_RCMC;	% 大斜视角下的RCM。这是矩阵，大小Naz*N_rg。

num_range = c/(2*Fr);   % 一个距离采样单元，对应的长度
delta_Rrd_fn_num = delta_Rrd_fn./num_range; % 每一个方位向频率，其RCM对应的距离采样单元数

R = 8;  % sinc插值核长度
S_rd_rcmc = zeros(NFFT_a,N_rg); % 用来存放RCMC后的值
for p = 1 : NFFT_a
    for q = 1 : N_rg   % 此时距离向的长度是 (Nrg-Nr+1)=N_rg        
        delta_Rrd_fn_p = delta_Rrd_fn_num(p,q);
        Rrd_fn_p = q + delta_Rrd_fn_p;
        
        Rrd_fn_p = rem(Rrd_fn_p,N_rg);  % 由于RCM的长度会超过N_rg，所以这样处理一下（取余）。
        
        Rrd_fn_p_zheng = ceil(Rrd_fn_p);        % ceil，向上取整。
        ii = ( Rrd_fn_p-(Rrd_fn_p_zheng-R/2):-1:Rrd_fn_p-(Rrd_fn_p_zheng+R/2-1)  );        
        rcmc_sinc = sinc(ii);
        rcmc_sinc = rcmc_sinc/sum(rcmc_sinc);   % 插值核的归一化
        % ii 是sinc插值过程的变量;
        % g(x)=sum(h(ii)*g_d(x-ii)) = sum(h(ii)*g_d(ll));
               
        % 由于S_rd只有整数点取值，且范围有限。因此插值中要考虑它的取值溢出边界问题。
        % 这里我采取循环移位的思想，用来解决取值溢出问题。
        if (Rrd_fn_p_zheng-R/2) > N_rg    % 全右溢
            ll = (Rrd_fn_p_zheng-R/2-N_rg:1:Rrd_fn_p_zheng+R/2-1-N_rg);
        else
            if (Rrd_fn_p_zheng+R/2-1) > N_rg    % 部分右溢
                ll_1 = (Rrd_fn_p_zheng-R/2:1:N_rg);
                ll_2 = (1:1:Rrd_fn_p_zheng+R/2-1-N_rg);
                ll = [ll_1,ll_2];
            else
                if (Rrd_fn_p_zheng+R/2-1) < 1    % 全左溢（不可能发生，但还是要考虑）
                    ll = (Rrd_fn_p_zheng-R/2+N_rg:1:Rrd_fn_p_zheng+R/2-1+N_rg);
                else
                    if (Rrd_fn_p_zheng-R/2) < 1       % 部分左溢
                        ll_1 = (Rrd_fn_p_zheng-R/2+N_rg:1:N_rg);
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

% 作图
% 图6——距离多普勒域，RCMC后的结果
figure;
subplot(1,2,1);
imagesc(real(S_rd_rcmc));
title('（a）实部');
xlabel('距离时域（采样点）');
ylabel('方位频域（采样点）');
text(1350,-60,'图6，距离多普勒域');       % 给图6进行文字说明
text(1550,-20,'已RCMC');       

subplot(1,2,2);
imagesc(abs(S_rd_rcmc));
title('（b）幅度');
xlabel('距离时域（采样点）');
ylabel('方位频域（采样点）');

%%
% --------------------------------------------------------------------
% 方位压缩
% --------------------------------------------------------------------
fa_azimuth_MF = fa;         % 方位频率轴，采用和RCMC中所用的频率轴相同。
Haz = exp(1j*4*pi.*(D_fn_Vr*R0_RCMC).*f0./c);   % 大斜视角下，改进的方位向MF
% 这里要注意，生成的MF的零频既不是在两端，也不是在中心的。
% 考虑下频率轴是什么样的，间断点在哪里。注意fa的构成。
% 这里的频率轴和距离多普勒域的方位频谱是对应的。

S_rd_c = S_rd_rcmc.*Haz;            % 乘以匹配滤波器
s_ac = ifft(S_rd_c,[],1);       	% 完成方位压缩，变到图像域。结束。

% 作图
% 图7——成像结果
figure;
imagesc(abs(s_ac));
title('点目标成像');
xlabel('距离时域（采样点）');
ylabel('方位时域（采样点）');     

%%
% 下面通过调用函数，得到三个点目标各自的切片，并进行升采样
% 同时对点目标中心做距离向切片，方位向切片
% 计算出相应的指标：PSLR，ISLR，IRW

NN = 20;
% 分别得到每个点目标的切片放大；行切片、列切片；和相应的指标
% 目标1，点目标中心在 （ tg_1_x，tg_1_y ）
tg_1_x = 219;
tg_1_y = round(N_rg/2);
target_1 = target_analysis_2( s_ac(tg_1_x-NN:tg_1_x+NN,tg_1_y-NN:tg_1_y+NN),Fr,Fa,Vr);

% 目标2，点目标中心在 （tg_2_x，target_2_y）
tg_2_x = tg_1_x + delta_R1/Vr*Fa;
tg_2_y = tg_1_y;
% target_2 = target_analysis_2( s_ac(tg_2_x-NN:tg_2_x+NN,tg_2_y-NN:tg_2_y+NN),Fr,Fa,Vr);

% 目标3，点目标中心在（tg_3_x，tg_3_y）
tg_3_x = tg_2_x + delta_R2*tan(sita_r_c)/Vr*Fa;
tg_3_x = fix(tg_3_x);
tg_3_y = tg_2_y + 2*delta_R2/c*Fr;
% target_3 = target_analysis_2( s_ac(tg_3_x-NN:tg_3_x+NN,tg_3_y-NN:tg_3_y+NN),Fr,Fa,Vr);




