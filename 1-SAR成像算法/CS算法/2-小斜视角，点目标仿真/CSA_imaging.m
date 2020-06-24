%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            小斜视角（3.5°）
%                CSA
%              点目标仿真
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 该程序截止到 2014.10.29. 15:50 p.m.
%
% 2014.10.29. 修改：初步解决了“方位向定位问题”，即绝对位置已经能计算得出。

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
Nrg = 320;              % 距离线采样点数（即数据矩阵，列数）
sita_r_c = (3.5*pi)/180;	% 波束斜视角，3.5 度，这里转换为弧度
c = 3e8;                % 光速

R0 = R_nc*cos(sita_r_c);	% 与R_nc相对应的最近斜距，记为R0
Nr = Tr*Fr;             % 线性调频信号采样点数
BW_range = Kr*Tr;       % 距离向带宽
lamda = c/f0;           % 波长
fnc = 2*Vr*sin(sita_r_c)/lamda;     % 多普勒中心频率，根据公式（4.33）计算。
La_real = 0.886*2*Vr*cos(sita_r_c)/BW_dop;  % 方位向天线长度，根据公式（4.36）
beta_bw = 0.886*lamda/La_real;              % 雷达3dB波束
La = beta_bw*R0;        % 合成孔径长度
a_sr = Fr / BW_range;   % 距离向过采样因子
a_sa = Fa / BW_dop;     % 方位向过采样因子

Mamb = round(fnc/Fa);   % 多普勒模糊

NFFT_r = Nrg;           % 距离向FFT长度
NFFT_a = Naz;           % 方位向FFT长度

R_ref = R0;             % 参考目标选在场景中心，其最近斜距为 R_ref
fn_ref = fnc;        	% 参考目标的多普勒中心频率

% --------------------------------------------------------------------
% 设定仿真点目标的位置
% 以距离向作为x轴正方向
% 以方位向作为y轴正方向
% -------------------------------------------------------------------- 
delta_R0 = 0;       % 将目标1的波束中心穿越时刻，定义为方位向时间零点。
delta_R1 = 120; 	% 目标1和目标2的方位向距离差，120m
delta_R2 = 80;      % 目标2和目标3的距离向距离差，80m

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
fa = fnc + fftshift( -NFFT_a/2 : NFFT_a/2-1 )*( Fa/NFFT_a );	% 方位频率轴

% 生成距离（方位）时间（频率）矩阵
tr_mtx = ones(Naz,1)*tr;    % 距离时间轴矩阵，大小：Naz*Nrg
ta_mtx = ta.'*ones(1,Nrg);  % 方位时间轴矩阵，大小：Naz*Nrg
fr_mtx = ones(Naz,1)*fr;    % 距离频率轴矩阵，大小：Naz*Nrg
fa_mtx = fa.'*ones(1,Nrg);  % 方位频率轴矩阵，大小：Naz*Nrg

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
    % w_azimuth1是天线双程方向图。
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
text(300,-70,'图1，原始数据');       % 给图1进行文字说明 
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

% 作图
figure;
imagesc(abs(S_RD));
title('原始数据变换到距离多普勒域，幅度');
figure;
imagesc(abs(S_RD_1));
title('距离多普勒域，补余RCMC后，幅度');

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

% 作图
figure;
imagesc(abs(S_2df_1));
title('变换到二维频域');
figure;
imagesc(abs(S_2df_2));
title('相位相乘，实现距离压缩，SRC，一致RCMC后，二维频域');

figure;
imagesc(abs(S_RD_2));
title('完成距离压缩，SRC，一致RCMC后，距离多普勒域');

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

% 作图
figure;
imagesc(abs(S_RD_3));
title('距离多普勒域，进行了相位相乘后（方位MF和附加相位校正）');

figure;
imagesc(abs(s_image));
title('成像结果');

%%
% 下面通过调用函数，得到三个点目标各自的切片，并进行升采样
% 同时对点目标中心做距离向切片，方位向切片
% 计算出相应的指标：PSLR，ISLR，IRW
NN = 20;
% 分别得到每个点目标的切片放大；行切片、列切片；和相应的指标

% 目标1，点目标中心在 （ tg_1_x，tg_1_y ）
% =========================================================================
% 现在的点目标位置计算如下：
tg_1_x = rem( R0*tan(sita_r_c)/Vr*Fa , Naz );
if tg_1_x < Naz/2
    tg_1_x = tg_1_x + (Naz/2+1);
else
    tg_1_x = tg_1_x - (Naz/2+1);
end
tg_1_x = round(tg_1_x);    	% 四舍五入，得到整数值，作为点目标的方位中心坐标。
% 这里得到的 tg_1_x 即是点目标中心方位向的位置（坐标）。
% =========================================================================
% 下面计算目标1的距离向位置:
% 由于CSA的变标作用，从原来的压至零多普勒（R0），变为压至方位参考频率（fn_ref）处
% 的距离单元（即 R0/D_fn_ref_Vr ），因此对应的目标1的y轴位置如下，为 tg_1_y ：
tg_1_y = round( (Nrg/2+1) + 2*(R0/D_fn_ref_Vr-R0)/c*Fr );
target_1 = target_analysis_2( s_image(tg_1_x-NN:tg_1_x+NN,tg_1_y-NN:tg_1_y+NN),Fr,Fa,Vr);


% 目标2，点目标中心在 （tg_2_x，target_2_y）
tg_2_x = tg_1_x + delta_R1/Vr*Fa;
tg_2_y = tg_1_y;
% target_2 = target_analysis_2( s_image(tg_2_x-NN:tg_2_x+NN,tg_2_y-NN:tg_2_y+NN),Fr,Fa,Vr);


% 目标3，点目标中心在（tg_3_x，tg_3_y）
tg_3_x = tg_2_x + delta_R2*tan(sita_r_c)/Vr*Fa;
tg_3_x = fix(tg_3_x);
tg_3_y = tg_2_y + 2*(delta_R2/D_fn_ref_Vr)/c*Fr;
% target_3 = target_analysis_2( s_image(tg_3_x-NN:tg_3_x+NN,tg_3_y-NN:tg_3_y+NN),Fr,Fa,Vr);



