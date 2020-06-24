%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               正侧视
%                CSA
%              点目标仿真
%
%            仿真“运动补偿”
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 该程序截止到 2014.11.21. 16:23 p.m.
%
% 点目标仿真
% 仿真内容：“运动补偿”
%
% 2014.11.21. 16:23 p.m. 修改了原来程序中的一些错误

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
sita_r_c = (0*pi)/180;	% 波束斜视角，0 度，这里转换为弧度
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

theta = pi/3;           % 天线下视角
H = R0*cos(theta);      % 理想载机平台高度

% --------------------------------------------------------------------
% 设定仿真点目标的位置
% 以距离向作为x轴正方向
% 以方位向作为y轴正方向
% -------------------------------------------------------------------- 
delta_R0 = 0;       % 将目标1的波束中心穿越时刻，定义为方位向时间零点。
delta_R1 = 120; 	% 目标1和目标2的方位向距离差，120m（地距y轴）
delta_R2 = 80;      % 目标2和目标3的地距x轴坐标差，80m（地距x轴）

% 目标1
r1 = R0;                % 目标1的最近斜距
x1 = r1*sin(theta);     % 地距平面上的x轴坐标
z1 = 0;                 % z轴坐标（设置为在地平面上）
y1 = delta_R0 + r1*tan(sita_r_c);	% 目标1的方位向距离（地距平面的y轴坐标）

% 目标2（相对于目标1，只有地距y轴坐标不同，即方位向位置不同）
x2 = x1;                % 地距平面上的x轴坐标
z2 = 0;                 % z轴坐标（设置为在地平面上）
r2 = r1;                % 目标2的最近斜距（与目标1相同）
y2 = y1 + delta_R1;     % 目标2的方位向距离（地距平面的y轴坐标）

% 目标3（相对于目标2，只有地距x轴坐标不同）
x3 = x2 + delta_R2;    	% 地距平面的x轴坐标
z3 = 0;               	% z轴坐标（设置为在地平面上）
r3 = sqrt(x3^2 + H^2);  % 目标3的最近斜距
y3 = y2;                % 目标3的方位向距离（地距平面的y轴坐标）

% 目标4（与目标2和目标3，只有地距x轴坐标不同）
x4 = x2 - delta_R2;    	% 地距平面的x轴坐标
z4 = 0;               	% z轴坐标（设置为在地平面上）
r4 = sqrt(x4^2 + H^2);  % 目标4的最近斜距
y4 = y3;                % 目标4的方位向距离（地距平面的y轴坐标） 

% 定义以下数组，便于处理
target_x = [x1,x2,x3,x4];      % 目标的地距 x 轴坐标
target_z = [z1,z2,z3,z4];      % 目标的 z 轴坐标
r_range = [r1,r2,r3,r3];       % 目标的最近斜距
y_azimuth = [y1,y2,y3,y4];     % 目标的地距 y 轴坐标（方位向距离）

% 计算三个目标各自的波束中心穿越时刻
nc_1 = (y1-r1*tan(sita_r_c))/Vr;    % 目标1的波束中心穿越时刻。
nc_2 = (y2-r2*tan(sita_r_c))/Vr;    % 目标2的波束中心穿越时刻。
nc_3 = (y3-r3*tan(sita_r_c))/Vr;    % 目标3的波束中心穿越时刻。
nc_4 = (y4-r4*tan(sita_r_c))/Vr;    % 目标4的波束中心穿越时刻
nc_target = [nc_1,nc_2,nc_3,nc_4];       % 定义该数组，便于处理。

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
%
% ==================================================================
% 生成载机运动平台的运动误差
% 沿地距 x 轴的运动误差
a = 4;
w = 8;     % 这个控制着正弦误差是属于低频误差，还是高频误差
% w < 1，是低频误差
% w > 1，是高频误差
delta_x_t = a*sin(2*pi*w/La*Vr.*ta_mtx);    % 这是沿地距 x 轴的运动误差
% delta_x_t = 0;
% 沿z轴（载机平台高度）的运动误差
% delta_z_t = a*sin(2*pi*w/La*Vr.*ta_mtx);   	% 这是沿 z 轴的运动误差
delta_z_t = 0;
% ==================================================================

% 下面生成点目标原始数据
s_echo = zeros(Naz,Nrg);    % 用来存放生成的回波数据
A0 = 1;                     % 目标回波幅度，都设置为1.
for k = 1:4                 % 生成k个目标的原始回波数据
    R_n = sqrt( (delta_x_t-target_x(k)).^2 + (Vr.*ta_mtx-y_azimuth(k)).^2 ...
        + (H.*ones(Naz,Nrg)+delta_z_t).^2 );    % 目标k的瞬时斜距
%     R_n = sqrt( (r_range(k).*ones(Naz,Nrg)).^2 + (Vr.*ta_mtx-y_azimuth(k).*ones(Naz,Nrg)).^2 );% 无运动误差时，目标k的瞬时斜距
    w_range = ((abs(tr_mtx-2.*R_n./c)) <= ((Tr/2).*ones(Naz,Nrg)));     % 距离向包络，即距离窗
    % =====================================================================    
    % 方位向包络，也就是 天线的双程方向图作用因子。
    %{
    % 方式1
    % sinc平方型函数，根据公式（4.31）计算    
    % 用每个目标对应的 波束中心穿越时刻 。
    sita = atan( Vr.*(ta_mtx-nc_target(k).*ones(Naz,Nrg))/r_range(k) );
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
    %{
    if k == 1
        s_1 = s_k;          % 目标1的回波信号
    end
    if k == 2   
        s_2 = s_k;          % 目标2的回波信号
    end
    if k == 3
        s_3 = s_k;          % 目标3的回波信号
    end
    if k == 4
        s_4 = s_k;          % 目标4的回波信号
    end
    %}
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
% 计算在LOS方向的运动误差，这是进行“运动补偿”的依据
% --------------------------------------------------------------------
r = tr_mtx.*c/2;        % 用于运动误差计算的斜距r
% 沿LOS方向的运动误差
delta_r = delta_x_t.*( sqrt(r.^2-H^2)./r ) - delta_z_t.*(H./r);% 沿LOS方向，总的运动误差
delta_r = -delta_r*cos(sita_r_c);
delta_r_R0 = delta_x_t.*( sqrt(R0^2-H^2)/R0 ) - delta_z_t.*(H/R0); % 沿LOS方向，场景中心处的运动误差（空不变误差）
delta_r_R0 = -delta_r_R0*cos(sita_r_c);

%%
% --------------------------------------------------------------------
% 变换到距离频域-方位时域，进行“包络校正，一次运动补偿”
% --------------------------------------------------------------------
S_range = fft(s_echo,NFFT_r,2);     % 进行距离向傅里叶变换，零频在两端。

%{
% 作图
% 图2——距离频域，方位时域，频谱（“未进行包络校正和一次运动补偿”）
figure;
subplot(1,2,1);
imagesc(real(S_range));
title('（a）实部');
xlabel('距离频域（采样点）');
ylabel('方位时域（采样点）');
text(280,-60,'图2，距离频域');       % 给图2进行文字说明
text(250,-15,'未进行包络校正和一次运动补偿');       

subplot(1,2,2);
imagesc(abs(S_range));
title('（b）幅度');
xlabel('距离频域（采样点）');
ylabel('方位时域（采样点）');
%}

% ====================================================
% 计算用于 包络校正 和 一次运动误差 的参考函数
% 包络校正，参考函数
He_fr = exp(1j*4*pi/c.*delta_r_R0.*fr_mtx); % 距离零频在中心
He_fr = fftshift(He_fr,2);                  % 距离零频在两端

% 一次运动补偿，参考函数
Hc1 = exp(1j*4*pi/lamda.*delta_r_R0);   % 距离零频在中心
Hc1 = fftshift(Hc1,2);                  % 距离零频在两端
% ====================================================
% 对距离频谱进行：包络校正 + 一次运动补偿
S_range_f = S_range.*He_fr.*Hc1;        % 零频在两端。    
s_rf = ifft(S_range_f,[],2);            % 完成包络校正和一次运补，回到二维时域。
% s_rf的长度为：Naz*Nrg

%{
% 作图
% 图3——距离频域，方位时域，频谱（已“包络校正 + 一次运动补偿”）
figure;
subplot(1,2,1);
imagesc(real(S_range_f));
title('（a）实部');
xlabel('距离频域（采样点）');
ylabel('方位时域（采样点）');
text(280,-60,'图3，距离频域');       % 给图3进行文字说明
text(250,-15,'已完成：包络校正 + 一次运动补偿');       

subplot(1,2,2);
imagesc(abs(S_range_f));
title('（b）幅度');
xlabel('距离频域（采样点）');
ylabel('方位时域（采样点）');
%}
%{
% 作图
% 图4——二维时域（已“包络校正 + 一次运动补偿”）
figure;
subplot(1,2,1);
imagesc(real(s_rf));
title('（a）实部');
xlabel('距离时域（采样点）');
ylabel('方位时域（采样点）');   
text(280,-60,'图4，二维时域');       % 给图4进行文字说明
text(250,-15,'已完成：包络校正 + 一次运动补偿');       

subplot(1,2,2);
imagesc(abs(s_rf));
title('（b）幅度');
xlabel('距离时域（采样点）');
ylabel('方位时域（采样点）');
%}
%}

%%
% --------------------------------------------------------------------
% 变换到距离多普勒域，进行“补余RCMC”
% --------------------------------------------------------------------
s_rf = s_rf.*exp(-1j*2*pi*fnc.*(ta.'*ones(1,Nrg))); 	% 数据搬移
S_RD = fft(s_rf,NFFT_a,1);  % 进行方位向傅里叶变换，得到距离多普勒域频谱

D_fn_Vr = sqrt(1-lamda^2.*(fa.').^2./(4*Vr^2));     % 徙动因子，列向量
D_fn_Vr_mtx = D_fn_Vr*ones(1,Nrg);  % 形成矩阵，大小：Nrg*Naz

D_fn_ref_Vr = sqrt(1-lamda^2*fn_ref^2/(4*Vr^2));    % 参考频率fn_ref处的徙动因子，是常数。

K_src = 2*Vr^2*f0^3.*D_fn_Vr.^3./(c*R_ref*(fa.').^2);   % 列向量，使用R_ref处的值 
K_src_mtx = K_src*ones(1,Nrg);  % 形成矩阵
Km = Kr./(1-Kr./K_src_mtx);     % 矩阵，这是变换到距离多普勒域的距离调频率。
                                % 使用 R_ref 处的值

% 下面生成 变标方程 s_sc
s_sc = exp(1j*pi.*Km.*(D_fn_ref_Vr./D_fn_Vr_mtx-1).*(tr_mtx-2*R_ref./(c.*D_fn_Vr_mtx)).^2);

% 下面将距离多普勒域的信号与变标方程相乘，实现“补余RCMC”
S_RD_1 = S_RD.*s_sc;        % 相位相乘，实现“补余RCMC”

% 作图
figure;
imagesc(abs(S_RD));
title('包络校正和一次运补后，变换到距离多普勒域，幅度');
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
%
% --------------------------------------------------------------------
% 距离多普勒域，完成“附加相位校正”
% --------------------------------------------------------------------
R0_RCMC = (c/2).*tr;   % 随距离线变化的R0，记为R0_RCMC。
% 附加相位校正项
H2 = exp(-1j*4*pi.*Km./(c^2).*(1-D_fn_Vr_mtx./D_fn_ref_Vr)...
    .*((1./D_fn_Vr)*R0_RCMC-R_ref./D_fn_Vr_mtx).^2); 	% 附加相位校正项
S_RD_3 = S_RD_2.*H2;           % 距离多普勒域，相位相乘，完成“附加相位校正”

%%
%
% --------------------------------------------------------------------
% 回到二维时域，进行二次运动补偿
% --------------------------------------------------------------------
% 下面在二维时域进行二次运动误差补偿
s_rd_rcmc = ifft(S_RD_3,[],1);      % 将距离多普勒域的结果变换到时域
Hc2 = exp(1j*4*pi/lamda.*(delta_r-delta_r_R0)); % 二次运动补偿，参考函数
s_rd_rcmc_MoCo = s_rd_rcmc.*Hc2;    % 进行二次运动补偿

S_rd_rcmc_MoCo = fft(s_rd_rcmc_MoCo,[],1);  % 二次运动补偿后，变换到距离多普勒域

% 作图
figure;
subplot(1,2,1);
imagesc(real(s_rd_rcmc));
title('（a）实部');
xlabel('距离时域（采样点）');
ylabel('方位频域（采样点）');
text(280,-60,'时域');      
text(250,-15,'未进行二次运动补偿前');

subplot(1,2,2);
imagesc(abs(s_rd_rcmc));
title('（b）幅度');
xlabel('距离时域（采样点）');
ylabel('方位频域（采样点）');

figure;
subplot(1,2,1);
imagesc(real(s_rd_rcmc_MoCo));
title('（a）实部');
xlabel('距离时域（采样点）');
ylabel('方位频域（采样点）');
text(280,-60,'时域');      
text(250,-15,'完成了二次运动补偿');

subplot(1,2,2);
imagesc(abs(s_rd_rcmc_MoCo));
title('（b）幅度');
xlabel('距离时域（采样点）');
ylabel('方位频域（采样点）');

figure;
subplot(1,2,1);
imagesc(real(S_rd_rcmc_MoCo));
title('（a）实部');
xlabel('距离时域（采样点）');
ylabel('方位频域（采样点）');
text(280,-60,'距离多普勒');
text(200,-15,'完成了二次运动补偿，还未进行方位MF');       

subplot(1,2,2);
imagesc(abs(S_rd_rcmc_MoCo));
title('（b）幅度');
xlabel('距离时域（采样点）');
ylabel('方位频域（采样点）');
%}

%%
% --------------------------------------------------------------------
% 距离多普勒域，完成“方位压缩”
% --------------------------------------------------------------------
% 生成方位向匹配滤波器
Haz = exp(1j*4*pi.*(D_fn_Vr*R0_RCMC).*f0./c);       % 方位MF

% 下面进行相位相乘，在距离多普勒域，完成方位MF
S_RD_4 = S_rd_rcmc_MoCo.*Haz;      	% 距离多普勒域，相位相乘

% 最后通过IFFT回到图像域，完成方位处理
s_image = ifft(S_RD_4,NFFT_a,1); 	% 完成成像过程，得到成像结果为：s_image

% 作图
figure;
imagesc(abs(S_RD_4));
title('距离多普勒域，进行了相位相乘后（方位MF）');

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
% 由于这里是正侧视，因此有两种方法来计算方位向的位置：
% =========================================================================
% 方式1
% tg_1_x = round(Naz/2);          % 正侧视，目标A的方位向位置就在图像中心。

% 方式2
% 在有斜视角时，方位向位置要通过以下的方式计算。以下程序对于正侧视也兼容，如下：
tg_1_x = rem( R0*tan(sita_r_c)/Vr*Fa , Naz );
if tg_1_x < Naz/2
    tg_1_x = tg_1_x + (Naz/2+1);
else
    tg_1_x = tg_1_x - (Naz/2+1);
end
tg_1_x = round(tg_1_x);    	% 四舍五入，得到整数值，作为点目标的方位中心坐标。
% 这里得到的 tg_1_x 即是点目标中心方位向的位置（坐标）。
% =========================================================================
tg_1_y = round(Nrg/2);
% target_1 = target_analysis( s_image(tg_1_x-NN:tg_1_x+NN,tg_1_y-NN:tg_1_y+NN),Fr,Fa,Vr);


% 目标2，点目标中心在 （tg_2_x，target_2_y）
tg_2_x = tg_1_x + delta_R1/Vr*Fa;
tg_2_y = tg_1_y;
% target_2 = target_analysis( s_image(tg_2_x-NN:tg_2_x+NN,tg_2_y-NN:tg_2_y+NN),Fr,Fa,Vr);


% 目标3，点目标中心在（tg_3_x，tg_3_y）
tg_3_x = tg_2_x + delta_R2*tan(sita_r_c)/Vr*Fa;
tg_3_x = fix(tg_3_x);
tg_3_y = tg_2_y + 2*delta_R2/c*Fr;
% target_3 = target_analysis( s_image(tg_3_x-NN:tg_3_x+NN,tg_3_y-NN:tg_3_y+NN),Fr,Fa,Vr);


% 目标4，点目标中心在（tg_4_x，tg_4_y）
tg_4_x = tg_3_x;
tg_4_x = fix(tg_4_x);
tg_4_y = tg_2_y + 2*(r4 - r2)/c*Fr;
% target_4 = target_analysis( s_image(tg_4_x-NN:tg_4_x+NN,tg_4_y-NN:tg_4_y+NN),Fr,Fa,Vr);






