%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                       Radarsat_1 光盘中数据
%                             CSA 成像
%
%
%                               WD
%                       2014.10.31. 23:40 p.m.
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
% 注：
%   这里还附加了一步：对原始数据进行补零，再进行后续处理。（补零是非常必要的！！）
%
% （2） 本文件夹中还有一个文件：CD_run_params
%           ——这里面是仿真中需要用的许多参数，直接载入即可。
%
% （3）成像程序说明：
%       由CSA的点目标程序修改而来；
%
% （4）成像流程：
%   ——原始数据
%   ——经过方位向FFT，变换到距离多普勒域，进行“补余RCMC”
%   ——经过距离向FFT，变换到二维频域，进行“距离压缩”、“SRC”、“一致RCMC”
%   ——经过距离向IFFT，变换到距离多普勒域，进行“方位压缩”和“附加相位校正”
%   ——经过方位向IFFT，回到图像域，成像结束。
%
% （5）分别采用了两种方式抑制相干斑：
%       a）“利用3*3的窗口，进行滑动平均，以抑制相干斑”；
%       b）进行“4视叠加”，来抑制相干斑；
%
% 本程序修改截止到： 2014.10.31. 23:40 p.m.

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
% clear S_RD_3;     % 后面进行多视叠加需要距离多普勒域频谱，因此这里不要清除。

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

if b == 1 || b == 2 
    % 对单一分区（比如，分区1或者分区2），使用这部分程序来显示图像
    figure;
    imagesc(((0:Nrg-1)+first_rg_cell)/Fr*c/2+R0,((0:Naz-1)+first_rg_line)/Fa*Vr,fftshift(s_tmp,1),clim);
    axis xy;
    title('RADARSAT-1数据，使用CS算法，成像结果')
    xlabel('Range(m)')
    ylabel('Azimuth(m)')
    colormap(gray)
end

if b == 3
    % 对两个分区一起成像时，使用这部分来显示图像。
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
    colormap(gray)
end


%%
% --------------------------------------------------------------------
% 利用 3*3 的窗口，进行滑动平均
% 以此来抑制相干斑噪声
% --------------------------------------------------------------------
s_image_look = zeros(Naz,Nrg);      % 用来存放滑动平均后的结果

h = waitbar(0,'利用3*3的窗口，进行滑动平均，please wait');
for p = 1 : Naz
    for q = 1 : Nrg
        count = 0;
        s_image_look(p,q) = s_image_look(p,q) + abs(s_image(p,q));
        count = count + 1;
        if p>2 && p<(Naz-1) && q>2 && q<(Nrg-1)
            s_image_look(p,q) = s_image_look(p,q)+abs(s_image(p-1,q-1))+...
                abs(s_image(p-1,q))+abs(s_image(p-1,q+1))+abs(s_image(p,q-1))+...
                abs(s_image(p,q+1))+abs(s_image(p+1,q-1))+abs(s_image(p+1,q+1))+...
                abs(s_image(p+1,q+1));
            count = 9;
        else
            if (p-1) > 0
                s_image_look(p,q) = s_image_look(p,q) + abs(s_image(p-1,q));
                count = count + 1;
                if (q-1) > 0
                    s_image_look(p,q) = s_image_look(p,q) + abs(s_image(p-1,q-1));
                    count = count+1;
                end
                if (q+1) <= Nrg
                    s_image_look(p,q) = s_image_look(p,q) + abs(s_image(p-1,q+1));
                    count = count+1;
                end
            end

            if (p+1) <= Naz
                s_image_look(p,q) = s_image_look(p,q) + abs(s_image(p+1,q));
                count = count + 1;
                if (q-1) > 0
                    s_image_look(p,q) = s_image_look(p,q) + abs(s_image(p+1,q-1));
                    count = count+1;
                end
                if (q+1) <= Nrg
                    s_image_look(p,q) = s_image_look(p,q) + abs(s_image(p+1,q+1));
                    count = count+1;
                end
            end
            
            if (q-1) > 0
                s_image_look(p,q) = s_image_look(p,q) + abs(s_image(p,q-1));
                count = count+1;
            end
            
            if (q+1) <= Nrg
                s_image_look(p,q) = s_image_look(p,q) + abs(s_image(p,q+1));
                count = count+1;
            end
        end
        s_image_look(p,q) = s_image_look(p,q)/count;    % 取平均
    end
    waitbar(p/Naz);
end
close(h);
% 到此为止，我们得到了滑动平均后的结果： s_image_look
disp('利用 3*3 的窗口，进行滑动平均，得到抑制相干斑后的结果');

% 下面进行显示
clear sout;clear G;clear clim;
clear tmp;clear s_tmp;clear ss_tmp
% 抑制相干斑后的结果
% 对亮度进行非线性变换，减小对比度，显示图像
sout = abs(s_image_look)/max(max(abs(s_image_look)));
G = 20*log10(sout+eps);             % dB显示
clim = [-55 0];                     % 动态显示范围

tmp = round(2*(R0/D_fn_ref_Vr-R0)/c*Fr);
s_tmp(:,1:Nrg-tmp+1) = G(:,tmp:end);
s_tmp(:,Nrg-tmp+1+1:Nrg) = G(:,1:tmp-1);
if b == 1 || b == 2
    % 对单一分区（比如，分区1或者分区2），使用这部分程序来显示图像
    figure;
    imagesc(((0:Nrg-1)+first_rg_cell)/Fr*c/2+R0,((0:Naz-1)+first_rg_line)/Fa*Vr,fftshift(s_tmp,1),clim);
    axis xy;
    title('滑动平均后的结果');
    xlabel('Range(m)')
    ylabel('Azimuth(m)')
    colormap(gray)
end

if b == 3
    % 对两个分区一起成像时，使用这部分来成像。
    % 作用是：将上下部分进行一定的移位
    %       （ 原来的图像的第2900行到最后一行应该在新图像的最开头 ）
    ss_tmp(1:Naz-2900+1,:) = s_tmp(2900:Naz,:);
    ss_tmp(Naz-2900+1+1:Naz,:) = s_tmp(1:2900-1,:);
    figure;
    imagesc(((0:Nrg-1)+first_rg_cell)/Fr*c/2+R0,((0:Naz-1)+first_rg_line)/Fa*Vr,ss_tmp,clim);
    axis xy;
    title('滑动平均后的结果');
    xlabel('Range(m)')
    ylabel('Azimuth(m)')
    colormap(gray)
end


%% 
% --------------------------------------------------------------------
% 进行“4视叠加”
% 抑制相干斑
% --------------------------------------------------------------------
S_rd_c_fftshift = fftshift(S_RD_3,1);

% 子视1
if b == 1 || b == 2
    S_rd_c_1 = S_rd_c_fftshift(1:569,:);
end
if b == 3
    S_rd_c_1 = S_rd_c_fftshift(1:1135,:);
end
s_ac1 = ifft(S_rd_c_1);
% figure;
% imagesc(abs(s_ac1));
% title('子视1');

% 子视2
if b == 1 || b == 2
    S_rd_c_2 = S_rd_c_fftshift(494:1062,:);
end
if b == 3
    S_rd_c_2 = S_rd_c_fftshift(988:2122,:);
end
s_ac2 = ifft(S_rd_c_2);
% figure;
% imagesc(abs(s_ac2));
% title('子视2');

% 子视3
if b == 1 || b == 2
    S_rd_c_3 = S_rd_c_fftshift(987:1555,:);
end
if b == 3
    S_rd_c_3 = S_rd_c_fftshift(1975:3109,:);
end
s_ac3 = ifft(S_rd_c_3);
% figure;
% imagesc(abs(s_ac3));
% title('子视3');

% 子视4
if b == 1 || b == 2
    S_rd_c_4 = S_rd_c_fftshift(1480:2048,:);
end
if b == 3
    S_rd_c_4 = S_rd_c_fftshift(2962:4096,:);
end
s_ac4 = ifft(S_rd_c_4);
% figure;
% imagesc(abs(s_ac4));
% title('子视4');

% 子视求和
s_ac_look = sqrt(abs(s_ac1).^2 + abs(s_ac2).^2 + abs(s_ac3).^2 + abs(s_ac4).^2);   
% 对子视幅度平方求和再开方
% 到此为止，我们得到了4视叠加后的结果： s_ac_look
disp('进行了4视叠加，得到抑制相干斑后的结果');


% 下面进行显示
clear sout;clear G;clear clim;
clear tmp;clear s_tmp;clear ss_tmp
% 抑制相干斑后的结果
% 对亮度进行非线性变换，减小对比度，显示图像
sout = abs(s_ac_look)/max(max(abs(s_ac_look)));
G = 20*log10(sout+eps);             % dB显示
clim = [-55 0];                     % 动态显示范围

[Naz1,Nrg1] = size(s_ac_look);
tmp = round(2*(R0/D_fn_ref_Vr-R0)/c*Fr);
s_tmp(:,1:Nrg1-tmp+1) = G(:,tmp:end);
s_tmp(:,Nrg1-tmp+1+1:Nrg1) = G(:,1:tmp-1);
if b == 1 || b == 2
    % 对单一分区（比如，分区1或者分区2），使用这部分程序来显示图像
    figure;
    imagesc(((0:Nrg1-1)+first_rg_cell)/Fr*c/2+R0,((0:Naz1-1)+first_rg_line)/Fa*Vr,fftshift(s_tmp,1),clim);
    axis xy;
    title('4视叠加后的结果');
    xlabel('Range(m)')
    ylabel('Azimuth(m)')
    colormap(gray)
end

if b == 3
    % 对两个分区一起成像时，使用这部分来成像。
    % 作用是：将上下部分进行一定的移位
    %       （ 原来的图像的第780行到最后一行应该在新图像的最开头 ）
    ss_tmp(1:Naz1-780+1,:) = s_tmp(780:Naz1,:);
    ss_tmp(Naz1-780+1+1:Naz1,:) = s_tmp(1:780-1,:);
    figure;
    imagesc(((0:Nrg1-1)+first_rg_cell)/Fr*c/2+R0,((0:Naz1-1)+first_rg_line)/Fa*Vr,ss_tmp,clim);
    axis xy;
    title('4视叠加后的结果');
    xlabel('Range(m)')
    ylabel('Azimuth(m)')
    colormap(gray)
end




