function [s_image,R0_RCMC,Parameter] = CSA_imaging(raw_data_type)
% 本函数用来对生成的原始数据成像
% 利用：
% 1）CSA 算法
% 2）可选择是否进行两步式运动误差补偿
%
% 载入原始数据时：
% 1）raw_data_type == 1，代表天线 A 对应的原始数据；
% 2）raw_data_type == 2，代表天线 B 对应的原始数据；
%
% 返回值：
% 1）s_ac        是成像结果；
% 2）R0_RCMC     是随距离线变化的斜距，用于后面计算平地相位；
% 2）Parameter   代表成像结果中的一些参数，用于后面进行平地相位计算；
%
% 程序版本截止到：
% 2014.12.23. 9:45 a.m.

%%
% 载入原始数据
if raw_data_type == 1
%     load raw_data1_A;           % 载入天线A对应的原始数据
    load raw_data2_A;
end
if raw_data_type == 2       
%     load raw_data1_B;           % 载入天线B对应的原始数据
    load raw_data2_B;
end

R_ref = R0;             % 参考目标选在场景中心，其最近斜距为 R_ref
fn_ref = fnc;        	% 参考目标的多普勒中心频率

%%
%{
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
s_rf = s_echo;
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

%{
% 作图
figure;
imagesc(abs(S_RD));
title('包络校正和一次运补后，变换到距离多普勒域，幅度');
figure;
imagesc(abs(S_RD_1));
title('距离多普勒域，补余RCMC后，幅度');
%}

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

%{
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
%}

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
%{
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
S_rd_rcmc_MoCo = S_RD_3;
% --------------------------------------------------------------------
% 距离多普勒域，完成“方位压缩”
% --------------------------------------------------------------------
% 生成方位向匹配滤波器
Haz = exp(1j*4*pi.*(D_fn_Vr*R0_RCMC).*f0./c);       % 方位MF

% 下面进行相位相乘，在距离多普勒域，完成方位MF
S_RD_4 = S_rd_rcmc_MoCo.*Haz;      	% 距离多普勒域，相位相乘

% 最后通过IFFT回到图像域，完成方位处理
s_image = ifft(S_RD_4,NFFT_a,1); 	% 完成成像过程，得到成像结果为：s_image

s_image = fftshift(s_image,1);

% 作图
%{
figure;
imagesc(abs(S_RD_4));
title('距离多普勒域，进行了相位相乘后（方位MF）');
%}

figure;
imagesc(abs(s_image));
title('成像结果');

%%
% 需要返回的一些参数
Parameter = [ Naz;          % Parameter 的第一行代表 Naz
            H;              % Parameter 的第二行代表 H
            lamda ];        % Parameter 的第三行代表 lamda　



end
