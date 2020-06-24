function [s_ac,R0_RCMC,Parameter] = RDA_imaging2_v3(raw_data_type)
% 对函数   RDA_imaging_v2（）    进行了以下修改：
% 1）将距离脉压方式改为方式2，并同时设置复制脉冲的参数，使得峰值点被压至脉冲中心；
% 2）该版程序不考虑去除弃置区问题；
% 3）修改方位向匹配滤波器，使得成像结果能够有更好的保相性；
% 将以上修改后的该版本记为：RDA_imaging2_v3（）
%
% 本函数用来对生成的原始数据成像
% 利用：
% 1）带有SRC（利用二维频域相位相乘的方法）的 RD 算法
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
% 2015.01.15. 17:58 p.m.

%%
% 载入原始数据
load test_raw_data_one_point;
% load test_raw_data_one_point2;
% load test_raw_data_100_100;
% load test_raw_data_100_100_2;
% load test_raw_data_3_points;
% load test_raw_data_5_points;
% load test_raw_data_3_5_points;

if raw_data_type == 1
    s_echo = s_echoA;           % 载入天线A对应的原始数据
%     load raw_data2_A;
end
if raw_data_type == 2       
    s_echo = s_echoB;           % 载入天线B对应的原始数据
%     load raw_data2_B;
end

%%
% --------------------------------------------------------------------
% 计算在LOS方向的运动误差，这是进行“运动补偿”的依据
% --------------------------------------------------------------------
r = tr_mtx.*c/2;        % 用于运动误差计算的斜距 r；
% 沿LOS方向的运动误差
delta_r = delta_x_t.*( sqrt(r.^2-H^2)./r ) - delta_z_t.*(H./r); % 沿LOS方向，总的运动误差
delta_r = -delta_r*cos(sita_r_c);
delta_r_R0 = delta_x_t.*( sqrt(R0^2-H^2)/R0 ) - delta_z_t.*(H/R0); % 沿LOS方向，场景中心处的运动误差（空不变误差）
delta_r_R0 = -delta_r_R0*cos(sita_r_c);

%%
% --------------------------------------------------------------------
% 距离压缩，包络校正，一次运动补偿
% --------------------------------------------------------------------
S_range = fft(s_echo,NFFT_r,2);     % 进行距离向傅里叶变换，零频在两端。

%{
% 作图
% 图2——距离频域，方位时域，频谱（未距离压缩）
figure;
subplot(1,2,1);
imagesc(real(S_range));
title('（a）实部');
xlabel('距离频域（采样点）');
ylabel('方位时域（采样点）');
% text(1500,-60,'图2，距离频域');       % 给图2进行文字说明
% text(1700,-20,'未压缩');       
text(280,-60,'图2，距离频域');       % 给图2进行文字说明
text(340,-10,'未压缩');   

subplot(1,2,2);
imagesc(abs(S_range));
title('（b）幅度');
xlabel('距离频域（采样点）');
ylabel('方位时域（采样点）');
%}

% 生成距离向匹配滤波器
% ====================================================
% 采用方式2
% 时域复制脉冲，补零fft，再取复共轭。
% 同时，在设计复制脉冲时，设计一个适当的时延，使得最终的脉压峰值位于LFM脉冲中心；
t_ref = ( -Nr/2 : (Nr/2-1) )/Fr;    % 用来生成距离MF的距离时间轴，长度为Nr（Nr<Nrg）；
t_ref_mtx = ones(Naz,1)*t_ref;      % 矩阵形式
w_range_ref = (abs(t_ref_mtx)) <= ((Tr/2).*ones(Naz,Nr));
% 距离向包络，即距离窗，限制脉冲长度
w_ref = kaiser(Nr,2.5);             % 距离向，构建Kaiser窗，此为列向量。
w_ref = ones(Naz,1)*(w_ref.');      % 构成矩阵形式，每一行都相同的加窗。

s_ref = w_range_ref.*exp((1j*pi*Kr).*((t_ref_mtx).^2)); % 复制（发射）脉冲，未加窗。
% s_ref = w_ref.*w_range_ref.*exp((1j*pi*Kr).*((t_ref_mtx).^2)); % 复制（发射）脉冲，加了窗。

% ************************************************************************
% 对s_ref设计一个时间延迟，使得峰值点被压至 LFM 脉冲中心；
% 同时，这部分也完成对s_ref的补零；
s_ref = [s_ref,zeros(Naz,Nrg-Nr)];      % 对复制脉冲，末端补零至长度 Nrg；
s_ref = circshift(s_ref,[0 -Nr/2]);     % 向左循环移位（Nr/2）；
% 要设计将峰值点压至LFM脉冲中心，既可以用上面的两行代码；
% 也可以用下面的这两行；
% s_ref = fftshift(s_ref,2);
% s_ref = [s_ref(:,1:Nr/2),zeros(Naz,Nrg-Nr),s_ref(:,Nr/2+1:end)];
% ************************************************************************

S_ref = fft(s_ref,NFFT_r,2);            % 复制脉冲的距离傅里叶变换，零频在两端。
H_range = conj(S_ref);                  % 距离向匹配滤波器，零频在两端。
% ====================================================
% 计算用于 包络校正 和 一次运动误差 的参考函数

% 包络校正，参考函数
He_fr = exp(1j*4*pi/c.*delta_r_R0.*fr_mtx); % 距离零频在中心
He_fr = fftshift(He_fr,2);  % 距离零频在两端

% 一次运动补偿，参考函数
Hc1 = exp(1j*4*pi/lamda.*delta_r_R0);   % 距离零频在中心
Hc1 = fftshift(Hc1,2);      % 距离零频在两端
% ====================================================
% 对距离频谱进行：距离压缩 + 包络校正 + 一次运动补偿
S_range_c = S_range.*H_range;
% S_range_c = S_range.*H_range.*He_fr.*Hc1;	% 零频在两端。      
s_rc = ifft(S_range_c,[],2);            % 完成距离压缩，包络校正和一次运补，回到二维时域。
% s_rc的长度为：Naz*Nrg。未去除弃置区。
% 这里都统一不考虑去除弃置区问题
% ====================================================

%{
% 作图
% 图3——距离频域，方位时域，频谱（已“距离压缩 + 包络校正 + 一次运动补偿”）
figure;
subplot(1,2,1);
imagesc(real(S_range_c));
title('（a）实部');
xlabel('距离频域（采样点）');
ylabel('方位时域（采样点）');
% text(1500,-60,'图3，距离频域');       % 给图3进行文字说明
% text(1300,-20,'已完成：距离压缩 + 包络校正 + 一次运动补偿');       
text(280,-60,'图3，距离频域');       % 给图3进行文字说明
text(230,-15,'已完成：距离压缩 + 包络校正 + 一次运动补偿');       

subplot(1,2,2);
imagesc(abs(S_range_c));
title('（b）幅度');
xlabel('距离频域（采样点）');
ylabel('方位时域（采样点）');
%}
%{
% 作图
% 图4——二维时域（已“距离压缩 + 包络校正 + 一次运动补偿”）
figure;
subplot(1,2,1);
imagesc(real(s_rc));  %　这及其以下，都直接使用去除弃置区后的结果
title('（a）实部');
xlabel('距离时域（采样点）');
ylabel('方位时域（采样点）');
% text(1350,-60,'图4，二维时域');       % 给图4进行文字说明
% text(1250,-20,'已完成：距离压缩 + 包络校正 + 一次运动补偿');       
text(150,-60,'图4，二维时域');       % 给图4进行文字说明
text(140,-15,'已完成：距离压缩 + 包络校正 + 一次运动补偿');       

subplot(1,2,2);
imagesc(abs(s_rc));
title('（b）幅度');
xlabel('距离时域（采样点）');
ylabel('方位时域（采样点）');
%}

%%
% --------------------------------------------------------------------
% 变换到二维频域，进行SRC
% --------------------------------------------------------------------
s_rc = s_rc.*exp(-1j*2*pi*fnc.*(ta.'*ones(1,Nrg)));    % 数据搬移
S_2df = fft(s_rc,NFFT_a,1);    	% 方位向傅里叶变换，到距离多普勒域。
S_2df = fft(S_2df,NFFT_r,2);   	% 距离向傅里叶变换，到二维频域
% ！！！注意：距离向零频在两端。！！！
% ====================================================================
% 设置方位频率轴——这是关键点
fa = fnc + fftshift(-NFFT_a/2:NFFT_a/2-1)/NFFT_a*Fa; 	% 方位频率轴如此设置。
% =====================================================================
D_fn_Vr = sqrt(1-lamda^2.*(fa.').^2./(4*Vr^2));         % 大斜视角下的徙动因子
K_src = 2*Vr^2*f0^3.*D_fn_Vr.^3./(c*R0*(fa.').^2);      % 列向量
K_src_1 = 1./K_src;     	% 列向量。为了后面能使用矩阵乘法，这里先求倒数
H_src = exp(-1j*pi.*K_src_1*(fr.^2)); 	% 二次距离压缩滤波器。距离向，零频在中间。
% 这是矩阵，大小Naz*Nrg
H_src = fftshift(H_src,2); 	% （左右半边互换）距离向，零频在两端。 ！！！这很关键！！！
S_2df_src = S_2df.*H_src;	% 这一步点乘时，要注意两者的距离向频率轴应该对应上，不然会出错！！
% 这就是为什么上面的 H_src 要 fftshift 的原因！！

S_rd = ifft(S_2df_src,[],2);    	% 完成二次距离压缩（SRC），回到距离多普勒域。

%{
% 作图
figure;
imagesc(abs(S_rd));
title('SRC后，距离多普勒域');
%}

%%
% --------------------------------------------------------------------
% 距离多普勒域，进行距离徙动校正
% --------------------------------------------------------------------
% 每一个最近斜距（R0）都随着距离门的不同而改变。
tr_RCMC = 2*R0/c + ( -Nrg/2 : (Nrg/2-1) )/Fr;   % 用于RCMC的时间轴（其实和tr相同）；
R0_RCMC = (c/2).*tr_RCMC;   % 随距离线变化的R0，记为R0_RCMC，用来计算RCM和Ka。
delta_Rrd_fn = ((1-D_fn_Vr)./D_fn_Vr)*R0_RCMC; 	% 更精确的 RCM 值。

num_range = c/(2*Fr);   % 一个距离采样单元，对应的长度
delta_Rrd_fn_num = delta_Rrd_fn./num_range; % 每一个方位向频率，其RCM对应的距离采样单元数

R = 8;  % sinc插值核长度
S_rd_rcmc = zeros(NFFT_a,Nrg); 	% 用来存放RCMC后的值

h = waitbar(0,'RCMC, please wait');
for p = 1 : NFFT_a
    for q = 1 : Nrg   % 此时距离向的长度是 Nrg。    
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
    waitbar(p/NFFT_a);
end
close(h);
% S_rd_rcmc 就是RCMC后的距离多普勒域频谱。

%{
% 作图
% 图5——距离多普勒域（未RCMC）
figure;
subplot(1,2,1);
imagesc(real(S_rd));
title('（a）实部');
xlabel('距离时域（采样点）');
ylabel('方位频域（采样点）');
% text(1350,-60,'图5，距离多普勒域');       % 给图5进行文字说明
% text(1550,-20,'未RCMC');     
text(150,-60,'图5，距离多普勒域');       % 给图5进行文字说明
text(172,-10,'未RCMC'); 

subplot(1,2,2);
imagesc(abs(S_rd));
title('（b）幅度');
xlabel('距离时域（采样点）');
ylabel('方位频域（采样点）');
%}
%{
% 作图
% 图6——距离多普勒域，RCMC后的结果
figure;
subplot(1,2,1);
imagesc(real(S_rd_rcmc));
title('（a）实部');
xlabel('距离时域（采样点）');
ylabel('方位频域（采样点）');
% text(1350,-60,'图6，距离多普勒域');       % 给图6进行文字说明
% text(1550,-20,'已RCMC');    
text(150,-60,'图6，距离多普勒域');       % 给图6进行文字说明
text(172,-10,'已RCMC');

subplot(1,2,2);
imagesc(abs(S_rd_rcmc));
title('（b）幅度');
xlabel('距离时域（采样点）');
ylabel('方位频域（采样点）');
%}

%%
%{
% --------------------------------------------------------------------
% 回到时域，进行二次运动补偿
% --------------------------------------------------------------------
% 由于去除了弃置区，因此这里用于二次运动补偿的斜距 r 要重新计算
r = ones(Naz,1)*R0_RCMC;
% 下面计算对应去除弃置区后，引入的运动误差
delta_x_t = a*sin(2*pi*w/La*Vr.*(ta.'*ones(1,Nrg)));   % 这是沿地距 x 轴的运动误差
% delta_x_t = 0;
% delta_z_t = a*sin(2*pi*w/La*Vr.*(ta.'*ones(1,Nrg)));   % 这是沿 z 轴的运动误差
delta_z_t = 0;

% 计算，沿LOS方向的运动误差
delta_r = delta_x_t.*( sqrt(r.^2-H^2)./r ) - delta_z_t.*(H./r);% 沿LOS方向，总的运动误差
delta_r = -delta_r*cos(sita_r_c);
delta_r_R0 = delta_x_t.*( sqrt(R0^2-H^2)/R0 ) - delta_z_t.*(H/R0); % 沿LOS方向，场景中心处的运动误差（空不变误差）
delta_r_R0 = -delta_r_R0*cos(sita_r_c);

% 下面在二维时域进行二次运动误差补偿
Hc2 = exp(1j*4*pi/lamda.*(delta_r-delta_r_R0)); % 二次运动补偿，参考函数
s_rd_rcmc = ifft(S_rd_rcmc,[],1);   % 将RCMC后的结果变换到时域
s_rd_rcmc_MoCo = s_rd_rcmc.*Hc2;    % 进行二次运动补偿

S_rd_rcmc_MoCo = fft(s_rd_rcmc_MoCo,[],1);  % 二次运动补偿后，变换到距离多普勒域

figure;
subplot(1,2,1);
imagesc(real(S_rd_rcmc_MoCo));
title('（a）实部');
xlabel('距离时域（采样点）');
ylabel('方位频域（采样点）');
% text(1350,-60,'距离多普勒');
% text(1250,-20,'完成了二次运动补偿');       
text(150,-60,'距离多普勒');
text(172,-15,'完成了二次运动补偿');      

subplot(1,2,2);
imagesc(abs(S_rd_rcmc_MoCo));
title('（b）幅度');
xlabel('距离时域（采样点）');
ylabel('方位频域（采样点）');
%}

%%
% --------------------------------------------------------------------
% 方位压缩
% --------------------------------------------------------------------
fa_azimuth_MF = fa;         % 方位频率轴，采用和RCMC中所用的频率轴相同。
% Haz = exp(1j*4*pi.*(D_fn_Vr*R0_RCMC).*f0./c);       % 更精确的的方位向脉压滤波器
Haz = exp(1j*4*pi.*((D_fn_Vr-1)*R0_RCMC).*f0./c);   % 将 D 换成（D-1）；
% 这里的频率轴和距离多普勒域的方位频谱是对应的。

S_rd_c = S_rd_rcmc.*Haz;    	% 乘以匹配滤波器（这是没有进行二次运动补偿时采用的）
% S_rd_c = S_rd_rcmc_MoCo.*Haz; 	% 乘以匹配滤波器
s_ac = ifft(S_rd_c,[],1);    	% 完成方位压缩，变到图像域。结束。

% 作图
% 图7——成像结果
figure;
imagesc(abs(s_ac));
title('点目标成像');
xlabel('距离时域（采样点）');
ylabel('方位时域（采样点）'); 

%%
% 需要返回的一些参数
Parameter = [ Naz;          % Parameter 的第一行代表 Naz
            H;              % Parameter 的第二行代表 H
            lamda ];        % Parameter 的第三行代表 lamda　


end



