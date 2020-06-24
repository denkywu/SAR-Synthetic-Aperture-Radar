function PHY_s_after_filtering = Filtering_local_frequency_estimate(s_after_flat_earth,fnc,NN,N_UP,window_M,window_N)
% 基于“局部频率估计”的自适应干涉图滤波——版本 1
% 该版本选取各自不重叠的分析窗；
% 
% 输入：
% 1）s_after_flat_earth  是去平地效应后的复数据，包括幅度和相位；
% 2) fnc         是多普勒中心频率；
% 3）NN          是对数据进行分块的大小，一般选择为16；
% 4）N_UP        是对分块数据进行局部频率估计时，二维补零FFT的大小；
%              	 当NN选择为16时，N_UP应该选择为256；
% 5） window_M 和 window_N    表示所选定的窗口大小；
%                注意：（2*window_M+1） 的大小要小于分析窗的大小（NN）；
%                注意：真正的滤波窗口大小是（2*window_M+1）；
%                      故如果 window_M = 5，则滤波窗口是 11；
%
% 输出：
% 1）PHY_s_after_avg_filtering 是滤波后的相位图
%
% 截止到：2015.02.04. 11:20 a.m.

%%
disp('正在进行基于“局部频率估计”的自适应干涉图滤波，请等待');
% -----------------------------------------------------------------------
%                  基于“局部频率估计”的自适应干涉图滤波                   
% -----------------------------------------------------------------------
[Naz,Nrg] = size(s_after_flat_earth);
% Naz 和 Nrg 必须是 NN 的整数倍；
s_after_filtering = zeros(Naz,Nrg);

for p = 1 : NN : Naz
    for q = 1 : NN : Nrg
        % 截取一段（从去平地相位后的干涉复数据中截取）
        tmp = s_after_flat_earth(p:p+NN-1,q:q+NN-1); % 大小 NN×NN；
        % tmp 称为一个分析窗中的复数据，是我们的处理对象
        
        % -----------------------------------------------------------------
        %             下面进行“局部频率估计”得到近似的线性相位条纹          
        % -----------------------------------------------------------------
        % 下面首先进行“局部频率估计”得到近似的线性相位条纹，为后续自适应滤波作准备。
        % 补零，准备进行二维FFT
        tmp1 = zeros(N_UP,N_UP);
        tmp1(1:NN,1:NN) = tmp;      % 补零至 N_UP×N_UP（256×256）；
        TMP1 = fft2(tmp1);          % 进行二维 FFT
        TMP1 = fftshift(TMP1);      % 进行二维 fftshift，使得距离向和方位向的零频都位于中心；
        % 求幅度谱最大值位置
        [p1,q1] = max(abs(TMP1));
        [p2,q2] = max(max(abs(TMP1)));
        % 峰值位于：（ q1(q2)，(q2) ）；  
        clear tmp1;clear TMP1;
        % 下面根据峰值位置，得到所估计的方位向和斜距向频率；
        % 斜距向频率轴
        fr_tmp = (-N_UP/2:N_UP/2-1)*(1/N_UP);       % 采用归一化频率轴
        fr_max = fr_tmp(q2);            % 这是估计得到的距离向频率
        % 方位向频率轴
        fa_tmp = fnc + (-N_UP/2:N_UP/2-1)*(1/N_UP); % 采用归一化频率轴
        fa_max = fa_tmp(q1(q2));        % 这是估计得到的方位向频率
        clear fr_tmp;clear fa_tmp;
        % 以上面得到的距离向、方位向频率为基础；
        % 认为这一段小数据的相位近似是线性相位，其频率是同一点频，即上面估计所得；
        % 然后下面构建线性相位条纹：
        tmp2 = zeros(NN,NN);   % 用来存放该小块数据对应的线性相位条纹数据（复数据）
        jj = [1:1+NN-1].'*ones(1,NN);
        kk = ones(NN,1)*[1:1+NN-1];
        tmp2 = exp(1j*2*pi.*(fa_max*jj + fr_max*kk));
        % 计算初始相位
        phy0 = angle(sum(sum(tmp.*exp(-1j*2*pi.*(fa_max*jj + fr_max*kk)))));
        tmp2 = tmp2.*exp(1j*phy0);   
        % 这就是我们需要的线性相位条纹；
        % tmp2 是复数据，我们实际需要的是其相位 angle(tmp2)；
        % 但用复数据便于后面的计算；
        
        % -----------------------------------------------------------------
        %                       下面进行自适应滤波          
        % -----------------------------------------------------------------        
        % 一个分析窗中的复数据，我们认为由三部分组成：
        % 1）是反映条纹的大致架构的线性相位，由上面的局部频率估计并重构得到；
        % 2）反映地形细节的相位信息；
        % 3）相位噪声；
        % 我们在得到了1）后，从tmp中去除线性相位并进行滤波，即认为得到了2）；
        M = window_M;              	% 距离向滤波窗口大小；
        N = window_N;              	% 方位向滤波窗口大小；
        tmp3 = tmp.*conj(tmp2);     % 这是去除了线性相位后待滤波的数据；
        % 下面进行滤波
        for pp = 1:NN
            for qq = 1:NN
                % 首先进行条件判断，看窗口window是否超过了矩阵的边界：
                if pp<(N+1) || pp>(NN-N) || qq<(M+1) || qq>(NN-M)
                    % 若满足上述条件中的任何一个，说明窗口位于矩阵边界，进行以下进一步判断
                    if (pp-N) < 1
                        x_min = 1;
                    else
                        x_min = pp - N;
                    end
                    if (pp+N) > NN
                        x_max = NN;
                    else
                        x_max = pp + N;
                    end
                    if (qq-M) < 1
                        y_min = 1;
                    else
                        y_min = qq - M;
                    end
                    if (qq+M) > NN
                        y_max = NN;
                    else
                        y_max = qq + M;
                    end
                    s_window = tmp3(x_min:x_max,y_min:y_max);
                else
                    % 若上述四个条件都不满足，说明窗口不位于矩阵边界，则可以取到全部
                    % （2N+1）*（2M+1）个点，因此直接用以下命令即可
                    s_window = tmp3(pp-N:pp+N,qq-M:qq+M);
                end
                tmp3_after_filtering(pp,qq) = exp(1j*angle(sum(sum(s_window))));
            end
        end
        % 对 tmp3 的滤波完成，得到滤波后的结果：tmp3_after_filtering
        % 下面将 tmp2 和 tmp3_after_filtering 相乘，就可以得到最终对tmp的滤波结果
        s_after_filtering(p:p+NN-1,q:q+NN-1) = tmp2.*tmp3_after_filtering; 
    end
    p
end

PHY_s_after_filtering = angle(s_after_filtering);
% 这就是滤波后的相位图，也是该函数的返回值； 
disp('相位滤波已完成');

end







