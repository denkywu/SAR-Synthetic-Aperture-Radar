function [PSLR_r,ISLR_r,IRW_r, PSLR_a,ISLR_a,IRW_a] = target_analysis_2(s_ac,Fr,Fa,Vr)
    % 截止到 2014.10.10 10:34 a.m.
    % 
    % 2014.10.10. 修改了： 指标计算函数采用更改后的: zhibiao_2( )
    % 2014.10.10. 修改了： 二维补零方式更改为第 2 种方法：先对每列补零，再对每行补零。
    %
    % ===============================================================
    % 该程序考虑了大斜视角下，对点目标的扭曲
    % 分别考虑距离向和方位向
    % 考虑距离向旋转角度，将距离向轴旋转到平行于水平轴，再取出行切片
    % 考虑方位向旋转角度，将方位向轴旋转到平行于垂直轴，再取出列切片
    % ===============================================================
    %
    % 输入变量： s_ac  ——是需要进行指标分析的点目标点k，s_ac是它所在的一块矩阵数据。
    % 输入变量： Fr  —— 距离向采样率
    % 输入变量： Fa  —— 方位向采样率
    % 输入变量： Vr  —— 平台速度
    % 该函数用来求解点目标k的各项指标
    % 点目标中心的行切片（距离向）：峰值旁瓣比(PSLR)，积分旁瓣比(ISLR),距离分辨率（IRW）。
    % 点目标中心的列切片（方位向）：峰值旁瓣比(PSLR)，积分旁瓣比(ISLR),距离分辨率（IRW）。
    % 输出值分别是行切片和列切片的各项指标
    %
    % 补零方式很关键
    % 这里采用的方法是逐行，逐列判断数据最小值位置，然后将在最小值位置处进行补零
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %       取点目标中心 NN*NN 切片，进行升采样
    %   再取点目标中心，距离向切片，方位向切片，进行升采样
    %                   计算各项性能指标
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%
    c = 3e8;                % 光速
    
    % 取点目标的中心切片，NN*NN
    % 进行二维升采样
    NN = 32;        % 切片总长度，NN*NN
    [row,column] = size(s_ac);      % s_ac的矩阵大小
    [aa,p] = max(abs(s_ac));           
    [bb,q] = max(max(abs(s_ac)));
     
    row_max = p(q);    	% 二维矩阵最大值，所在第几行—— row_max。
    column_max = q;     % 二维矩阵最大值，所在第几列—— column_max。
    s_ac_max = bb;      % 矩阵最大值是—— x_max。
    
    s_ac_test = s_ac(row_max-NN/2:row_max+NN/2-1,column_max-NN/2:column_max+NN/2-1);
    % 得到NN*NN的切片
    
    % 下面进行二维升采样
    S_ac_test_1 = fft(s_ac_test,[],1);     % 方位向fft
    S_ac_test_2 = fft(S_ac_test_1,[],2);   % 距离向fft
    
    % 接下来进行二维补零
    % =========================================================================
    % 利用两个for循环，在每行和每列最小值的位置处补上7*NN的零，实现高频补零。
    %------------------------------------------------------------------
    %{
    % 第 1 种方法： 先对每行补零，再对每列补零。
    S_ac_test_buling_1 = zeros(NN,8*NN);  % 中间变量
    S_ac_test_buling = zeros(8*NN,8*NN);
    % 下面先对每行补零，再对每列补零
    for pp = 1:NN           % 在每行的最小值位置补零
        [C,I] = min(S_ac_test_2(pp,:));
        S_ac_test_buling_1(pp,1:I) = S_ac_test_2(pp,1:I);
        S_ac_test_buling_1(pp,8*NN-(NN-I)+1:8*NN) = S_ac_test_2(pp,I+1:NN);
    end
    % for qq = 1:8*NN         % 在每列的最小值位置补零
    %     [C,I] = min(S_ac_test_buling_1(:,qq));
    %     S_ac_test_buling(1:I,qq) = S_ac_test_buling_1(1:I,qq);
    %     S_ac_test_buling(8*NN-(NN-I)+1:8*NN,qq) = S_ac_test_buling_1(I+1:NN,qq);
    % end
    % 由于上面的在每列的最小值位置补零，得到的结果是有问题的，因此下面进行修改
    % 修改方法是：
    %       直接对每一列，在中间补零
    S_ac_test_buling(1:NN/2,:) = S_ac_test_buling_1(1:NN/2,:);
    S_ac_test_buling(8*NN-NN/2+1:8*NN,:) = S_ac_test_buling_1(NN/2+1:NN,:);
    %}
    %------------------------------------------------------------------
    %
    % 第 2 种方法： 先对每列补零，再对每行补零。
    S_ac_test_buling_1 = zeros(8*NN,NN);  % 中间变量
    S_ac_test_buling = zeros(8*NN,8*NN);
    % 下面先对每列补零，再对每行补零——注意：和上面的先行后列相比，这是有区别的。
    for pp = 1:NN           % 在每行的最小值位置补零
        [C,I] = min(S_ac_test_2(:,pp));
        S_ac_test_buling_1(1:I,pp) = S_ac_test_2(1:I,pp);
        S_ac_test_buling_1(8*NN-(NN-I)+1:8*NN,pp) = S_ac_test_2(I+1:NN,pp);
    end
    for qq = 1:8*NN         % 在每列的最小值位置补零
        [C,I] = min(S_ac_test_buling_1(qq,:));
        S_ac_test_buling(qq,1:I) = S_ac_test_buling_1(qq,1:I);
        S_ac_test_buling(qq,8*NN-(NN-I)+1:8*NN) = S_ac_test_buling_1(qq,I+1:NN);
    end   
    %}
    %------------------------------------------------------------------
    % =========================================================================
    S_ac_test_1 = ifft(S_ac_test_buling,[],2);
    s_ac_test = ifft(S_ac_test_1,[],1);         % 完成二维升采样。
    
    % 作图
    figure;
    imagesc(abs(s_ac_test));
    title('将成像结果做升采样，看效果如何');
    
    %% 
    % 下面分别对点目标中心（二维最大值）做行切片，和列切片。
    % 每一个切片，都做 16倍 升采样。
    % 并分别作出归一化的对数幅度图。
    % -------------------------------------------------------
    % 考虑图像扭曲
    % 将目标进行旋转和扭曲，以使大部分旁瓣对齐至水平轴和垂直轴
    % -------------------------------------------------------
    
    %% 
    % 取出行切片
    
    % 第一步，找出升采样后点目标中心的位置和大小
    [row_test,column_test] = size(s_ac_test);      % s_ac_test的矩阵大小
    [aa_test,p_test] = max(abs(s_ac_test));           
    [bb_test,q_test] = max(max(abs(s_ac_test)));
    
    row_test_max = p_test(q_test);      % 二维矩阵最大值，所在第几行—— row_test_max。
    column_test_max = q_test;           % 二维矩阵最大值，所在第几列—— column_test_max。
    s_ac_test_max = bb_test;            % 矩阵最大值是—— x_max。
    
    % 第二步，找出点目标中心左侧2/3长度内的最大值，以此来计算距离旁瓣的扭曲角度
    % 注意到，如果直接在左侧2/3长度内寻找最大值，有可能会误判到方位旁瓣上的点；
    % 因此，由于斜视角的关系，左侧的距离旁瓣一定是向上倾斜的；
    % 因此下面的寻找中，行直接规定为在点目标中心所在行的上方。
    [aa_test_2,p_test_2] = max(abs(s_ac_test(1:row_test_max,1:2*column_test_max/3)));           
    [bb_test_2,q_test_2] = max(max(abs(s_ac_test(1:row_test_max,1:2*column_test_max/3))));
    
    row_test_max_2 = p_test_2(q_test_2);  % 二维矩阵点目标中心左侧2/3长度的最大值，所在第几行
    column_test_max_2 = q_test_2;% 二维矩阵点目标中心左侧2/3长度的最大值，所在第几列
    
    % 第三步，计算距离旁瓣扭曲的转角
    range_theta = atan( abs((row_test_max_2-row_test_max)/(column_test_max_2-column_test_max)) );
    % 结果是弧度，下面转成角度
    range_theta = range_theta/pi*180;   % 这是距离旁瓣扭曲的角度
    
    % 第四步，将升采样的结果s_ac_test以角range_theta进行旋转（这里是逆时针旋转）
    s_ac_range = imrotate(s_ac_test,range_theta,'bilinear');% 采用'bilinear'，双线性插值的方式
    % s_ac_range是旋转后的图像
    % 距离向切片就以此来计算。
    
    % 作图
    figure;
    imagesc(abs(s_ac_range));
    title('将距离向旁瓣旋转到平行于水平轴后的成像结果');
    
    % 第五步，找出旋转后的最大值中心，并取出相应的行切片
    [aa_test_range,p_test_range] = max(abs(s_ac_range));           
    [bb_test_range,q_test_range] = max(max(abs(s_ac_range)));
    row_test_max_range = p_test_range(q_test_range); % 旋转后，点目标中心，所在第几行
    column_test_max_range = q_test_range;  % 旋转后，点目标中心，所在第几列
    
    s_ac_test_row_max = s_ac_range(row_test_max_range,column_test_max_range/3:5*column_test_max_range/3);
    % s_ac_test_row_max是取出的点目标中心行切片。
    % 其中，这里并没有把s_ac_range的一行全部取出来，而是将最左右两侧的一部分去除了。
    
    % 下面进行16倍升采样
    S_AC_test_row_max_1 = fft(s_ac_test_row_max);   % 变换到频域
    
    % -----------------------------------------------------------------------
    % 以下分别是两种补零方式：
    % 方式 1 ：在高频（中间）补零
     S_AC_test_row_max = [S_AC_test_row_max_1(1,1:length(S_AC_test_row_max_1)/2),...
                            zeros(1, 15*length(S_AC_test_row_max_1)),...
                       S_AC_test_row_max_1(1,length(S_AC_test_row_max_1)/2+1:end)];
    
    % 方式 2 ：在最小值处补15倍的零
    % S_AC_test_row_max = zeros(1,16*length(S_AC_test_row_max_1));
    % [C1,I1] = min(S_AC_test_row_max_1);
    % S_AC_test_row_max(1,1:I1) = S_AC_test_row_max_1(1,1:I1);
    % S_AC_test_row_max(1,16*length(S_AC_test_row_max_1)+1-(length(S_AC_test_row_max_1)-I1):end) = S_AC_test_row_max_1(1,I1+1:end);
    % ------------------------------------------------------------------------- 
    % 到这里为止，频域补零完成    
    s_ac_test_row_max = ifft(S_AC_test_row_max);  % 进行IFFT，回到时域，得到升采样过后的行切片。
    % 到这里为止，升采样完成
    
    % 下面做归一化的对数幅度图
    s_ac_test_row_max_abs = abs(s_ac_test_row_max);             % 取幅度
    s_ac_test_row_max_abs = 20*log10(s_ac_test_row_max_abs);    % 取对数
    s_ac_test_row_max_abs = s_ac_test_row_max_abs - max(s_ac_test_row_max_abs); % 归一化
    
    % 作图
    figure;
    plot(s_ac_test_row_max_abs);
    title('点目标中心，行切片');
    
    
    %%
    % 取出列切片
    
    % 第一步，找出升采样后点目标中心的位置和大小
    % （这和取出行切片时相同）
    
    % 第二步，找出点目标中心上方（方位向）2/3长度内的最大值，以此来计算方位旁瓣的扭曲角度
    % 注意到，如果直接在上方2/3长度内寻找最大值，由于距离旁瓣较强，所以有可能会误判到距离旁瓣上的点；
    % 因此，由于斜视角的关系，上方的方位旁瓣一定是向右倾斜的；
    % 因此下面的寻找中，列直接规定为在点目标中心所在列的右侧。
    [aa_test_3,p_test_3] = max(abs(s_ac_test(1:2*row_test_max/3,column_test_max:end-10)));           
    [bb_test_3,q_test_3] = max(max(abs(s_ac_test(1:2*row_test_max/3,column_test_max:end-10))));
    
    row_test_max_3 = p_test_3(q_test_3);  % 二维矩阵点目标中心上方2/3长度的最大值，所在第几行
    column_test_max_3 = q_test_3;% 二维矩阵点目标中心上方2/3长度的最大值，所在第几列，这里要注意！！
    % 注意：column_test_max_3 由于是直接在点目标中心所在列的右侧进行判断的，因此这里的值是相对于点目标中心所在列的相对值。
    
    % 第三步，计算方位旁瓣扭曲的转角
    azimuth_theta = atan( abs((column_test_max_3)/(row_test_max_3-row_test_max)) );
    % 结果是弧度，下面转成角度
    azimuth_theta = azimuth_theta/pi*180;   % 这是方位旁瓣扭曲的角度
    
    % 第四步，将升采样的结果s_ac_test以角azimitu_theta进行旋转（这里是逆时针旋转）
    s_ac_azimuth = imrotate(s_ac_test,azimuth_theta,'bilinear');% 采用'bilinear'，双线性插值的方式
    % s_ac_azimuth是旋转后的图像
    % 方位向切片就以此来计算。
    
    % 作图
    figure;
    imagesc(abs(s_ac_azimuth));
    title('将方位向旁瓣旋转到平行于垂直轴后的成像结果');
    
    % 第五步，找出旋转后的最大值中心，并取出相应的列切片
    [aa_test_azimuth,p_test_azimuth] = max(abs(s_ac_azimuth));           
    [bb_test_azimuth,q_test_azimuth] = max(max(abs(s_ac_azimuth)));
    row_test_max_azimuth = p_test_azimuth(q_test_azimuth); % 旋转后，点目标中心，所在第几行
    column_test_max_azimuth = q_test_azimuth;  % 旋转后，点目标中心，所在第几列。
    
    s_ac_test_column_max = s_ac_azimuth(row_test_max_azimuth/3:5*row_test_max_azimuth/3,column_test_max_azimuth);
    % s_ac_test_column_max是取出的点目标中心列切片。
    % 其中，这里并没有把s_ac_azimuth的一列全部取出来，而是将最上下两侧的一部分去除了。
    s_ac_test_column_max = s_ac_test_column_max.';  % 转为行向量，便于统一处理
    
    % 下面进行16倍升采样
    S_AC_test_column_max_1 = fft(s_ac_test_column_max);
    % 在最小值处补15倍的零
    S_AC_test_column_max = zeros(1,16*length(S_AC_test_column_max_1));
    [C2,I2] = min(S_AC_test_column_max_1);
    S_AC_test_column_max(1,1:I2) = S_AC_test_column_max_1(1,1:I2);
    S_AC_test_column_max(1,16*length(S_AC_test_column_max_1)+1-(length(S_AC_test_column_max_1)-I2):end) = S_AC_test_column_max_1(1,I2+1:end);
    s_ac_test_column_max = ifft(S_AC_test_column_max);    	% 升采样过后的列切片。
    % 完成16倍升采样   
    
    % 下面做归一化的对数幅度图
    s_ac_test_column_max_abs = abs(s_ac_test_column_max);           % 取幅度   
    s_ac_test_column_max_abs = 20*log10(s_ac_test_column_max_abs);  % 取对数
    s_ac_test_column_max_abs = s_ac_test_column_max_abs - max(s_ac_test_column_max_abs);% 归一化
    
    % 作图
    figure;
    plot(s_ac_test_column_max_abs);
    title('点目标中心，列切片');
    
    %%
    % 行切片，点目标中心距离向指标
    [PSLR_r,ISLR_r,IRW_r] = zhibiao_2(s_ac_test_row_max,16*8*NN,NN/Fr);
    parameter_r = [PSLR_r,ISLR_r,IRW_r];
    
    % 列切片，点目标中心方位向指标
    [PSLR_a,ISLR_a,IRW_a] = zhibiao_2(s_ac_test_column_max,16*8*NN,NN/Fa);
    IRW_a = (IRW_a/c*2)*Vr;
    parameter_a = [PSLR_a,ISLR_a,IRW_a];
    
    disp('------------------------------------------------------');
    disp('行切片，点目标中心距离向指标');
    disp('      PSLR     ISLR       IRW');
    disp(parameter_r);
    disp('------------------------------------------------------');
    disp('列切片，点目标中心方位向指标');
    disp('      PSLR     ISLR       IRW');
    disp(parameter_a);
    disp('------------------------------------------------------');
    
    