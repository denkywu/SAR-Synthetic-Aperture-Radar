function [PSLR_r,ISLR_r,IRW_r, PSLR_a,ISLR_a,IRW_a] = target_analysis(s_ac,Fr,Fa,Vr)
% 补零方式很关键
% 这里采用的方法是逐行，逐列判断数据最小值位置，然后将在最小值位置处进行补零

% 输入变量： s_ac  ——是需要进行指标分析的点目标点k，s_ac是它所在的一块矩阵数据。
% 输入变量： Fr  —— 距离向采样率
% 输入变量： Fa  —— 方位向采样率
% 输入变量： Vr  —— 平台速度
% 该函数用来求解点目标k的各项指标
% 点目标中心的行切片（距离向）：峰值旁瓣比(PSLR)，积分旁瓣比(ISLR),距离分辨率（IRW）。
% 点目标中心的列切片（方位向）：峰值旁瓣比(PSLR)，积分旁瓣比(ISLR),距离分辨率（IRW）。
% 输出值分别是行切片和列切片的各项指标


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       取点目标中心 NN*NN 切片，进行升采样
%   再取点目标中心，距离向切片，方位向切片，进行升采样
%                   计算各项性能指标
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
% 取点目标的中心切片，NN*NN
% 进行二维升采样——如何进行二维升采样？？？
% 高频补零！！！！！——因此，判断哪里是高频，很重要！！
c = 3e8;                % 光速

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
S_ac_test_buling_1 = zeros(NN,8*NN);  % 中间变量
S_ac_test_buling = zeros(8*NN,8*NN);
% ====================================================================
% 利用两个for循环，在每行和每列最小值的位置处补上7*NN的零，实现高频补零。
for pp = 1:NN           % 在每行的最小值位置补零
    [C,I] = min(S_ac_test_2(pp,:));
    S_ac_test_buling_1(pp,1:I) = S_ac_test_2(pp,1:I);
    S_ac_test_buling_1(pp,8*NN-(NN-I)+1:8*NN) = S_ac_test_2(pp,I+1:NN);
end
for qq = 1:8*NN         % 在每列的最小值位置补零
    [C,I] = min(S_ac_test_buling_1(:,qq));
    S_ac_test_buling(1:I,qq) = S_ac_test_buling_1(1:I,qq);
    S_ac_test_buling(8*NN-(NN-I)+1:8*NN,qq) = S_ac_test_buling_1(I+1:NN,qq);
end     
% ====================================================================
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

% 这里没有考虑点目标中心轴的旋转
% 没有将真正的最大值所在的行和列旋转到x轴和y轴方向

[row_test,column_test] = size(s_ac_test);      % s_ac_test的矩阵大小
[aa_test,p_test] = max(abs(s_ac_test));           
[bb_test,q_test] = max(max(abs(s_ac_test)));

row_test_max = p_test(q_test);      % 二维矩阵最大值，所在第几行—— row_max。
column_test_max = q_test;           % 二维矩阵最大值，所在第几列—— column_max。
s_ac_test_max = bb_test;            % 矩阵最大值是—— x_max。

% 做行切片，升采样，做归一化对数幅度图
s_ac_test_row_max = (s_ac_test(row_test_max,:)); 	% 取出最大值所在的行，做行切片。
S_AC_test_row_max_1 = fft(s_ac_test_row_max);

% 在最小值处补15倍的零
S_AC_test_row_max = zeros(1,16*length(S_AC_test_row_max_1));
[C1,I1] = min(S_AC_test_row_max_1);
S_AC_test_row_max(1,1:I1) = S_AC_test_row_max_1(1,1:I1);
S_AC_test_row_max(1,16*length(S_AC_test_row_max_1)+1-(length(S_AC_test_row_max_1)-I1):end) = S_AC_test_row_max_1(1,I1+1:end);
% 完成16倍升采样     

s_ac_test_row_max = ifft(S_AC_test_row_max);        % 升采样过后的行切片。
s_ac_test_row_max_abs = abs(s_ac_test_row_max);             % 取幅度
s_ac_test_row_max_abs = 20*log10(s_ac_test_row_max_abs);    % 取对数
s_ac_test_row_max_abs = s_ac_test_row_max_abs - max(s_ac_test_row_max_abs); % 归一化


% 做列切片，升采样，做归一化对数幅度图
s_ac_test_column_max = (s_ac_test(:,column_test_max)); 	% 取出最大值所在的列，做列切片。
s_ac_test_column_max = s_ac_test_column_max.';          % 转置，表示成行向量（便于统一处理）
S_AC_test_column_max_1 = fft(s_ac_test_column_max);
    % 对列切片做升采样时，和之前的成像切片升采样相同，要考虑真正的高频位置。
    % 因此要分两种情况，对应多普勒模糊的两种处理方法，分别采用相应的补零方法
% ====================================================================
% 对应于处理多普勒模糊的方法二~~~移动频率轴
% 高频就在fa_1和fa_2相接的那一段

% 在最小值处补15倍的零
S_AC_test_column_max = zeros(1,16*length(S_AC_test_column_max_1));
[C2,I2] = min(S_AC_test_column_max_1);
S_AC_test_column_max(1,1:I2) = S_AC_test_column_max_1(1,1:I2);
S_AC_test_column_max(1,16*length(S_AC_test_column_max_1)+1-(length(S_AC_test_column_max_1)-I2):end) = S_AC_test_column_max_1(1,I2+1:end);
% 完成16倍升采样   
% ====================================================================

s_ac_test_column_max = ifft(S_AC_test_column_max);    	% 升采样过后的列切片。
s_ac_test_column_max_abs = abs(s_ac_test_column_max);           % 取幅度   
s_ac_test_column_max_abs = 20*log10(s_ac_test_column_max_abs);  % 取对数
s_ac_test_column_max_abs = s_ac_test_column_max_abs - max(s_ac_test_column_max_abs);% 归一化

% 作图
figure;
plot(s_ac_test_row_max_abs);
title('点目标中心，行切片');
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

