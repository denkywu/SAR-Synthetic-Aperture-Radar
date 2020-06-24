function PHY_after_unwrapping = LS_unwrapping(PHY_s_after_X_filtering)
% 本函数用于相位解缠绕，采用 “最小二乘方法”；
% 残差点个数为 0 时，或者不等于 0 时都可以使用；
%           ——参考保铮等著的《雷达成像技术》第8.6.3节，详见第311页到第315页；
%
% 输入数据：
%   PHY_s_after_X_filtering  表示经过某种滤波方法后的干涉相位图；
% 输出数据：
%   PHY_after_unwrapping  表示相位解缠绕后的相位图；
%
% 程序截止至： 2014.12.22. 15:21 a.m.

%%
disp('正在采用“最小二乘方法”进行相位解缠绕，请等待');
% ------------------------------------------------------------------------
%                               相位解缠绕
% ------------------------------------------------------------------------
% 步骤 1 ：
% 对“PHY_s_after_X_filtering”作镜像反射。
% 在二维平面，将“PHY_s_after_X_filtering”先以 m = Naz+1/2 为轴作镜像反射；
%                                         再以 n = Nrg+1/2 为轴作镜像反射；
% 得到所谓“半采样对称图像”；
[Naz,Nrg] = size(PHY_s_after_X_filtering);  % Naz*Nrg 大小的矩阵
PHY = zeros(2*Naz,2*Nrg);

PHY(1:Naz,1:Nrg) = PHY_s_after_X_filtering;
PHY(Naz+1:end,1:Nrg) = flipud(PHY_s_after_X_filtering);
PHY(1:Naz,Nrg+1:end) = fliplr(PHY_s_after_X_filtering);
PHY(Naz+1:end,Nrg+1:end) = fliplr(flipud(PHY_s_after_X_filtering));

% 步骤 2 ：
% 分别沿两个方向求出一阶差分；
% 下面首先求沿方位向的差分：
delta_x = zeros(2*Naz,2*Nrg);
delta_x = circshift(PHY,[-1 0]) - PHY;
% 接下来条件判断
L_delta_x_pi = delta_x > pi;
delta_x(L_delta_x_pi) = delta_x(L_delta_x_pi) - 2*pi;
L_delta_x_minus_pi = delta_x < -pi;
delta_x(L_delta_x_minus_pi) = delta_x(L_delta_x_minus_pi) + 2*pi;

% 下面求沿斜距向的差分
delta_y = zeros(2*Naz,2*Nrg);
delta_y = circshift(PHY,[0 -1]) - PHY;
% 接下来条件判断
L_delta_y_pi = delta_y > pi;
delta_y(L_delta_y_pi) = delta_y(L_delta_y_pi) - 2*pi;
L_delta_y_minus_pi = delta_y < -pi;
delta_y(L_delta_y_minus_pi) = delta_y(L_delta_y_minus_pi) + 2*pi;

% 步骤 3 ：
% 根据式（8.57）求出 ruo ；
ruo = zeros(2*Naz,2*Nrg);
ruo = (delta_x - circshift(delta_x,[1 0 ])) + (delta_y - circshift(delta_y,[0 1]));
% 由此得到 ruo，矩阵大小 2Naz * 2Nrg

% 步骤 4 ：
% 对 ruo 作二维FFT，得到 P；
P = fft2(ruo);

% 步骤 5 ：
% 用式（8.60）求出 PHY_P 
k = 1:2*Naz;
% P1 = cos(pi*k/Naz);
P1 = cos(pi*(k-1)/Naz);
P1_mtx = (P1.')*ones(1,2*Nrg);

l = 1:2*Nrg;
% P2 = cos(pi*l/Nrg);
P2 = cos(pi*(l-1)/Nrg);
P2_mtx = ones(2*Naz,1)*P2;

PHY_P = P./(2*P1_mtx + 2*P2_mtx -4);
% PHY_P(2*Naz,2*Nrg) = 0;     % 因为分母为零，所以此项的公式计算结果为NaN，故这里重新赋值；
PHY_P(1,1) = 0;

% 步骤 6 ：
% 对 PHY_P 作二维 IFFT；
phy_p = ifft2(PHY_P);

% 步骤 7 ：
% 取出 phy_p 中原始区间部分的值，即得到相位解缠绕的结果
PHY_after_unwrapping = phy_p(1:Naz,1:Nrg);

disp('完成二维相位解缠绕');


end