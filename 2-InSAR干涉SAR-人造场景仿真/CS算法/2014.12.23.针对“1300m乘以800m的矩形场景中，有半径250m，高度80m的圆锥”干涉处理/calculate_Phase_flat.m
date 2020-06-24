function phase_flat_earth = calculate_Phase_flat(R0_RCMC,Parameter,B,theta_B)
% 该函数用来计算平地相位
% 输入量：
% 1）R0_RCMC     是随距离线变化的斜距；
% 2）Parameter   代表计算平地相位需要的一些参数，包括：Naz，H，lamda；
% 3）B           代表基线距离；
% 4）theta_B     代表基线倾角；
% 返回值：phase_flat_earth 是所需要的平地相位。
%
% 程序截止到：2014.12.17. 19:57 p.m.

%%
% ------------------------------------------------------------------------
%                   载入所需要的一些参数
% ------------------------------------------------------------------------
Naz = Parameter(1,1);       % Parameter 的第一行代表 Naz
H = Parameter(2,1);         % Parameter 的第二行代表 H
lamda = Parameter(3,1);     % Parameter 的第三行代表 lamda　

%%
% ------------------------------------------------------------------------
%                               计算平地相位
% ------------------------------------------------------------------------
R0_RCMC_mtx = ones(Naz,1)*R0_RCMC; % 形成矩阵
r_range_A = R0_RCMC_mtx;         % 天线 A 到平地场景的最近斜距

% 下面计算对应于天线 A 的最近斜距时，天线 B 所对应的最近斜距
cos_theta = H./r_range_A;         % 对应于每一个天线A的最近斜距时，波束下视角的余弦；
theta_range_A = acos(cos_theta);  % 计算出波束下视角；  

% 天线 B 到平地场景的最近斜距，如下：
r_range_B = sqrt(r_range_A.^2 + B^2 + 2*B.*r_range_A.*cos(pi/2-theta_B + theta_range_A)); 

phase_flat_earth = -4*pi/lamda.*(r_range_A-r_range_B);  % 计算得到的平地相位;
                                                        % 注意，负号很关键！！这由SAR的回波特点决定

end
