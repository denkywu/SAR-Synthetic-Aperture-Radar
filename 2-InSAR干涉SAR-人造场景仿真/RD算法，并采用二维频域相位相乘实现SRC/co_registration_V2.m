function s_imag_B_after_CoRe = co_registration_V2(s_imag_A,s_imag_B,R0_RCMC,Parameter,B,theta_B)
% 本函数用来对生成的两幅 SLC 进行“图像配准”
% ************************************************************************
%
%                                思想
%
% 1）根据场景设计，逐点计算每一点的相应的斜距平移量；
% 2）利用sinc插值逐点进行配准：
%   a）对图像 B 先进行斜距向16倍升采样；
%   b）然后对升采样后的结果，按照计算出来的斜距平移量，逐点用sinc插值计算；
%   c）对插值后的结果相应降采样，回到原来的矩阵大小；
%      实际过程中，b）和c）同时完成；
%   d）进行幅值恢复；
%   e）完成对图像 B 的配准；
% ************************************************************************
%
% 数据说明：
% 输入变量：
%   1）s_imag_A    是天线 A 的SLC-A；
%   2）s_imag_B    是天线 B 的 SLC-B；
%   3）R0_RCMC     是随距离线变化的斜距；
%   4）Parameter   代表函数计算需要的一些参数，包括：Naz，H，lamda；
%   5）B           代表基线距离；
%   6）theta_B     代表基线倾角；
% 输出变量：
%   1）s_imag_B_after_CoRe 配准后的SLC-B图像；
%
% 该程序截止到：2015.04.30. 9:15 a.m.

%%
disp('正在进行“图像配准”，请等待')
% -------------------------------------------------------------------------
%                           载入一些需要的参数
% -------------------------------------------------------------------------
s_A = s_imag_A;                   % 天线 A 的成像结果
s_B = s_imag_B;                   % 天线 B 的成像结果

c = 299792458;              % 光速

Naz = Parameter(1,1);       % Parameter 的第一行代表 Naz
H = Parameter(2,1);         % Parameter 的第二行代表 H
lamda = Parameter(3,1);     % Parameter 的第三行代表 lamda
Fr = Parameter(4,1);        % Parameter 的第四行代表 Fr

%%
% -------------------------------------------------------------------------
%               对天线 B 的成像结果，首先进行斜距向的16倍升采样
% -------------------------------------------------------------------------
[Naz,Nrg] = size(s_B);

S_B = fft(s_B,[],2);
S_B_BU_LING = [S_B(:,1:Nrg/2),zeros(Naz,15*Nrg),S_B(:,Nrg/2+1:end)];

s_B_UP = ifft(S_B_BU_LING,[],2);    % 完成16倍升采样，回到时域
% figure;imagesc(abs(s_B_UP));title('天线 B，升采样后的成像结果');

%%
% -------------------------------------------------------------------------
%               下面对天线 B 的升采样结果，进行逐点 sinc 插值配准
% -------------------------------------------------------------------------
%
% 首先计算移动量
% -------------------------------------------------------------------------
%                               计算平移量
% -------------------------------------------------------------------------
R0_RCMC_mtx = ones(Naz,1)*R0_RCMC; % 形成矩阵
r_range_A = R0_RCMC_mtx;         % 天线 A 到平地场景的最近斜距

% 下面计算对应于天线 A 的最近斜距时，天线 B 所对应的最近斜距
cos_theta = H./r_range_A;         % 对应于每一个天线A的最近斜距时，波束下视角的余弦；
theta_range_A = acos(cos_theta);  % 计算出波束下视角；  

% 天线 B 到平地场景的最近斜距，如下：
r_range_B = sqrt(r_range_A.^2 + B^2 + 2*B.*r_range_A.*cos(pi/2-theta_B + theta_range_A)); 

% 下面逐点计算最近斜距差
delta_r = r_range_A - r_range_B; 

% 计算“图像配准”需要沿斜距向移动的点数
% P.S 因为仿真条件下，方位向是完全配准的，因此只需要考虑沿斜距向的配准
delta_N = 2*delta_r/c*Fr; 
delta_N = -delta_N;     	% 需要沿斜距向移动的点数。
                        	%  P.S. “负号”很重要，代表是将副图像向右移位。
           
% 然后进行图像配准（使用sinc插值）
% -----------------------------------------------------------------------
%                    对天线 B 的成像结果进行图像配准
% -----------------------------------------------------------------------
% 利用 sinc 插值
R = 16;              % sinc插值核长度
[num_azi,num_rag] = size(s_B_UP);
h = waitbar(0,'正在进行插值处理，请等待');
for pp = 1 : num_azi
    kk = 1;
    for qq = 1:16: num_rag
        N_pp_qq = qq + 16*delta_N(pp,kk);  
        N_pp_qq_zheng = ceil(N_pp_qq);        % ceil，向上取整。
        ii = ( N_pp_qq-(N_pp_qq_zheng-R/2):-1:N_pp_qq-(N_pp_qq_zheng+R/2-1)  );
        registration_sinc = sinc(ii);
        registration_sinc = registration_sinc/sum(registration_sinc);   % 插值核的归一化
        % ii 是sinc插值过程的变量;
        % g(x)=sum(h(ii)*g_d(x-ii)) = sum(h(ii)*g_d(ll));
        % 由于s_imag_B只有整数点取值，且范围有限。因此插值中要考虑它的取值溢出边界问题。
        % 这里我采取循环移位的思想，用来解决取值溢出问题。
        if (N_pp_qq_zheng-R/2) > num_rag    % 全右溢
            ll = (N_pp_qq_zheng-R/2-num_rag:1:N_pp_qq_zheng+R/2-1-num_rag);
        else
            if (N_pp_qq_zheng+R/2-1) > num_rag    % 部分右溢
                ll_1 = (N_pp_qq_zheng-R/2:1:num_rag);
                ll_2 = (1:1:N_pp_qq_zheng+R/2-1-num_rag);
                ll = [ll_1,ll_2];
            else
                if (N_pp_qq_zheng+R/2-1) < 1    % 全左溢（不可能发生，但还是要考虑）
                    ll = (N_pp_qq_zheng-R/2+num_rag:1:N_pp_qq_zheng+R/2-1+num_rag);
                else
                    if (N_pp_qq_zheng-R/2) < 1       % 部分左溢
                        ll_1 = (N_pp_qq_zheng-R/2+num_rag:1:num_rag);
                        ll_2 = (1:1:N_pp_qq_zheng+R/2-1);
                        ll = [ll_1,ll_2];
                    else
                        ll = (N_pp_qq_zheng-R/2:1:N_pp_qq_zheng+R/2-1);
                    end
                end
            end
        end
        rcmc_r_imag_B = s_B_UP(pp,ll);
        s_B_ReCo(pp,kk) = sum( registration_sinc.*rcmc_r_imag_B );
        kk = kk + 1;
    end
    waitbar(pp/num_azi);
end
close(h);
% 得到的 s_B_ReCo 就已经是经过升采样，然后插值配准，再降采样的数据了；
% 下面再恢复幅值：
s_B_ReCo = 16.*s_B_ReCo;            % 恢复幅值

% 此时的 s_B_ReCo 就是我们需要的最终结果：
% 1）经过升采样；
% 2）然后进行了采用sinc插值的配准；
% 3）降采样回到原来的大小；
% 4）恢复幅值

s_imag_B_after_CoRe = s_B_ReCo;     % 返回值 

disp('完成对天线 B 的 SLC-B 的“图像配准”');

end










