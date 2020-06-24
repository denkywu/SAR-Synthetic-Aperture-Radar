function  [s_imag_B_after_CoRe,R] = co_registration(s_imag_A,s_imag_B)
% 本函数用来对生成的两幅 SLC 进行“图像配准”
% 以天线 A 的SLC-A作为参考，将天线 B 的SLC-B（副图像）进行“配准”；
% 实现方法：
%   分为两步：粗配准和精配准；
%   1）粗配准，采用最大相关系数法：
%      分别计算SLC-B沿不同方位和距离平移后，与SLC-A的相关系数，当相关系数最大
%      时，即为粗配准需要采用的平移量；
%   2）精配准——（未完成）
%
% 输入变量：
% 1）s_imag_A是天线 A 的SLC-A；
% 2）s_imag_B 是天线 B 的 SLC-B；
% 输出变量：
% 1）s_imag_B_after_CoRegistration 配准后的SLC-B图像；
% 2）R 是粗配准时计算得到的不同平移量时的相关系数；
%
% 该程序截止到：2014.12.22. 17:05 p.m.

%%
disp('正在进行“图像配准”，请等待')
%*************************************************************************
%                        下面首先进行“粗配准”
%*************************************************************************

% -----------------------------------------------------------------------
%           计算两幅图像在不同方位和距离偏移处的互相关系数 R
% -----------------------------------------------------------------------
% 利用幅度计算，得到幅度互相关为：
for pp = 1:5            % pp 代表沿方位方向的偏移
    for qq = 1: 5       % qq 代表沿斜距方向的偏移
        tmp = circshift(s_imag_B,[pp-1,qq-1]);
        R(pp,qq) = sum(sum(abs(s_imag_A).*abs(tmp)))/...
            ( sqrt(sum(sum(abs(s_imag_A).^2)))*sqrt(sum(sum(abs(tmp).^2))) );
        clear tmp;
    end
end
[cor_R1,cor_p] = max(R);
[cor_R2,cor_q] = max(max(R));
% 注意：
%   1）我们得到的（cor_p(cor_q)，cor_q）是相关系数 R 最大值的位置；
%   2）而粗配准需要移动的大小与此密切相关；
%   3）需要移动的大小是：（cor_p(cor_q)-1，cor_q-1）

% -----------------------------------------------------------------------
%                       对天线 B 的成像结果进行粗配准
% -----------------------------------------------------------------------
% s_imag_B_2 就是经过粗配准后的结果，如下：
s_imag_B_2 = circshift(s_imag_B,[cor_p(cor_q)-1,cor_q-1]); 

s_imag_B_after_CoRe = s_imag_B_2;
%%
%*************************************************************************
%                               “精配准”
%*************************************************************************
% 下面在粗配准的基础上进行“精配准”




%{
% 利用 sinc 插值，将第二幅图像进行“整体配准”
%       ——类似进行RCMC，只是这里的整幅图像所需的平移量都是相同的，为 delta_N；
R = 16;              % sinc插值核长度
[num_azi,num_rag] = size(s_imag_B);
s_imag_B_after_CoRegistration = zeros(num_azi,num_rag); % 用来存放整体配准后的值
disp('正在进行插值计算，请等待');
h = waitbar(0,'正在进行插值计算，请等待');
for pp = 1 : num_azi
    for qq = 1 : num_rag
        N_pp_qq = qq + delta_N;  
        N_pp_qq = rem(N_pp_qq,num_rag);   
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
        rcmc_r_imag_B = s_imag_B(pp,ll);
        s_imag_B_after_CoRegistration(pp,qq) = sum( registration_sinc.*rcmc_r_imag_B );
    end
    waitbar(pp/num_azi);
end
close(h);
% s_imag_B_after_CoRegistration 即是 s_imag_B 经过整体配准后的结果
%}


disp('完成对天线 B 的 SLC-B 的“图像配准”');

end