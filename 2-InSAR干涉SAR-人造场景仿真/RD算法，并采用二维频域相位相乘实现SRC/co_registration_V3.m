function  s_imag_B_after_CoRe = co_registration_V3(s_imag_A,s_imag_B)
% 本函数用来对生成的两幅 SLC 进行“图像配准”
% 以天线 A 的SLC-A作为参考，将天线 B 的SLC-B（副图像）进行“配准”；
% 实现方法：
%   分为两步：粗配准和精配准；
%   1）粗配准，采用最大相关系数法：
%      分别计算SLC-B沿不同方位和距离平移后，与SLC-A的相关系数，当相关系数最大
%      时，即为粗配准需要采用的平移量；
%   2）精配准
%      同样采取最大相关系数法。
%      在粗配准后，对粗配准的结果首先进行二维的升采样（20倍升采样，使得精度达到
%       0.05像素）。然后以升采样后的SLC-A为参照，对升采样后的SLC-B逐点进行相应
%      的平移，并以相关系数最大为标准得到精配准的结果。
%      详细说明见我的笔记。
%
% 输入变量：
% 1）s_imag_A 是天线 A 的SLC-A；
% 2）s_imag_B 是天线 B 的 SLC-B；
% 输出变量：
% 1）s_imag_B_after_CoRegistration 是配准后的SLC-B图像；
%
% 该程序截止到：2015.01.05. 14:30 p.m.

%%
disp('正在进行“图像配准”，请等待');

disp('首先进行“粗配准”，请等待');
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

disp('粗配准，已完成');

%%
disp('下面进行“精配准”，请等待');
%*************************************************************************
%                               “精配准”
%*************************************************************************
% 下面在粗配准的基础上进行“精配准”
% 选择20倍升采样；
N_UP = 20;

[Naz,Nrg] = size(s_imag_A);
% 下面首先将 s_imag_A 和 s_imag_B_2 分为四块；
% 然后调用函数 co_registration_V3_need（） 将分块结果分别进行精配准；
% 最后再合并即得到我们的最终结果。
% 分块 1：
s_B_new_1 = co_registration_V3_need(s_imag_A(1:Naz/2+1,1:Nrg/2+1),s_imag_B_2(1:Naz/2+1,1:Nrg/2+1),N_UP);
save('s_B_new_1.mat','s_B_new_1');
% 分块 2：
s_B_new_2 = co_registration_V3_need(s_imag_A(Naz/2:end,1:Nrg/2+1),s_imag_B_2(Naz/2:end,1:Nrg/2+1),N_UP);
save('s_B_new_2.mat','s_B_new_2');
% 分块 3：
s_B_new_3 = co_registration_V3_need(s_imag_A(1:Naz/2+1,Nrg/2:end),s_imag_B_2(1:Naz/2+1,Nrg/2:end),N_UP);
save('s_B_new_3.mat','s_B_new_3');
% 分块 4：
s_B_new_4 = co_registration_V3_need(s_imag_A(Naz/2:end,Nrg/2:end),s_imag_B_2(Naz/2:end,Nrg/2:end),N_UP); 
save('s_B_new_4.mat','s_B_new_4');

% 对上面的分块结果进行合并，得到图像 B 最终的配准结果：
s_imag_B_after_CoRe = [ s_B_new_1(1:Naz/2,1:Nrg/2) , s_B_new_3(1:Naz/2,2:end);
                        s_B_new_2(2:end,1:Nrg/2)   , s_B_new_4(2:end,2:end) ];
% s_imag_B_after_CoRe 就是 SLC-B 的最终的配准结果。

disp('精配准，已完成');
disp('完成对天线 B 的 SLC-B 的“图像配准”');

end