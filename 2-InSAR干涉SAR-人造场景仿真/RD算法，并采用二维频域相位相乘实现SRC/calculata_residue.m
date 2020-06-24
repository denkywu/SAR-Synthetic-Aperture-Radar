function [PHY_residue,residue_count] = calculata_residue(PHY_s_after_X_filtering)
% 本函数用来计算经过前面步骤后所得到的干涉相位图的残差点；
% 输入数据：
%   PHY_s_after_X_filtering  表示经过某种滤波方法后的干涉相位图；
% 输出数据：
%   1）PHY_residue  表示计算得到的残差点；
%   2）residue_count 表示残差点个数（包括正负残差点）；
%
% 思路：
%   以输入的干涉相位图为基础，（逐点进行）以（m,n）为第一个值的2*2窗口，共四个点的
%   环路积分计算。计算结果进行以下标注：
%   1）如果顺时针环路积分计算结果为2π，则为正残点，将(m,n)标记为 +1；
%   2）如果顺时针环路积分计算结果为-2π，则为负残点，将(m,n)标记为 -1；
%   3）如果环路积分为0，则说明没有残点，将(m,n)标记为 0；
%   全部计算完成后，返回 PHY_residue，（这里面记录了所有标记结果）以此作为残差点
%   计算结果。
%
% 程序截止到： 2014.12.18. 15:16 p.m.

%%
% ------------------------------------------------------------------------
%                               计算残差点
% ------------------------------------------------------------------------
% 思路：
%   以输入的干涉相位图为基础，（逐点进行）以（m,n）为中心的2*2窗口，共四个点的
%   环路积分计算。计算结果进行以下标注：
%   1）如果顺时针环路积分计算结果为2π，则为正残点，将(m,n)标记为 +1；
%   2）如果顺时针环路积分计算结果为-2π，则为负残点，将(m,n)标记为 -1；
%   3）如果环路积分为0，则说明没有残点，将(m,n)标记为 0；
%   全部计算完成后，返回 PHY_residue，（这里面记录了所有标记结果）以此作为残差点
%   计算结果。

[Naz,Nrg] = size(PHY_s_after_X_filtering);   % 数据大小
PHY_residue = zeros(Naz,Nrg);   % 赋 0 初始值
residue_count = 0;              % 对残差点（包括正负残差点）计数

disp('正在进行残差点计算，请等待');
h = waitbar(0,'正在进行残差点计算，请等待');
for pp = 1:Naz-1
    for qq = 1:Nrg-1
        delta_1 = PHY_s_after_X_filtering(pp,qq+1) - PHY_s_after_X_filtering(pp,qq);
        if delta_1 > pi
            delta_1 = delta_1 - 2*pi;
        end
        if delta_1 < -1*pi
            delta_1 = delta_1 + 2*pi;
        end
        
        delta_2 = PHY_s_after_X_filtering(pp+1,qq+1) - PHY_s_after_X_filtering(pp,qq+1);
        if delta_2 > pi
            delta_2 = delta_2 - 2*pi;
        end
        if delta_2 < -1*pi
            delta_2 = delta_2 + 2*pi;
        end
        
        delta_3 = PHY_s_after_X_filtering(pp+1,qq) - PHY_s_after_X_filtering(pp+1,qq+1);
        if delta_3 > pi
            delta_3 = delta_3 - 2*pi;
        end
        if delta_3 < -1*pi
            delta_3 = delta_3 + 2*pi;
        end
        
        delta_4 = PHY_s_after_X_filtering(pp,qq) - PHY_s_after_X_filtering(pp+1,qq);
        if delta_4 > pi
            delta_4 = delta_4 - 2*pi;
        end
        if delta_4 < -1*pi
            delta_4 = delta_4 + 2*pi;
        end
        
        delta = delta_1 + delta_2 + delta_3 + delta_4;
        if delta == 2*pi
            PHY_residue(pp,qq) = 1;
            residue_count = residue_count + 1;
        end
        if delta == -2*pi
            PHY_residue(pp,qq) = -1;
            residue_count = residue_count + 1;
        end
    end
    waitbar(pp/(Naz-1));
end
close(h);

disp('完成残差点计算');


end