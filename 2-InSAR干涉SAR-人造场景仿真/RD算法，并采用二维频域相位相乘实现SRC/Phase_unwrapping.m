function PHY_after_unwrapping = Phase_unwrapping(PHY_s_after_X_filtering)
% 本函数用于在残差点个数为 0 时进行相位解缠绕；
% ！！！！！    注意   ！！！！！
% 当残差点个数为 0 时，才可以使用本函数进行相位解缠绕；
% 输入数据：
%   PHY_s_after_X_filtering  表示经过某种滤波方法后的干涉相位图；
% 输出数据：
%   PHY_after_unwrapping  表示相位解缠绕后的相位图；
%
% 思路：
%   此时积分结果与路径无关，故可以选择积分路径如下：
%       a）先从左至右解缠绕第一行，再从上向下分别解缠绕各列；
%       b）先从上到下解缠绕第一列，再从左向右分别解缠绕各行；
%   下面可以选择用这两种积分路径中的某一种进行相位解缠绕
%  注意：
%   我经过仿真验证，这两种积分路径的解缠绕结果只有最大约为 10^(-14) 次方的误差；
%   因此，可以认为这两种积分路径下的解缠绕结果是相同的。
%   也同时说明，残差点个数为 0 时，积分结果是与积分路径无关的。
%
% 程序截止到： 2014.12.18. 16:35 p.m.

%%
% ------------------------------------------------------------------------
%                               相位解缠绕
% ------------------------------------------------------------------------
% 注意：
% 下面的两种路径计算结果应该是一致的。
% 选择其中的一种路径进行相位解缠绕的计算即可。
% 不要同时将两种路径都选上！！

[Naz,Nrg] = size(PHY_s_after_X_filtering);   % 数据大小

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 路径1： 
%   先从左至右解缠绕第一行，再从上向下分别解缠绕各列；
disp('正在进行二维相位解缠绕，请等待');

% 下面首先将第一个值赋值给解缠绕后的第一个值，以此为基础
PHY_after_unwrapping(1,1) = PHY_s_after_X_filtering(1,1);
% 下面从左至右解缠绕第一行
for qq = 2:Nrg
    delta_qq = PHY_s_after_X_filtering(1,qq) - PHY_s_after_X_filtering(1,qq-1);
    if delta_qq > pi
        delta_qq = delta_qq - 2*pi;
    end
    if delta_qq < -1*pi
        delta_qq = delta_qq +2*pi;
    end
    PHY_after_unwrapping(1,qq) = PHY_after_unwrapping(1,qq-1) + delta_qq;
end
% 下面分别从上向下分别解缠绕各列
for qq = 1:Nrg
    for pp = 2:Naz
        delta_qq_pp = PHY_s_after_X_filtering(pp,qq) - PHY_s_after_X_filtering(pp-1,qq);
        if delta_qq_pp > pi
            delta_qq_pp = delta_qq_pp - 2*pi;
        end
        if delta_qq_pp < -1*pi
            delta_qq_pp = delta_qq_pp + 2*pi;
        end 
        PHY_after_unwrapping(pp,qq) = PHY_after_unwrapping(pp-1,qq) + delta_qq_pp;
    end
end
disp('完成二维相位解缠绕')
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
% 路径2
%   先从上到下解缠绕第一列，再从左向右分别解缠绕各行；
disp('正在进行二维相位解缠绕，请等待');

% 下面首先将第一个值赋值给解缠绕后的第一个值，以此为基础
PHY_after_unwrapping(1,1) = PHY_s_after_X_filtering(1,1);
% 先从上到下解缠绕第一列
for pp = 2:Naz
    delta_pp = PHY_s_after_X_filtering(pp,1) - PHY_s_after_X_filtering(pp-1,1);
    if delta_pp > pi
        delta_pp = delta_pp - 2*pi;
    end
    if delta_pp < -1*pi
        delta_pp = delta_pp + 2*pi;
    end
    PHY_after_unwrapping(pp,1) = PHY_after_unwrapping(pp-1,1) + delta_pp;
end
% 下面分别从左向右分别解缠绕各行
for pp = 1:Naz
    for qq = 2:Nrg
        delta_pp_qq = PHY_s_after_X_filtering(pp,qq) - PHY_s_after_X_filtering(pp,qq-1);
        if delta_pp_qq > pi
            delta_pp_qq = delta_pp_qq - 2*pi;
        end
        if delta_pp_qq < -1*pi
            delta_pp_qq = delta_pp_qq + 2*pi;
        end 
        PHY_after_unwrapping(pp,qq) = PHY_after_unwrapping(pp,qq-1) + delta_pp_qq;
    end
end
disp('完成二维相位解缠绕')
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end