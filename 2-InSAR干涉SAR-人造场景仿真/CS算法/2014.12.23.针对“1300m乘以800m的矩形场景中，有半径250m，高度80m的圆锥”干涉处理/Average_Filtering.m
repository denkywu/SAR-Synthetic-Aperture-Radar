function PHY_s_after_avg_filtering = Average_Filtering(PHY_s_after_flat_earth,window_M,window_N)
% 本函数原理基于：    “回转均值滤波”
%
% 本函数用来对去平地相位后的干涉相位做“回转均值滤波”
% 该方法参考：保铮《雷达成像技术》第8章，8.5.2节。
%            详见第 304 页 相关内容
%
% 函数输入值：
% 1） PHY_s_after_flat_earth      代表去平地相位后的干涉相位（即缠绕相位）；
% 2） window_M 和 window_N        表示“回转均值滤波”所选定的窗口大小；
% 函数返回值：
% 1）PHY_s_after_avg_filtering 	 表示“回转均值滤波”后的干涉相位（也是缠绕相位）；
%
% 本程序截止至：2014.12.18. 14:04 p.m.

%%
% ------------------------------------------------------------------------
%                               回转均值滤波
% ------------------------------------------------------------------------
M = window_M;      % 斜距向平滑窗口大小，由外界输入
N = window_N;      % 方位向平滑窗口大小，由外界输入

[Naz,Nrg] = size(PHY_s_after_flat_earth);   % 数据大小

disp('正在进行“回转均值滤波”，请等待');
h = waitbar(0,'正在进行“回转均值滤波”，请等待');
for pp = 1:Naz
    for qq = 1:Nrg
        % 首先进行条件判断，看窗口window是否超过了矩阵的边界：
        if pp<(N+1) || pp>(Naz-N) || qq<(M+1) || qq>(Nrg-M)
            % 若满足上述条件中的任何一个，说明窗口位于矩阵边界，进行以下进一步判断
            if (pp-N) < 1
                x_min = 1;
            else
                x_min = pp - N;
            end
            if (pp+N) > Naz
                x_max = Naz;
            else
                x_max = pp + N;
            end
            if (qq-M) < 1
                y_min = 1;
            else
                y_min = qq - M;
            end
            if (qq+M) > Nrg
                y_max = Nrg;
            else
                y_max = qq + M;
            end
            PHY_window = PHY_s_after_flat_earth(x_min:x_max,y_min:y_max);
        else
            % 若上述四个条件都不满足，说明窗口不位于矩阵边界，则可以取到全部
            % （2N+1）*（2M+1）个点，因此直接用以下命令即可
            PHY_window = PHY_s_after_flat_earth(pp-N:pp+N,qq-M:qq+M);
        end
        
        % 下面根据“回转均值滤波”的方法进行处理：
        f_pp_qq_window = sum(sum(exp(1j*PHY_window)));
        PHY_s_after_avg_filtering(pp,qq) = angle(f_pp_qq_window) + ...
            1/((2*M+1)*(2*N+1))*sum(sum(angle(exp(1j*PHY_window)./f_pp_qq_window)));
    end
    waitbar(pp/Naz);
end
close(h);

disp('“回转均值滤波”已完成');


end