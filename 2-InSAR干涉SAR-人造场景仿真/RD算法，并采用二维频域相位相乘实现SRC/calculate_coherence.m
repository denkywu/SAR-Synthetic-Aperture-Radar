function Coherence_imag = calculate_coherence(s_imag_A,s_imag_B)
% 本函数用来计算两幅图像逐点的相关性；
% 方法：
% 1）只采用幅度信息计算；因此得到的是幅度相关系数图；
% 2）采用3*3的窗口计算每一点的互相关系数；
% 输入数据：
% 1）s_imag_A 输入的图像 A；
% 2）s_imag_B 输入的图像 B；
% 输出数据：
% Cohorence_imag 是返回的幅度相关系数图
%
% 该程序截止到 2015.01.04. 22:15 p.m.

%%
% ------------------------------------------------------------------------
%                           计算“幅度相关系数图”
% ------------------------------------------------------------------------
% 计算相关系数图选择的窗口大小：（2N+1）*（2M+1），即为：3*3，如下：
N = 1;
M = 1;

[Naz,Nrg] = size(s_imag_A);   % 数据大小

disp('正在计算“幅度相关系数图”，请等待');
h = waitbar(0,'正在计算“幅度相关系数图”，请等待');
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
            tmp_A = s_imag_A(x_min:x_max,y_min:y_max);
            tmp_B = s_imag_B(x_min:x_max,y_min:y_max);
        else
            % 若上述四个条件都不满足，说明窗口不位于矩阵边界，则可以取到全部
            % （2N+1）*（2M+1）个点，因此直接用以下命令即可
            tmp_A = s_imag_A(pp-N:pp+N,qq-M:qq+M);
            tmp_B = s_imag_B(pp-N:pp+N,qq-M:qq+M);
        end
        
        % 下面计算该点（pp,qq）的幅度相关系数值：
        Coherence_imag(pp,qq) = sum(sum(abs(tmp_A).*abs(tmp_B)))/...
            ( sqrt(sum(sum(abs(tmp_A).^2)))*sqrt(sum(sum(abs(tmp_B).^2))) );
        clear tmp_A;clear tmp_B;
    end
    waitbar(pp/Naz);
end
close(h);

disp('完成“幅度相关系数图”计算');


end