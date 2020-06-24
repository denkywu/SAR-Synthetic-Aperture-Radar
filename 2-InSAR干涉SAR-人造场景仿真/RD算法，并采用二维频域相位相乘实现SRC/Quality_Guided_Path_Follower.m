function PHY_after_unwrapping = Quality_Guided_Path_Follower(PHY_s_after_X_filtering,Coherence_imag)
% 本函数采用“质量指导的路径跟踪法”进行相位解缠绕
% 参考：
% 1）王超等的《星载合成孔径雷达干涉测量》第111页，第4.4.2节“质量指导的路径跟踪法”；
% 2）《Two Dimensional Phase Unwrapping Theory Algorithms and Software》；
%
% 输入数据：
%   1）PHY_s_after_X_filtering  表示经过某种滤波方法后的干涉相位图；
%   2）Cohorence_imag 质量图（我采用的是幅度相关系数图）
% 输出数据：
%   PHY_after_unwrapping  表示相位解缠绕后的相位图；
%
% 程序截止到： 2015.01.07. 20:33 p.m.

%%
disp('正在使用“质量指导的路径跟踪法”进行相位解缠绕，请等待');

[Naz,Nrg] = size(PHY_s_after_X_filtering);      % 数据大小
PHY_after_unwrapping = zeros(Naz,Nrg);          % 存放解缠绕相位结果  
% ------------------------------------------------------------------------
%                           选择起始点
%               并从该起始点开始相位解缠绕的第一步
% ------------------------------------------------------------------------
% **********************************************************************
%                              选择起始点
% **********************************************************************
Coherence_imag2 = zeros(Naz,Nrg);
Coherence_imag2(2:end-1,2:end-1) = Coherence_imag(2:end-1,2:end-1);

[corR1,cor_p] = max(Coherence_imag2);
[corR2,cor_q] = max(max(Coherence_imag2));
% 则最大值位置是：（cor_p(cor_q)，cor_q）
clear Coherence_imag2;

% **********************************************************************
%                         以此起始点开始相位解缠绕
% **********************************************************************
% （1）首先存下其坐标
unwrap_start_x = cor_p(cor_q);
unwrap_start_y = cor_q;
% 这是相位解缠起点坐标，也是质量图最大值位置

% （2）直接将待解缠绕相位的，质量图最大值位置所对应的值存放在解缠绕后的结果中
%      作为起点
PHY_after_unwrapping(unwrap_start_x,unwrap_start_y) = ...
    PHY_s_after_X_filtering(unwrap_start_x,unwrap_start_y);

% （3）定义一个新的矩阵，用来标记状态（待解缠绕，已解缠绕，邻接列等）
Unwrap_Mode = zeros(Naz,Nrg);   % 初始值为0，用来表示待解缠绕
Unwrap_Mode(unwrap_start_x,unwrap_start_y) = 1;
% 将起点对应的位置赋值为1，代表已解缠绕

% **********************************************************************
%                           下面更新“邻接列”
% **********************************************************************
% 对于第一个点而言，更新邻接列很简单，只需将周围四个点都标记为 NaN 即可
Unwrap_Mode(unwrap_start_x-1,unwrap_start_y) = NaN;
Unwrap_Mode(unwrap_start_x,unwrap_start_y-1) = NaN;
Unwrap_Mode(unwrap_start_x,unwrap_start_y+1) = NaN;
Unwrap_Mode(unwrap_start_x+1,unwrap_start_y) = NaN;

%%
% ------------------------------------------------------------------------
%             下面进行 while 循环，以对整个矩阵执行相位解缠任务
% ------------------------------------------------------------------------
NUM_while = 1;
while( any(any(isnan(Unwrap_Mode))) )
    % （1）这时要做的是在邻接列中（标记NaN的），寻找质量图最大的位置
    %      因此，首先定位NaN——采用寻访NaN的办法
    LR = isnan(Unwrap_Mode);
    % 此时，在LR中，所有NaN的位置被标记为1，其余的标记为0；
    % ===========================================================
    % 下面这部分程序代码，是用来求找在标记了NaN的数据中，即邻接列中，
    % 最大值的位置，以便进行后续解缠绕。
    % 但是，下面这部分代码直接用for循环，效率太低，必须改进。
    %{
    [rj,cj] = find(LR);
    
    m = length(rj);
    tmp = zeros(1,m);
    for pp = 1:m
        tmp(1,pp) = Coherence_imag(rj(pp),cj(pp));
    end
    [max_tmp,max_tmp_num] = max(tmp);
    %}
    % 我将上面的结果进行改进，不使用for循环，如下：
    [rj,cj] = find(LR);
    [max_tmp,max_tmp_num] = max(Coherence_imag(find(LR)));
    % ===========================================================
    
    % 则标记 NaN 的邻接列中，相关系数最大的，其坐标为：
    %   （rj(max_tmp_num),cj(max_tmp_num)）
    % 下面通过一个判断语句，使得while（）循环体中跳过矩阵最外层
    if rj(max_tmp_num) == 1 || rj(max_tmp_num) == Naz || ...
        cj(max_tmp_num) == 1 || cj(max_tmp_num) == Nrg
        % 将状态从 NaN 改为 0，排除出邻接列
        Unwrap_Mode(rj(max_tmp_num),cj(max_tmp_num)) = 0;
        % 跳过此次循环
        continue;
    end
    
    % （2）邻接列中，对应质量图最大的元素，其待解缠绕相位是：
    wait_unwrap = PHY_s_after_X_filtering(rj(max_tmp_num),cj(max_tmp_num));
    % 位置在：（ rj(max_tmp_num),cj(max_tmp_num) ）
    
    % （3）按照一维解缠绕方法直接解缠
    %      这是要首先判断由哪一点进行解缠
    %      我目前的方法是：依次判断该点周围四点，只要有一个点的状态（Unwrap_Mode）
    %      其值是1，就用来解缠即可     ——是否有更好的方法？？
    if Unwrap_Mode(rj(max_tmp_num)-1,cj(max_tmp_num)) == 1
        delta_unwrap = wait_unwrap - PHY_s_after_X_filtering(rj(max_tmp_num)-1,cj(max_tmp_num));
        if delta_unwrap > pi
            delta_unwrap = delta_unwrap - 2*pi;
        end
        if delta_unwrap < -1*pi
            delta_unwrap = delta_unwrap + 2*pi;
        end
        PHY_after_unwrapping(rj(max_tmp_num),cj(max_tmp_num)) = ...
            PHY_after_unwrapping(rj(max_tmp_num)-1,cj(max_tmp_num))...
            + delta_unwrap;
    else
        if Unwrap_Mode(rj(max_tmp_num),cj(max_tmp_num)-1) == 1
            delta_unwrap = wait_unwrap - PHY_s_after_X_filtering(rj(max_tmp_num),cj(max_tmp_num)-1);
            if delta_unwrap > pi
                delta_unwrap = delta_unwrap - 2*pi;
            end
            if delta_unwrap < -1*pi
                delta_unwrap = delta_unwrap + 2*pi;
            end            
            PHY_after_unwrapping(rj(max_tmp_num),cj(max_tmp_num)) = ...
                PHY_after_unwrapping(rj(max_tmp_num),cj(max_tmp_num)-1)...
                + delta_unwrap;
        else
            if Unwrap_Mode(rj(max_tmp_num),cj(max_tmp_num)+1) == 1
                delta_unwrap = wait_unwrap - PHY_s_after_X_filtering(rj(max_tmp_num),cj(max_tmp_num)+1);
                if delta_unwrap > pi
                    delta_unwrap = delta_unwrap - 2*pi;
                end
                if delta_unwrap < -1*pi
                    delta_unwrap = delta_unwrap + 2*pi;
                end  
                PHY_after_unwrapping(rj(max_tmp_num),cj(max_tmp_num)) = ...
                    PHY_after_unwrapping(rj(max_tmp_num),cj(max_tmp_num)+1)...
                    + delta_unwrap;
            else
                if Unwrap_Mode(rj(max_tmp_num)+1,cj(max_tmp_num)) == 1
                    delta_unwrap = wait_unwrap - PHY_s_after_X_filtering(rj(max_tmp_num)+1,cj(max_tmp_num));
                    if delta_unwrap > pi
                        delta_unwrap = delta_unwrap - 2*pi;
                    end
                    if delta_unwrap < -1*pi
                        delta_unwrap = delta_unwrap + 2*pi;
                    end    
                    PHY_after_unwrapping(rj(max_tmp_num),cj(max_tmp_num)) = ...
                        PHY_after_unwrapping(rj(max_tmp_num)+1,cj(max_tmp_num))...
                        + delta_unwrap;
                else
                    disp('解缠绕错误，该点周围没有已解缠的点');
                    return;
                end
            end
        end
    end
    
    % （4）更新邻接列
    %   1）首先令上面刚解缠的结果，其状态更新为1；
    Unwrap_Mode(rj(max_tmp_num),cj(max_tmp_num)) = 1;
    %   2）将周围的点更新进邻接列，同时要加上条件判断
    %      因为在循环中，我不知道这些点是如何“生长的”，因此，要对该点周围四个
    %      点都进行状态判断。若是0，则更新为 NaN；
    %                       若是 1 或者 NaN，保持不变即可
    if Unwrap_Mode(rj(max_tmp_num)-1,cj(max_tmp_num)) == 0
        Unwrap_Mode(rj(max_tmp_num)-1,cj(max_tmp_num)) = NaN;
    end
    if Unwrap_Mode(rj(max_tmp_num),cj(max_tmp_num)-1) == 0
        Unwrap_Mode(rj(max_tmp_num),cj(max_tmp_num)-1) = NaN;
    end
    if Unwrap_Mode(rj(max_tmp_num),cj(max_tmp_num)+1) == 0
        Unwrap_Mode(rj(max_tmp_num),cj(max_tmp_num)+1) = NaN;
    end
    if Unwrap_Mode(rj(max_tmp_num)+1,cj(max_tmp_num)) == 0
        Unwrap_Mode(rj(max_tmp_num)+1,cj(max_tmp_num)) = NaN;
    end
    
    % （5）进行以下清除操作，避免在下一次循环中引起错误
    %      这些值在每一次循环中都只是作为中间变量存在。
    clear LR;
    clear rj;clear cj;
    clear tmp;clear max_tmp;clear max_tmp_num;
    clear wait_unwrap;
    clear delta_unwrap;
    
    NUM_while
    NUM_while = NUM_while + 1;      % 这个用来显示循环次数。便于了解程序进行情况。
    
    % （6）到此，while循环体结束，下面回到while（）开头，继续进行；
    %     直到除了最外层外，都完成解缠绕

end

%%
% ------------------------------------------------------------------------
%                   最后，对最外层执行相位解缠绕
% ------------------------------------------------------------------------
% （1）由第二行解缠绕第一行
for jj = 2:(Nrg-1)
    delta_jj = PHY_s_after_X_filtering(1,jj) - PHY_s_after_X_filtering(2,jj);
    if delta_jj > pi
        delta_jj = delta_jj - 2*pi;
    end
    if delta_jj < -1*pi
        delta_jj = delta_jj + 2*pi;
    end
    PHY_after_unwrapping(1,jj) = PHY_after_unwrapping(2,jj) + delta_jj;
end

% （2）由倒数第二行解缠最后一行
for jj = 2:(Nrg-1)
    delta_jj = PHY_s_after_X_filtering(Naz,jj) - PHY_s_after_X_filtering(Naz-1,jj);
    if delta_jj > pi
        delta_jj = delta_jj - 2*pi;
    end
    if delta_jj < -1*pi
        delta_jj = delta_jj + 2*pi;
    end
    PHY_after_unwrapping(Naz,jj) = PHY_after_unwrapping(Naz-1,jj) + delta_jj;
end

% （3）由第二列解缠第一列
for jj = 1:Naz
    delta_jj = PHY_s_after_X_filtering(jj,1) - PHY_s_after_X_filtering(jj,2);
    if delta_jj > pi
        delta_jj = delta_jj - 2*pi;
    end
    if delta_jj < -1*pi
        delta_jj = delta_jj + 2*pi;
    end
    PHY_after_unwrapping(jj,1) = PHY_after_unwrapping(jj,2) + delta_jj;
end

% （4）由倒数第二列解缠最后一列
for jj = 1:Naz
    delta_jj = PHY_s_after_X_filtering(jj,Nrg) - PHY_s_after_X_filtering(jj,Nrg-1);
    if delta_jj > pi
        delta_jj = delta_jj - 2*pi;
    end
    if delta_jj < -1*pi
        delta_jj = delta_jj + 2*pi;
    end
    PHY_after_unwrapping(jj,Nrg) = PHY_after_unwrapping(jj,Nrg-1) + delta_jj;
end

disp('已完成相位解缠绕');

end
