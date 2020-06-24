function s_B_new = co_registration_V3_need(s_A,s_B,N_UP)
% 本函数用来完成“精配准”任务；
% 这是主程序 co_registration_V3（）函数中用来进行“精配准”的主要程序；
% 由于图像A和粗配准后的图像B较大，直接进行FFT等运算不容易实现，因此，我首先
% 进行了分块，将分块后的图像送入这个函数进行精配准，然后将配准结果返回；
% 在主程序co_registration_V3（）中，只需要调用该函数，最终再将几个分块结果合并即可。
% P.S.
%   由于进行了升采样，所以在配准取点的时候要记住进行“幅值恢复”！！！
%
% 输入变量：
% 1）s_A 是分块后，天线 A 的SLC-A；
% 2）s_B 是分块后，天线 B 的 SLC-B；
% 3）N_UP 是升采样的倍数；
%
% 输出变量：
% s_B_new 是与分块后的SLC-B相对应的，配准后的图像；
%
% 该程序截止至： 2015.02.02. 17:44 p.m.
 
%%
% ------------------------------------------------------------------------
%                   依次对 s_A 和 s_B 进行二维升采样
% ------------------------------------------------------------------------
[Naz,Nrg] = size(s_A);

% 对 s_A 进行二维升采样
S_IMAG_A = fft2(s_A);      % 对SLC-A进行二维FFT，下面进行二维频域补零。
S_IMAG_A_UP = [S_IMAG_A(1:round(Naz/2),:);
               zeros((N_UP-1)*Naz,Nrg);
               S_IMAG_A(round(Naz/2)+1:end,:)];
S_IMAG_A_UP2 = [S_IMAG_A_UP(:,1:round(Nrg/2)),zeros(N_UP*Naz,(N_UP-1)*Nrg),S_IMAG_A_UP(:,round(Nrg/2)+1:end)];
s_A_UP = ifft2(S_IMAG_A_UP2);   % ifft回到时域，完成对SLC-A的升采样。
clear S_IMAG_A;
clear S_IMAG_A_UP;
clear S_IMAG_A_UP2;

% 对 s_B 进行二维升采样
S_IMAG_B = fft2(s_B);      % 对SLC-B进行二维FFT，下面进行二维频域补零。
S_IMAG_B_UP = [S_IMAG_B(1:round(Naz/2),:);
               zeros((N_UP-1)*Naz,Nrg);
               S_IMAG_B(round(Naz/2)+1:end,:)];
S_IMAG_B_UP2 = [S_IMAG_B_UP(:,1:round(Nrg/2)),zeros(N_UP*Naz,(N_UP-1)*Nrg),S_IMAG_B_UP(:,round(Nrg/2)+1:end)];
s_B_UP = ifft2(S_IMAG_B_UP2);   % ifft回到时域，完成对SLC-B的升采样。
clear S_IMAG_B;
clear S_IMAG_B_UP;
clear S_IMAG_B_UP2;

%%
% ------------------------------------------------------------------------
%                               下面进行精配准
% ------------------------------------------------------------------------
s_B_new = zeros(Naz,Nrg);       % s_B_new 中存放该分块的配准结果。
% 为了程序的简洁，下面首先将矩阵最外层的数据直接存入 s_B_new 中。
% 注意，如果这些数据不太重要（比如属于不完全成像区），则无需进行进一步处理；
%      但如果这些数据也很重要，那么还要对这些最外层数据进行单独的配准操作；
%      在这里，我忽略了这部分数据的处理。
s_B_new(1,:) = s_B(1,:);
s_B_new(end,:) = s_B(end,:);
s_B_new(:,1) = s_B(:,1);
s_B_new(:,end) = s_B(:,end);

[num_azi,num_rag] = size(s_B_UP);           % 升采样后的矩阵大小

mm = 1;
for pp = (1+N_UP) : N_UP : num_azi-N_UP      % 第1行和最后一行不用计算，上面已经赋值了
    nn = 2;
    mm = mm + 1;
    for qq = (1+N_UP) : N_UP : num_rag-N_UP  % 第1列和最后一列不用计算，上面已经赋值了
        s_A_UP_pp_qq = s_A_UP(pp-N_UP:pp+N_UP,qq-N_UP:qq+N_UP);
        s_B_UP_pp_qq = s_B_UP(pp-N_UP:pp+N_UP,qq-N_UP:qq+N_UP);
        % =================================================================
        %{
        % 下面是严格考虑了沿方位向和斜距向两个方向的平移可能性，逐个平移量进行计算。
        % 计算量非常大
        R_pp_qq = zeros(2*N_UP+1,2*N_UP+1);
        
        for mov_azi = 1 : (2*N_UP+1)
            for mov_rag = 1 : (2*N_UP+1)
                tmp = circshift(s_B_UP_pp_qq,[mov_azi-21,mov_rag-21]);
                R_pp_qq(mov_azi,mov_rag) =  sum(sum(abs(s_A_UP_pp_qq).*abs(tmp)))/...
                    ( sqrt(sum(sum(abs(s_A_UP_pp_qq).^2)))*sqrt(sum(sum(abs(tmp).^2))) );
                clear tmp;
            end
        end
        
        [cor_R1,cor_p] = max(R_pp_qq);
        [cor_R2,cor_q] = max(max(R_pp_qq));
        % 最大值位置是 （ cor_p(cor_q) , cor_q ）
        % 对应的平移量应该是： （ cor_p(cor_q)-21 , cor_q-21 ）
        s_B_new(mm,nn) = s_B_UP(pp-(cor_p(cor_q)-21),qq-(cor_q-21));        
        %}
        %
        % 由于在我的仿真条件下，不存在方位向的移位，因此我只考虑沿斜距向的平移；
        % 简化后的程序如下：
        R_pp_qq = zeros(1,2*N_UP+1);
        for mov_rag = 1 : (2*N_UP+1)           
            tmp = circshift(s_B_UP_pp_qq,[0,mov_rag-(N_UP+1)]);   
            R_pp_qq(1,mov_rag) =  sum(sum(abs(s_A_UP_pp_qq).*abs(tmp)))/...
                ( sqrt(sum(sum(abs(s_A_UP_pp_qq).^2)))*sqrt(sum(sum(abs(tmp).^2))) );
            clear tmp;
        end
        
        [cor_R1,cor_p] = max(R_pp_qq);
        % 最大值位置是 （ cor_p ）
        % 对应的平移量应该是： （ 0 , cor_p-(N_UP+1) ）
        s_B_new(mm,nn) = N_UP*N_UP*s_B_UP(pp,qq-(cor_p-(N_UP+1)));
        % 这里乘以 N_UP*N_UP 的目的，是为了“幅值恢复”；
        % 因为进行了N_UP*N_UP的升采样，使得每个点的幅值降低到原来的1/(N_UP*N_UP)；
        % 这里在循环的同时进行了降采样，故一定还要同时进行幅值恢复，只需要通过
        % 乘以 N_UP*N_UP 即可；
        % =================================================================
        nn = nn+1;
    end
    mm              % 这是为了显示运行的进度
end

% s_B_new 即为所求；
% 将 s_B_new 作为函数返回值，得到该分块 s_B 的精配准结果。

end