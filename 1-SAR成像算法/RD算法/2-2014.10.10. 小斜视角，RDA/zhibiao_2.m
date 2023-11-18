function [PSLR,ISLR,IRW] = zhibiao_2(x,s_number,T)
    %
    % 针对函数 zhibiao(x,s_number,T) 进行改进，改进的内容主要为：
    % 在计算3dB点所对应的坐标时——用来计算 分辨率（IRW）
    % 函数 zhibiao() 中采用的是临近取整的办法，这不准确。
    % 下面利用的方法是将离 3dB 最近的两个点进行线性插值，来得到更准确的3dB点所对应的坐标。
    %
    % 输入变量：信号x，采样点数 s_number，T是信号的时域长度。
    % 该函数用来求解 x 的峰值旁瓣比(PSLR)，积分旁瓣比(ISLR),距离分辨率（IRW）。
    soo = x;                % 信号x
    N_buling = s_number;    % 采样点数
    
    soo_abs = abs(soo);     % soo的模
    [C,I] = max(soo_abs);   % 求输出 soo的模 中的最大值 C，位置 I；
    y = soo_abs.^2;         % 输出的平方， y = soo^2。
    
    x1 = 0;
    while (soo_abs(I-x1-1)-soo_abs(I-x1))<0
        M1 = x1;
        x1 = x1+1;
    end
    x2 = 0;
    while (soo_abs(I+x2+1)-soo_abs(I+x2))<0
        M2 = x2;
        x2 = x2+1;
    end
    
    P1 = I-1-M1;            % 主瓣和旁瓣分界点，左边的坐标是 P1。
    P2 = I+1+M2;            % 主瓣和旁瓣分界点，右边的左边是 P2。
    
    [D_left,Q_left] = max(soo_abs(1,1:P1));     % 最大旁瓣，值为 D_left，位置为 Q_left。（左边的那一个）。
    [D_right,Q_right] = max(soo_abs(1,P2:end)); % 最大旁瓣，值为 D_right，位置为 Q_right。（右边的那一个）。
    D = max(D_left,D_right);    % 比较左边和右边两者中的最大值，得到两侧旁瓣中最大的旁瓣，值为 D。
    
    PSLR = 20*log10(D/C);                       % 峰值旁瓣比
    ISLR = 10*log10((sum(y(1,P1/20:P1))+sum(y(1,P2:end)))/sum(y(1,P1:P2)));% 积分旁瓣比。
    
    %%%%%%%%%%%%%%%%%%%%%%%  以下是求 IRW  %%%%%%%%%%%%%%%%%%%%%%%%%
    M = ( 10^(-3/20) )*C;       % 3dB 带宽处的函数取值。
    % 下面是为了求找出与该函数值最接近的值的大小和坐标。
    z1 = abs(soo_abs(P1)-M);
    x1 = 1;
    z1_x1 = 0;
    for k1 = P1:I
        cha1 = abs(soo_abs(P1+x1)-M);
        if cha1<z1
            z1 = cha1;
            z1_x1 = x1; % z1_x1 是我们需要的，它的值是所求坐标与 P1 的偏移量。（左侧的）
        end
        x1 = x1+1;
    end
    
    z2 = abs(soo_abs(I)-M);
    x2 =1;
    z2_x2 = 0;
    for k2 = I:P2
        cha2 = abs(soo_abs(I+x2)-M);
        if cha2<z2
            z2 = cha2;
            z2_x2 = x2;% z2_x2 是我们需要的，它的值是所求坐标与 I 的偏移量。（右侧的）
        end
        x2 = x2+1;
    end
    
    Th_x1 = P1+z1_x1;% Th_x1 就是我们所求3dB带宽点左侧那个点的坐标（距离3dB最近的整数点）
    Th_x2 = I+z2_x2; % Th_x2 就是····3dB带宽点右侧······（距离3dB最近的整数点）
    % ------------------------------------------------------------------
    % 通过上面求得的距离3dB最近的整数点 Th_x1 和 Th_x2
    % 利用线性插值，得到3dB带宽所对应点数的更准确的值
    % 计算左边3dB的那个坐标
    if soo_abs(Th_x1)-M < 0
        x0_linear = Th_x1;
        x1_linear = Th_x1+1;
    else
        x0_linear = Th_x1-1;
        x1_linear = Th_x1;
    end
    Th_x1_real = ...
      (M-soo_abs(x1_linear))/(soo_abs(x0_linear)-soo_abs(x1_linear))*x0_linear...              
    + (M-soo_abs(x0_linear))/(soo_abs(x1_linear)-soo_abs(x0_linear))*x1_linear;
    
    % 计算右边3dB的那个坐标
    if soo_abs(Th_x2)-M > 0
        x0_linear = Th_x2;
        x1_linear = Th_x2+1;
    else
        x0_linear = Th_x2-1;
        x1_linear = Th_x2;
    end
    Th_x2_real = ...
      (M-soo_abs(x1_linear))/(soo_abs(x0_linear)-soo_abs(x1_linear))*x0_linear...              
    + (M-soo_abs(x0_linear))/(soo_abs(x1_linear)-soo_abs(x0_linear))*x1_linear;
    % ------------------------------------------------------------------
    width = Th_x2_real-Th_x1_real;  % width 就是通过此种方法求得的 3dB带宽。
    
    c = 3e8;% 光速 c=3e8 m/s。
    IRW = T/N_buling*width*c/2;% 注意在SAR中，分辨率是 C*T/2,其中T是脉冲宽度。
    % IRW_real为图像分辨率，原来的width的单位是采样间隔。一般图像分辨率单位取 m，要转换。
    % 注意到采样点数用的是 N_buling，因为频域补零后等效为升采样，采样率提高，采样点数应该跟正为 N_buling。
    