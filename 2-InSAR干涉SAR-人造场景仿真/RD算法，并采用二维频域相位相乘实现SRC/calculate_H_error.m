function [H_area_ref,H_error] = calculate_H_error(X_area,Y_area,H_area,Parameter_calculate_H_error,TYPE)
% 本函数用来计算与“干涉处理得到的地距平面高程图 H_area”相对应的理论高程图
% H_area_ref，以及干涉处理结果和理论高程图的误差 H_error；
%
% 方法：
%   类似在生成原始数据时，计算场景模型高程信息的方法，直接利用
%   （ X_area , Y_area , H_area）这三个信息计算得到相对应的理论高程图及误差；
%   缺点在于我依然没有将斜距平面的高程信息 H_area 按照地距坐标系进行投影，也就
%   是说地距平面依然不是等间隔的。
%
% 输入数据：
% 1）X_area  是干涉处理最终得到 H_area 后，计算得到的与之相对应的地面 x 轴坐标；
% 2）Y_area  同上，是与之相对应的地面 y 轴坐标；
% 3）H_area  计算得到的高程信息（是我干涉处理的结果）；
% 4）Parameter_calculate_H_error 是需要的一些参数：
%    a）（1,1）表示：x1      场景中心的地距 X 轴坐标；
%    b）（2,1）表示：y1      场景中心的地距 Y 轴坐标；
%    c）（3,1）表示：r_cone  这是所设置的圆锥场景的半径；
%              如果是平地场景，由于不涉及到该值，可令输入值为0，代入函数计算；
%    d）（4,1）表示：r_cone_Height 这是所设置的圆锥场景中心的高度（最大高度值）；
%              如果是平地场景，由于不涉及到该值，可令输入值为0，代入函数计算；
% 5）TYPE    有两种取值：
%    a）TYPE == 1，代表平地场景，此时场景高程信息为0，是最简单的情况；
%    b）TYPE == 2，代表圆锥场景，此时需要一定的方法计算得到所要的信息；
%
% 输出数据：
% 1）H_area_ref      是与 H_area 相对应的理论高程图 H_area_ref；
% 2）H_error         干涉处理结果和理论高程图的误差；
%
% 程序版本截止到：
% 2015.01.29. 11:28 a.m.

%%
% --------------------------------------------------------------------
%                    计算与 H_area 相对应的理论高程图
%                   计算干涉处理结果与理论高程图的误差
% --------------------------------------------------------------------
[num_y,num_x] = size(H_area);
x1 = Parameter_calculate_H_error(1,1);
y1 = Parameter_calculate_H_error(2,1);
r_cone = Parameter_calculate_H_error(3,1);
r_cone_Height = Parameter_calculate_H_error(4,1);


% ********************************************************************
%                               平地场景时
% ********************************************************************
if TYPE == 1    % 平地场景
    % 此时 H_area_ref == 0 ，因此不需要再进行计算，可以直接得到处理误差
    H_area_ref = 0.*zeros(num_y,num_x);
    H_error = H_area - H_area_ref;
end

% ********************************************************************
%                               圆锥场景时
% ********************************************************************
if TYPE == 2    % 圆锥场景
    % 此时 H_area_ref 需要逐点计算，以便与 H_area 的坐标相对应 
    % 首先计算每点的地面斜距（ 相对于场景中心(x1,y1) ）
    R_target_all = sqrt( ((X_area - x1).^2) + ((Y_area - y1).^2) );
    % 下面利用“逻辑1寻访”的功能，来生成与 H_area 相对应的理论高程图
    L = (R_target_all <= r_cone);       % 赋值号右边：进行关系比较，产生逻辑结果；
                                        % 产生与 R_target_all 维数大小相同的
                                        % “0，1”逻辑数组；1表示“真”
                                        % 在此 L 数组中取 1 的位置对应的 R_target_all
                                        % 数组元素小于等于 r_cone ；

    H_area_ref = zeros(num_y,num_x);    % 用来存放理论高程信息，矩阵；
                                        % 初始值都为0；
    % 下面计算与 H_area 相对应的理论高程图 H_area_ref：
    for kk = 1:num_y
        for ll = 1:num_x
            if L(kk,ll) == 0
                continue;
            else
                H_area_ref(kk,ll) = (r_cone - R_target_all(kk,ll))*(r_cone_Height/r_cone);
            end      
        end
    end
    % 计算结束，这样我们就得到了 H_area_ref；
    % 下面利用 H_area 和 H_area_ref 计算干涉处理结果与理论高程图的误差
    H_error = H_area - H_area_ref;
    
    % 至此，程序结束。
end

end