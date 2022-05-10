%% 建立网格图
% 原始数据
x_axis = 0 : X_use;
y_axis = 0 : Y_use;
% 两根线的数据
x1 = [x_axis(1) x_axis(end)]';
y1 = [y_axis(1) y_axis(end)]';
% 所有线的xData
x2 = repmat(x1, 1, length(y_axis)-2);
x3 = repmat(x_axis(2) : x_axis(end-1), 2, 1);
xData = [x2 x3];
% 所有线的yData
y2 = repmat(y1, 1, length(x_axis)-2);
y3 = repmat(y_axis(2) : y_axis(end-1), 2, 1);
yData = [y3 y2];
% 绘图
h = line(xData, yData);
box on;
set(h, 'Color', 'k');

% 涂色方法（涂第i,j个格子）
% i = 5;
% j = 2;
% x = [0 1 1 0] + i - 1;
% y = [0 0 1 1] + j - 1;
% patch('xData', x, 'yData', y, 'FaceColor', 'r');

%% 开始选格子
%grid_number = sort(x(2:end));
grid_number = t(2:end);
index_matrix = zeros(2*length(grid_number),round(max(grid_number)));
grid_select = ones(X_use,Y_use);
color_select = ['r','g','b','y','m','c','k','w'];
for i = 1:X_use
    for j = 1:Y_use
        grid_select(i,j) = sqrt(i^2 + j^2);
    end
end
%grid_select
index = 1;
while (index <= length(grid_number))
    this_number = 1;
    max_value = max(max(grid_select));
    [i,j] = find(grid_select==max_value);
    %select_i,select_j为本次选定的格子
    select_i = i(1);
    select_j = j(1);
    %记录每次选的格子索引，奇数行为行索引，偶数行记录列索引，用于划线
    index_matrix(2*index-1,this_number) = select_i;
    index_matrix(2*index,this_number) = select_j;
    x_axis = [0 1 1 0] + select_i - 1;
    y_axis = [0 0 1 1] + select_j - 1;
    patch('xData', x_axis, 'yData', y_axis, 'FaceColor', color_select(index));
    grid_select(select_i,select_j) = 0;
    this_number = this_number + 1;
    %this_number
    grid_this_select = grid_select;
    for i =1:X_use
        for j = 1:Y_use
            if grid_select(i,j) ~= 0
                grid_this_select(i,j) = -1;
            end
        end
    end
    while this_number <= int8(grid_number(index))
        if select_i + 1 <= X_use && grid_this_select(select_i+1,select_j) == -1
            grid_this_select(select_i+1,select_j) = sqrt((select_i+1)^2 + select_j^2);
        end
        if select_i - 1 >= 1 && grid_this_select(select_i-1,select_j) == -1
            grid_this_select(select_i-1,select_j) = sqrt((select_i-1)^2 + select_j^2);
        end
        if select_j + 1 <= Y_use && grid_this_select(select_i,select_j+1) == -1
            grid_this_select(select_i,select_j+1) = sqrt((select_j+1)^2 + select_i^2);
        end
        if select_j - 1 >= 1 && grid_this_select(select_i,select_j-1) == -1
            grid_this_select(select_i,select_j-1) = sqrt((select_j-1)^2 + select_i^2);
        end
        max_value = max(max(grid_this_select));
        [i,j] = find(grid_this_select==max_value);
        select_i = i(1);
        select_j = j(1);
        index_matrix(2*index-1,this_number) = select_i;
        index_matrix(2*index,this_number) = select_j;
        x_axis = [0 1 1 0] + select_i - 1;
        y_axis = [0 0 1 1] + select_j - 1;
        patch('xData', x_axis, 'yData', y_axis, 'FaceColor', color_select(index));
        grid_this_select(select_i,select_j) = 0;
        grid_select(select_i,select_j) = 0;
        this_number = this_number + 1;
    end
    %this_number
    %grid_select
    %index
    index = index + 1;
    %index
end