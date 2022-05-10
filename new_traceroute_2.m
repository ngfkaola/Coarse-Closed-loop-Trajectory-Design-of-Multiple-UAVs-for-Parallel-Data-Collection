x_point = [];
y_point = [];
for i = 1:X_use
    for j = 1:Y_use
        x_point(end+1) = i-0.25;
        x_point(end+1) = i-0.25;
        x_point(end+1) = i-0.75;
        x_point(end+1) = i-0.75;
        y_point(end+1) = j - 0.25;
        y_point(end+1) = j - 0.75;
        y_point(end+1) = j - 0.25;
        y_point(end+1) = j - 0.75;
    end
end
plot(x_point,y_point,'+');
axis([0 X_use 0 Y_use]);

hold on
%1:length(grid_number)
for species = 1:length(grid_number)
    available_number = 1;
    start_index = 1;
    max_index = size(index_matrix,2);
    for this_index = 2:max_index
        if index_matrix(2*species-1,this_index) ~= 0
            available_number = available_number + 1;
            %if sqrt((index_matrix(2*species-1,this_index)-0.75)^2 + (index_matrix(2*species,this_index)-0.75)^2) < sqrt((index_matrix(2*species-1,start_index)-0.75)^2 + (index_matrix(2*species,start_index)-0.75)^2)
            if index_matrix(2*species-1,this_index) + index_matrix(2*species,this_index) < index_matrix(2*species-1,start_index) + index_matrix(2*species,start_index)
                start_index = this_index;
            end
        end
    end
%     index_matrix(2*species-1,start_index) = 0;
%     index_matrix(2*species,start_index) = 0;
%     available_number = available_number - 1;
    start_x = index_matrix(2*species-1,start_index);
    start_y = index_matrix(2*species,start_index);
    plot(start_x-0.75,start_y-0.75,'ok','MarkerSize',10,'MarkerFaceColor',[0.5,0.5,0.5]);
    drawLine_Arrow([start_x-0.75 start_y-0.75],[start_x-0.25,start_y - 0.75],1,species);
    drawLine_Arrow([start_x-0.25 start_y-0.75],[start_x-0.25,start_y - 0.25],1,species);
    drawLine_Arrow([start_x-0.25 start_y-0.25],[start_x-0.75,start_y - 0.25],1,species);
    drawLine_Arrow([start_x-0.75 start_y-0.25],[start_x-0.75,start_y - 0.75],1,species);
    index_matrix(2*species-1,start_index) = 0;
    index_matrix(2*species,start_index) = 0;
    available_number = available_number - 1;
    [index_matrix,available_number] = draw(index_matrix,species,max_index,start_x,start_y,available_number);
end     
hold off

function drawLine_Arrow(start_point,end_point,select,species)
% 绘制带箭头的直线
% 从start_point到end_point画一箭头,arrColor箭头颜色,arrSize，箭头大小
%判断参数多少
switch select
    case 1
        fill_color_select = ['r','g','b','y','m','c','k','w'];
        arrColor  = fill_color_select(species);
        lineColor = 'b';
        arrowSize = 12;
        lineWidth = 1;
    case 2
        arrColor  = 'w';
        lineColor = 'w';
        arrowSize = 12;
        lineWidth = 1;
end
K= 0.5;          % 箭头比例系数
theta= pi / 8;    % 箭头角度
A1 = [cos(theta), -sin(theta);sin(theta), cos(theta)]; % 箭头左侧
theta = -theta;
A2 = [cos(theta), -sin(theta);sin(theta), cos(theta)]; % 箭头右侧
arrow = start_point' - end_point';
%使得箭头跟直线长短无关(固定值)
arrow_1(1)    = end_point(1,1)-arrowSize*cos(theta)+K*sin(theta);
arrow_2(1)    = end_point(1,1)-arrowSize*cos(theta)-K*sin(theta);
arrow_1(2)    = end_point(1,2)-arrowSize*sin(theta)-K*cos(theta);
arrow_2(2)    = end_point(1,2)-arrowSize*sin(theta)+K*cos(theta);
arrow_1= A1 * arrow;
arrow_2= A2 * arrow;
arrow_1= K * arrow_1 + end_point'; % 箭头的边的x坐标
arrow_2= K * arrow_2 + end_point'; % 箭头的边的y坐标
hold on;
grid on;
axis equal;
plot([start_point(1), end_point(1)], [start_point(2), end_point(2)],lineColor,'lineWidth',lineWidth);
% 三角箭头(填充)
triangle_x= [end_point(1),arrow_1(1),arrow_2(1),end_point(1)];
triangle_y= [end_point(2),arrow_1(2),arrow_2(2),end_point(2)];
this_fill = fill(triangle_x,triangle_y,arrColor);
if arrColor == 'w'
    set(this_fill,'edgecolor',arrColor); 
end
% 线段箭头(不填充)
% plot([arrow_1(1), end_point(1)], [arrow_1(2), end_point(2)],color,'lineWidth',arrowSize);
% plot([arrow_2(1), end_point(1)], [arrow_2(2), end_point(2)], color,'lineWidth',arrowSize);
hold on;
end 

%%递归函数画路线
function[index_matrix,available_number] = draw(index_matrix,species,max_index,last_x,last_y,available_number)
    if available_number == 0
        return 
    end
    for this_index = 1:max_index
        x = index_matrix(2*species-1,this_index);
        y = index_matrix(2*species,this_index);
        if x == last_x&&y == last_y+1 || x == last_x && y == last_y-1 || x == last_x-1&&y==last_y || x == last_x+1&&y == last_y
            if x == last_x + 1 && y == last_y 
                drawLine_Arrow([x-0.75 y-0.75],[x-0.25,y - 0.75],1,species);
                drawLine_Arrow([x-0.25 y-0.75],[x-0.25,y - 0.25],1,species);
                drawLine_Arrow([x-0.25 y-0.25],[x-0.75,y - 0.25],1,species);
                %drawLine_Arrow([x-0.75 y-0.25],[x-0.75,y - 0.75],1,species);
                drawLine_Arrow([last_x-0.25 last_y-0.75],[x-0.75,y - 0.75],1,species);
                drawLine_Arrow([x-0.75 y-0.25],[last_x-0.25,last_y - 0.25],1,species);
                drawLine_Arrow([last_x-0.25 last_y-0.75],[last_x-0.25,last_y - 0.25],2,species);
            elseif x == last_x && y == last_y + 1
                %drawLine_Arrow([x-0.75 y-0.75],[x-0.25,y- 0.75],1,species);
                drawLine_Arrow([x-0.25 y-0.75],[x-0.25,y - 0.25],1,species);
                drawLine_Arrow([x-0.25 y-0.25],[x-0.75,y - 0.25],1,species);
                drawLine_Arrow([x-0.75 y-0.25],[x-0.75,y - 0.75],1,species);
                drawLine_Arrow([last_x-0.25 last_y-0.25],[x-0.25,y - 0.75],1,species);
                drawLine_Arrow([x-0.75 y-0.75],[last_x-0.75,last_y - 0.25],1,species);
                drawLine_Arrow([last_x-0.25 last_y-0.25],[last_x-0.75,last_y - 0.25],2,species);
            elseif x == last_x-1 && y == last_y 
                drawLine_Arrow([x-0.75 y-0.75],[x-0.25,y- 0.75],1,species);
                %drawLine_Arrow([x-0.25 y-0.75],[x-0.25,y - 0.25],1,species);
                drawLine_Arrow([x-0.25 y-0.25],[x-0.75,y - 0.25],1,species);
                drawLine_Arrow([x-0.75 y-0.25],[x-0.75,y - 0.75],1,species);
                drawLine_Arrow([last_x-0.75 last_y-0.25],[x-0.25,y - 0.25],1,species);
                drawLine_Arrow([x-0.25 y-0.75],[last_x-0.75,last_y - 0.75],1,species);
                drawLine_Arrow([last_x-0.75 last_y-0.25],[last_x-0.75,last_y - 0.75],2,species);
            elseif x == last_x && y == last_y - 1
                drawLine_Arrow([x-0.75 y-0.75],[x-0.25,y- 0.75],1,species);
                drawLine_Arrow([x-0.25 y-0.75],[x-0.25,y - 0.25],1,species);
                %drawLine_Arrow([x-0.25 y-0.25],[x-0.75,y - 0.25],1,species);
                drawLine_Arrow([x-0.75 y-0.25],[x-0.75,y - 0.75],1,species);
                drawLine_Arrow([last_x-0.75 last_y-0.75],[x-0.75,y - 0.25],1,species);
                drawLine_Arrow([x-0.25 y-0.25],[last_x-0.25,last_y - 0.75],1,species);
                drawLine_Arrow([last_x-0.75 last_y-0.75],[last_x-0.25,last_y - 0.75],2,species);          
            end
            index_matrix(2*species-1,this_index) = 0;
            index_matrix(2*species,this_index) = 0;
            available_number = available_number - 1;
            last_x = x;
            last_y = y;
            [index_matrix,available_number] = draw(index_matrix,species,max_index,last_x,last_y,available_number);
        end
    end
end