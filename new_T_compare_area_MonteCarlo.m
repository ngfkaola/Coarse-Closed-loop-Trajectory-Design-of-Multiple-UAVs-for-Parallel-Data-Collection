w_min = 400;
w_max = 4000;
w_gap = 400;
l_min = 400;
l_max = 4000;
l_gap = 400;
row = (w_max - w_min)/w_gap + 1;
column = (l_max - l_min)/l_gap + 1;
T_array = [];
X_array = [];
Y_array = [];
R_array = [];
User_number_array = [];
for w = w_min:w_gap:w_max
    User_number_array(end+1) = w^2*density;
    for l = l_min:l_gap:l_max
        R_max = sqrt(distance_max^2 - h_min^2);
        T_value = inf;
        for X= ceil(2^(-3/2)*w/sqrt(distance_max^2 - h_min^2)):1:ceil(2^(-3/2)*w/sqrt(distance_max^2 - h_max^2))
            for Y = ceil(2^(-3/2)*l/sqrt(distance_max^2 - h_min^2)):1:ceil(2^(-3/2)*l/sqrt(distance_max^2 - h_max^2))
                R_min = max([sqrt(distance_max^2 - h_max^2),w/(2*sqrt(2)*X),l/(2*sqrt(2)*Y)]);
                R_number = select_R(B,distance_max,Alpha,D,noise,P,d0,K,X,Y,R_max,R_min,10^-5);
                value = The_integral(B,distance_max,Alpha,D,noise,P,d0,K,R_number,X,Y,density);
                if value < T_value
                    T_value = value;
                    X_value = X;
                    Y_value = Y;
                    R_value = R_number;
                end
            end
        end
    T_array(end+1) = T_value;
    X_array(end+1) = X_value;
    Y_array(end+1) = Y_value;
    R_array(end+1) = R_value;
    end
end
T_array = reshape(T_array,row,column);
X_array = reshape(X_array,row,column);
Y_array = reshape(Y_array,row,column);
R_array = reshape(R_array,row,column);
R_max = sqrt(distance_max^2 - h_min^2);
Un_max = 8;
Un_min = 4;
Un_gap = Un_min;
UAV_compare_t = zeros(2,length(X_array));
UAV_compare_t_lower = zeros(2,length(X_array));
UAV_compare_t_MonteCarlo = zeros(2,length(X_array));
for Un = Un_min:Un_gap:8 % UAV_number
    f = zeros(1,Un+1);
    f(1) = 1;
    intcon = 2:Un;
    Aeq = ones(1,Un+1);
    Aeq(1) = 0 ; 
    s = ones(1,Un);
    for i = 1:length(s)
%       if i <= 2
%           s(i) = 5;
%       elseif i <= 4
%           s(i) = 10;
%       else
%           s(i) = 15;
%       end
      s(i) = 10;
    end
    A = zeros(Un,Un+1);
    b = zeros(Un,1);
    %length(X_array)
    for select_row = 1:length(X_array)
        select_colmun = select_row;%暂且设定区域为正方形，X_use与Y_use均取矩阵对角线数据
        X_use = X_array(select_row,select_colmun);
        Y_use = Y_array(select_row,select_colmun);
        R_use = R_array(select_row,select_colmun);
        T_use = T_array(select_row,select_colmun);
        beq = X_use*Y_use;   
        for i= 1:size(A)
            A(i,1) = -1;
            A(i,i+1) = 4*sqrt(2)*R_use/s(i) + T_use/(X_use*Y_use);
        end
        %x = intlinprog(f,intcon,A,b,Aeq,beq);
        [t,fval] = intlinprog(f,intcon,A,b,Aeq,beq);
        UAV_compare_t(Un/Un_gap,select_row) = fval;
        UAV_compare_t_lower(Un/Un_gap,select_row) = calculate(Un,X_use,Y_use,R_use,s(i),T_use);
        grid_number = t(2:end);
        index_matrix = zeros(2*length(grid_number),round(max(grid_number)));
        grid_select = ones(X_use,Y_use);
        %color_select = ['r','g','b','y','m','c','k','w'];
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
            %patch('xData', x_axis, 'yData', y_axis, 'FaceColor', color_select(index));
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
                if select_i + 1 <= Y_use && grid_this_select(select_i+1,select_j) == -1
                    grid_this_select(select_i+1,select_j) = sqrt((select_i+1)^2 + select_j^2);
                end
                if select_i - 1 >= 1 && grid_this_select(select_i-1,select_j) == -1
                    grid_this_select(select_i-1,select_j) = sqrt((select_i-1)^2 + select_j^2);
                end
                if select_j + 1 <= X_use && grid_this_select(select_i,select_j+1) == -1
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
                %patch('xData', x_axis, 'yData', y_axis, 'FaceColor', color_select(index));
                grid_this_select(select_i,select_j) = 0;
                grid_select(select_i,select_j) = 0;
                this_number = this_number + 1;
            end
            index = index + 1;
        end
        %随着w,l的不断改变，后续w = w_min+w_gap*(select_row-1),l = l_min+l_gap*(select_row-1)
        User_number = User_number_array(select_row);
        grid_x_gap = (w_min+w_gap*(select_row-1))/X_use;
        grid_y_gap = (l_min+l_gap*(select_row-1))/Y_use;
        h = sqrt(distance_max^2 - R_use^2);
%         user_number_test = 0;
        Monte_number = 1; %撒点模拟1000次
        for Monte = 1:Monte_number
            %蒙特卡洛撒点
            User_x_position = randi([0,w_min+w_gap*(select_row-1)],1,User_number);
            User_y_position = randi([0,l_min+l_gap*(select_row-1)],1,User_number);
            UAV_time = zeros(1,Un);
            for select_UAV_index = 1:size(index_matrix,1)/2
                uav_time = 0;
                for select_grid_index = 1:size(index_matrix,2)
                    select_x = index_matrix(2*select_UAV_index-1,select_grid_index);
                    select_y = index_matrix(2*select_UAV_index,select_grid_index);
                    for user_index = 1:User_number
                        x0 = ceil(User_x_position(user_index)/grid_x_gap);
                        y0 = ceil(User_y_position(user_index)/grid_y_gap);
                        if x0 == 0
                            x0 = 1;
                        end
                        if y0 == 0
                            y0 = 1;
                        end
                        if x0 == select_x && y0 == select_y
                            %uav_time
                            uav_time = uav_time + MonteCarlo_calculate(x0,y0,P,K,Alpha,noise,B,D,d0,grid_x_gap,grid_y_gap,h);                   
%                             if select_UAV_index == 1
%                                 user_number_test = user_number_test + 1;
%                             end
                        end
                    end
                end
                UAV_time(select_UAV_index) = uav_time + grid_number(select_UAV_index)*4*(grid_x_gap/2)/s(select_UAV_index); %默认长宽间隔相等（正方形区域）
               
            end
            UAV_compare_t_MonteCarlo(Un/Un_gap,select_row) = UAV_compare_t_MonteCarlo(Un/Un_gap,select_row) + max(UAV_time);
%             UAV_compare_t_MonteCarlo
        end
        UAV_compare_t_MonteCarlo(Un/Un_gap,select_row) = UAV_compare_t_MonteCarlo(Un/Un_gap,select_row)/Monte_number;
    end
end
make_config(w_min,w_gap,w_max,UAV_compare_t,UAV_compare_t_lower,UAV_compare_t_MonteCarlo);

function[time] = calculate(Un,X_use,Y_use,r,s_i,T_use)
    z = X_use*Y_use/Un;
    time = 4*sqrt(2)*r*z/s_i + T_use*z/(X_use*Y_use);
end

function[time] = MonteCarlo_calculate(x0,y0,P,K,Alpha,noise,B,D,d0,grid_x_gap,grid_y_gap,h)
    x_accurate = mod(x0,grid_x_gap);
    y_accurate = mod(y0,grid_y_gap);
    if x_accurate <= grid_x_gap/2
        if y_accurate <= grid_y_gap/2
            distance = sqrt((x0-grid_x_gap/4)^2 + (y0-grid_y_gap/4)^2+h^2);
        else
            distance = sqrt((x0-grid_x_gap/4)^2 + (y0-grid_y_gap*3/4)^2+h^2);
        end
    else
        if y_accurate <= grid_y_gap/2
            distance = sqrt((x0-grid_x_gap*3/4)^2 + (y0-grid_y_gap/4)^2+h^2);
        else
            distance = sqrt((x0-grid_x_gap*3/4)^2 + (y0-grid_y_gap*3/4)^2+h^2);
        end
    end
    Pr = P*K*(d0/distance)^Alpha;
    SNR = Pr/noise;
    rate = B*log2(1+SNR);
    time = D/rate;  
    %rate,time
end

function[] = make_config(w_min,w_gap,w_max,UAV_compare_t,UAV_compare_t_lower,UAV_compare_t_MonteCarlo)
x = w_min:w_gap:w_max;
plot(x,UAV_compare_t(1,:),'-*b',x,UAV_compare_t_lower(1,:),'-*r',x,UAV_compare_t_MonteCarlo(1,:),'-*g',x,UAV_compare_t(2,:),'--*b',x,UAV_compare_t_lower(2,:),'--*r',x,UAV_compare_t_MonteCarlo(2,:),'--*g');
% scatter(x,flip(t_max_h(1:6)),'*','r',x,flip(t_min_h(1:6)),'*','r',x,flip(t_max_h(7:12)),'*','r',x,flip(t_min_h(7:12)),'or',x,flip(t_max_h(13:18)),'xr',x,flip(t_min_h(13:18)),'xr');
legend('UAV=4,T','UAV=4,T.lowerbound','UAV=4,T.MonentCarlo','UAV=8,T','UAV=8,T.lowerbound','UAV=8,T.MonentCarlo');
xlabel('区域边长')
ylabel('完成时间/s')
end

function[value] = The_integral(B,distance_max,Alpha,D,noise,P,d0,K,R_number,X,Y,density)
prefix = 32*D*density*X*Y*log(2)/(Alpha*B);
p = 2/Alpha*log(P*K*d0^Alpha/noise);
func =  @(theta,x) 1./(p-log(x));
x_max = @(theta) distance_max.^2 - (2.*cos(theta).^2-1)./(2.*cos(theta).^2).*R_number.^2;
x_min = distance_max.^2 - R_number.^2;
value = prefix * integral2(func,0,pi/4,x_min,x_max);
end

function[value] = select_R(B,distance_max,Alpha,D,noise,P,d0,K,X,Y,R_max,R_min,threshold)
prefix = 32*D*X*Y*log(2)/(Alpha*B);
p = 2/Alpha*log(P*K*d0^Alpha/noise);
%lambda = (2.*cos(theta).^2-1)./(2.*cos(theta).^2);
func_R_max = @(theta) -2.*(2.*cos(theta).^2-1)./(2.*cos(theta).^2).*R_max./(p-log(distance_max.^2-(2.*cos(theta).^2-1)./(2.*cos(theta).^2).*R_max.^2))+2.*R_max/(p-log(distance_max.^2-R_max.^2));
func_R_min = @(theta) -2.*(2.*cos(theta).^2-1)./(2.*cos(theta).^2).*R_min./(p-log(distance_max.^2-(2.*cos(theta).^2-1)./(2.*cos(theta).^2).*R_min.^2))+2.*R_min/(p-log(distance_max.^2-R_min.^2));
slope_R_max = prefix*integral(func_R_max,0,pi/4);
slope_R_min = prefix*integral(func_R_min,0,pi/4);
if slope_R_max == 0 && slope_R_min < 0
    value = R_min;
elseif slope_R_max == 0 && slope_R_min > 0
    value = R_max;
elseif slope_R_min == 0 && slope_R_max < 0
    value = R_max;
elseif slope_R_min == 0 && slope_R_max > 0
    value = R_min;
elseif slope_R_min > 0 && slope_R_max > 0
    value = R_min;
elseif slope_R_min < 0 && slope_R_max < 0
    value = R_max;
else
    R_number = (R_max+R_min)/2;
    func_R_number = @(theta) -2.*(2.*cos(theta).^2-1)./(2.*cos(theta).^2).*R_number./(p-log(distance_max.^2-(2.*cos(theta).^2-1)./(2.*cos(theta).^2)*R_number.^2))+2.*R_number/(p-log(distance_max.^2-R_number.^2));
    slope_R_number = prefix*integral(func_R_number,0,pi/4);
    if slope_R_min < 0 && slope_R_number < 0
        if slope_R_max - slope_R_number <= threshold
            value = R_number;
        else
            R_min = R_number;
            value = select_R(B,distance_max,Alpha,D,noise,P,d0,K,X,Y,R_max,R_min,threshold);
        end
    elseif slope_R_min < 0 && slope_R_number > 0
        if slope_R_number - slope_R_min <= threshold
            value = R_number;
        else
            R_max = R_number;
            value = select_R(B,distance_max,Alpha,D,noise,P,d0,K,X,Y,R_max,R_min,threshold);
        end
    end
end
end