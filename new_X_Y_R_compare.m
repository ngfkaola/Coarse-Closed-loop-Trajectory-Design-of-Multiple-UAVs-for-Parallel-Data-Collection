w_min = 100;
w_max = 4000;
w_gap = 50;
l_min = 100;
l_max = 4000;
l_gap = 50;
row = (w_max - w_min)/w_gap + 1;
column = (l_max - l_min)/l_gap + 1;
T_array = [];
X_array = [];%X数组存储w,l对应的X值用于compare2
Y_array = [];%Y数组存储w,l对应的Y值用于compare2
R_array = [];%R数组存储w,l对应的R值用于compare2
w_MonteCarlo = [];
l_MonteCarlo = [];
for w = w_min:w_gap:w_max
    for l = l_min:l_gap:l_max
        w_MonteCarlo(end+1) = w;
        l_MonteCarlo(end+1) = l;
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
[W,L] = meshgrid(w_min:w_gap:w_max,l_min:l_gap:l_max);
C = zeros(size(W));
surf(W,L,T_array,C);
colormap = ([0,0,1]);
xlabel('w');
ylabel('l');
zlabel('T');
% text(X_value,Y_value,R_value,' \leftarrow R最大值');
title('w-l-T');

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

function[value] = The_integral(B,distance_max,Alpha,D,noise,P,d0,K,R_number,X,Y,density)
prefix = 32*D*density*X*Y*log(2)/(Alpha*B);
p = 2/Alpha*log(P*K*d0^Alpha/noise);
func =  @(theta,x) 1./(p-log(x));
x_max = @(theta) distance_max.^2 - (2.*cos(theta).^2-1)./(2.*cos(theta).^2).*R_number.^2;
x_min = distance_max.^2 - R_number.^2;
value = prefix * integral2(func,0,pi/4,x_min,x_max);
end