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
Un = 6; % UAV_number
r = R_value;
T_use = T_value; %T_use为P5问题中的T*
X_use = X_value;%X_use为P5问题中的X*
Y_use = Y_value;%Y_use为P5问题中的Y*
f = zeros(1,Un+1);
f(1) = 1;
intcon = 2:Un;
Aeq = ones(1,Un+1);
Aeq(1) = 0 ; 
s = ones(1,Un);
for i = 1:length(s)
    if i <= 2
        s(i) = 8;
    elseif i <= 4
        s(i) = 10;
    else
        s(i) = 12;
    end
end
A = zeros(Un,Un+1);
b = zeros(Un,1);
beq = X_use*Y_use;   
for i= 1:size(A)
    A(i,1) = -1;
    A(i,i+1) = 4*sqrt(2)*r/s(i) + T_use/(X_use*Y_use);
end
%x = intlinprog(f,intcon,A,b,Aeq,beq);
[t,fval] = intlinprog(f,intcon,A,b,Aeq,beq);

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