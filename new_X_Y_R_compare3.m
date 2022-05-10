w_min = 100;
w_max = 4000;
w_gap = 100;
l_min = 100;
l_max = 4000;
l_gap = 100;
row = (w_max - w_min)/w_gap + 1;
column = (l_max - l_min)/l_gap + 1;
T_square = [];
T_square_compare_min = [];
T_square_compare_max = [];
T_square_average = [];
theta_all = [];
for w = w_min:w_gap:w_max
    l = w;
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
    theta = acos(R_value/distance_max)*180/pi;
    theta_all(end+1) = theta;
    %compare组的X,Y
    X1 = ceil(2^(-3/2)*w/sqrt(distance_max^2 - h_max^2));
    Y1 = ceil(2^(-3/2)*l/sqrt(distance_max^2 - h_max^2));
    X2 = ceil(2^(-3/2)*w/sqrt(distance_max^2 - h_min^2));
    Y2 = ceil(2^(-3/2)*w/sqrt(distance_max^2 - h_min^2));
    T_square(end+1) = T_value;
    %T_square_compare(end+1) = The_integral(B,distance_max,Alpha,D,noise,P,d0,K,R_max,X1,Y1,density);
    T_square_compare_min(end+1) = The_integral(B,distance_max,Alpha,D,noise,P,d0,K,R_min,X1,Y1,density);
    T_square_compare_max(end+1) = The_integral(B,distance_max,Alpha,D,noise,P,d0,K,R_max,X2,Y2,density);
    T_square_average(end+1) = Time_average(B,Alpha,D,noise,P,d0,K,sqrt(distance_max^2-R_value^2),R_value,X_value,Y_value,density);
end

make_config(T_square,T_square_compare_min,T_square_compare_max,T_square_average,w_min,w_gap,w_max);

function[T_average] = Time_average(B,Alpha,D,noise,P,d0,K,h,r,X,Y,density)
time = [];
loop_number = 100;
UE_real_number = ceil(2*r^2*density);
for i = 1:loop_number
    this_time = 0;
    for j = 1:UE_real_number
    %theta = rand(1,1)*2*pi;
    %R = unifrnd(0,r^2,1);
    %x = cos(theta)*sqrt(R);
    %y = sin(theta)*sqrt(R);
        x = unifrnd(sqrt(2)/(-2)*r,sqrt(2)/2*r);
        y = unifrnd(sqrt(2)/(-2)*r,sqrt(2)/2*r);
        distance = sqrt(x^2 + y^2+h^2);
        Pr = P*K*(d0/distance)^Alpha;
        SNR = Pr/noise;
        rate = B*log2(1+SNR);
        t = D/rate;
        this_time = this_time + t;
    end
    time(end+1) = this_time;
end
%T_average = mean(time(:))*(8*r^2*(X-1)*(Y-1)+(X-1)*2*sqrt(2)*r*(l-(Y-1)*2*sqrt(2)*r)+(Y-1)*2*sqrt(2)*r*(w-(X-1)*2*sqrt(2)*r)+(l-(Y-1)*2*sqrt(2)*r)*(w-(X-1)*2*sqrt(2)*r))*density;
T_average = mean(time(:))*4*X*Y;
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
        


function[] = make_config(T_square,T_square_compare_min,T_square_compare_max,T_square_average,w_min,w_gap,w_max)
x = w_min:w_gap:w_max;
plot(x,T_square,'-*b',x,T_square_compare_min,'-*r',x,T_square_compare_max,'-*k',x,T_square_average,'-*g');
% scatter(x,flip(t_max_h(1:6)),'*','r',x,flip(t_min_h(1:6)),'*','r',x,flip(t_max_h(7:12)),'*','r',x,flip(t_min_h(7:12)),'or',x,flip(t_max_h(13:18)),'xr',x,flip(t_min_h(13:18)),'xr');
legend('T.square','T.square.compare.Rmin','T.square.compare.Rmax','T.square.average');
xlabel('长度/m')
ylabel('时间/s')
end