t_compare_h = [];
t_average_h = [];
t_max_h = [];
t_min_h = [];
%theta_all = [];
for threshold = 5:5:15
    for h = 110:-10:50
        distance_max = (P*K/(threshold*noise))^(1/Alpha)*d0;
        r = sqrt(distance_max^2 - h^2);
        %theta = asin(h/distance_max)*180/pi;
        %theta_all(end+1) = theta;
        T_compare = Time_compare(B,Alpha,D,noise,P,d0,K,distance_max,r);
        t_compare_h(end+1) = T_compare;
        T = Time_average(B,Alpha,D,noise,P,d0,K,h,r);
        T_max = T(1);
        T_average = T(2);
        T_min = T(3);
        t_max_h(end+1) = T_max;
        t_average_h(end+1) = T_average;
        t_min_h(end+1) = T_min;
    end
end

make_config(t_compare_h,t_max_h,t_average_h,t_min_h)

function[T_compare] = Time_compare(B,Alpha,D,noise,P,d0,K,distance_max,r)
% h = 100;
p = 2/Alpha*log(P*K*d0^Alpha/noise);
% 积分函数
prefix = 2*D*log(2)/(Alpha*B)/r^2;
func =  @(x) 1./(p-log(x));
T_compare = prefix * integral(func,distance_max^2-r^2,distance_max^2);
%*int(1/(p-log(x)*(1/r^2)),distance_max^2-r^2,distance_max^2)
end

function[T] = Time_average(B,Alpha,D,noise,P,d0,K,h,r)
x_point = [];
y_point = [];
time = [];
UE_number = 10000;
for i = 1:UE_number
    theta = rand(1,1)*2*pi;
    R = unifrnd(0,r^2,1);
    x = cos(theta)*sqrt(R);
    y = sin(theta)*sqrt(R);
    x_point(end+1) = x;
    y_point(end+1) = y;
end
for i = 1:UE_number
    x0 = x_point(i);
    y0 = y_point(i);
    distance = sqrt(x0^2 + y0^2+h^2);
    Pr = P*K*(d0/distance)^Alpha;
    SNR = Pr/noise;
    rate = B*log2(1+SNR);
    t = D/rate;
    time(end+1) = t;
end
T_max = max(time);
T_average = mean(time(:));
T_min = min(time);
T = [T_max,T_average,T_min];
end

function[] = make_config(t_compare_h,t_max_h,t_average_h,t_min_h)
x = 50:10:110;
plot(x,flip(t_compare_h(1:7)),'-*b',x,flip(t_average_h(1:7)),'--*r',x,flip(t_compare_h(8:14)),'-ob',x,flip(t_average_h(8:14)),'--or',x,flip(t_compare_h(15:21)),'-xb',x,flip(t_average_h(15:21)),'--xr');
% scatter(x,flip(t_max_h(1:6)),'*','r',x,flip(t_min_h(1:6)),'*','r',x,flip(t_max_h(7:12)),'*','r',x,flip(t_min_h(7:12)),'or',x,flip(t_max_h(13:18)),'xr',x,flip(t_min_h(13:18)),'xr');
legend('\gamma=5,T_3','\gamma=5,T_\alpha','\gamma=10,T_3','\gamma=10,T_\alpha','\gamma=15,T_3','\gamma=15,T_\alpha');
xlabel('高度/m')
ylabel('时间/s')
end