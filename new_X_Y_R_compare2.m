w_min = 100;
w_max = 4000;
w_gap = 50;
l_min = 100;
l_max = 4000;
l_gap = 50;
row = (w_max - w_min)/w_gap + 1;
column = (l_max - l_min)/l_gap + 1;
T_average = [];
for i = 1:length(X_array(:))
    r = R_array(i);
    h = sqrt(distance_max^2 - r^2);
    X = X_array(i);
    Y = Y_array(i);
    %t_average = Time_average(w_MonteCarlo,l_MonteCarlo,i,B,Alpha,D,noise,P,d0,K,h,r,X,Y,density);
    t_average = Time_average(B,Alpha,D,noise,P,d0,K,h,r,X,Y,density);
    T_average(end+1) = t_average;
end
T_average = reshape(T_average,row,column);
[W,L] = meshgrid(w_min:w_gap:w_max,l_min:l_gap:l_max);
C = zeros(size(W));
surf(W,L,T_average,C);
colormap = ([1,0,0]);
%alpha(0.5);
xlabel('w');
ylabel('l');
zlabel('T');
% text(X_value,Y_value,R_value,' \leftarrow R最大值');
title('w-l-T');

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
%w = w_MonteCarlo(number);
%l = l_MonteCarlo(number);
T_average = mean(time(:))*4*X*Y;
end