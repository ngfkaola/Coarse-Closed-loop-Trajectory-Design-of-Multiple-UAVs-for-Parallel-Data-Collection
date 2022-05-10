w_min = 100;
w_max = 4000;
w_gap = 50;
l_min = 100;
l_max = 4000;
l_gap = 50;
row = (w_max - w_min)/w_gap + 1;
colmun = (l_max - l_min)/l_gap + 1;
T = [];
for w = w_min:w_gap:w_max
    for l = l_min:l_gap:l_max
        R_max = sqrt(distance_max^2 - h_min^2);
        X = ceil(2^(-3/2)*w/sqrt(distance_max^2 - h_min^2));
        Y = ceil(2^(-3/2)*l/sqrt(distance_max^2 - h_min^2));
        % R_value = max([sqrt(distance_max^2 - h_max^2),w/(2*sqrt(2)*X),l/(2*sqrt(2)*Y)]);
        prefix = 32*D*density*X*Y*log(2)/(Alpha*B);
        p = 2/Alpha*log(P*K*d0^Alpha/noise);
        func =  @(theta,x) 1./(p-log(x));
        x_max = @(theta) distance_max.^2 - (2.*cos(theta).^2-1)./(2.*cos(theta).^2).*R_max.^2;
        x_min = distance_max.^2 - R_max.^2;
        T_value = prefix * integral2(func,0,pi/4,x_min,x_max);  
        T(end+1) = T_value;
    end
end
T = reshape(T,row,colmun);
[W,L] = meshgrid(w_min:w_gap:w_max,l_min:l_gap:l_max);
C = zeros(size(W));
surf(W,L,T,C);
colormap = ([1,0,0]);
alpha(0.5);
xlabel('w');
ylabel('l');
zlabel('T');
% text(X_value,Y_value,R_value,' \leftarrow R最大值');
title('w-l-T');