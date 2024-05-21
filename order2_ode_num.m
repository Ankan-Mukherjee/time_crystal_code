%% Setup

%1: m_x
%2: m_y
%3: m_z
%4: ch_xx
%5: ch_xy
%6: ch_xz
%7: ch_yy
%8: ch_yz
%9: ch_zz

w_num=1.0;
k_num=0.40;
fs = 10001;


%% State Initializations
state_number=1;
vec0_filename = sprintf('data/init/%03d_so.csv',state_number);
vec0 = readmatrix(vec0_filename,'Whitespace','()');


%% ODE
t_end = round(100/k_num);
tSpan = linspace(0, t_end, fs*t_end);
%tSpan = [0, 1000]; 
[T, Y] = ode15s(@(t,vec) order2(t,vec,k_num,w_num), tSpan, vec0);
Y=real(Y);


data = horzcat(k_num*T,Y);

filename = sprintf('data/evolution/state_%03d_k_%0.2f_N_inf_matlab.csv',state_number,k_num);
writematrix(data, filename);

function diff = order2(t,vec,k_num,w_num)
    diff = [...
        k_num*vec(6) + k_num*vec(1)*vec(3); ...
        vec(8)*k_num + k_num*vec(2)*vec(3) - w_num*vec(3);...
        -k_num*vec(4) - k_num*vec(7) - k_num*vec(1)^2 - k_num*vec(2)^2 + w_num*vec(2);...
        2*vec(4)*k_num*vec(3) + 2*vec(6)*k_num*vec(1);...
        2*vec(5)*k_num*vec(3) + vec(6)*k_num*vec(2) - w_num*vec(6) + vec(8)*k_num*vec(1);...
        -2*vec(4)*k_num*vec(1) - 2*vec(5)*k_num*vec(2) + w_num*vec(5) + vec(6)*k_num*vec(3) + vec(9)*k_num*vec(1);...
        2*vec(7)*k_num*vec(3) + 2*vec(8)*k_num*vec(2) - 2*vec(8)*w_num;...
        -2*vec(5)*k_num*vec(1) - 2*vec(7)*k_num*vec(2) + vec(7)*w_num + vec(8)*k_num*vec(3) + vec(9)*k_num*vec(2) - vec(9)*w_num;...
        -4*vec(6)*k_num*vec(1) - 4*vec(8)*k_num*vec(2) + 2*vec(8)*w_num...
        ];
end
