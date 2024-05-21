%% Setup

%1: m_x
%2: m_y
%3: m_z

w_num=1.0;

k_num=0.40;

fs = 10001;


%% State Initializations

state_number=1;
vec0_filename = sprintf('data/init/%03d_mf.csv',state_number);
vec0 = readmatrix(vec0_filename,'Whitespace','()');



%% ODE
t_end = round(100/k_num);
tSpan = linspace(0, t_end, fs*t_end);
[T, Y] = ode15s(@(t,vec) order1(t,vec,k_num,w_num), tSpan, vec0);
Y=real(Y);


data = horzcat(k_num*T,Y);

filename = sprintf('data/evolution/state_%03d_k_%0.2f_N_inf_matlab_mf.csv',state_number,k_num);
writematrix(data, filename);


function diff = order1(t,vec,k_num,w_num)
    diff = [...
        k_num*vec(1)*vec(3); ...
        k_num*vec(2)*vec(3) - w_num*vec(3);...
        -k_num*vec(1)^2 - k_num*vec(2)^2 + w_num*vec(2);...
        ];
end
