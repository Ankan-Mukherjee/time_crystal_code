syms m_x(t) m_y(t) m_z(t) ch_xx(t) ch_xy(t) ch_xz(t) ch_yy(t) ch_yz(t) ch_zz(t) k w;

%% Setup

%1: m_x
%2: m_y
%3: m_z

w_num=1.0;
k_list = [0.5:0.01:0.99, 1.0:0.05:5.0];
theta = 0.79;
phi = 0;
fs = 1001;


%% State Initializations
m_x_0 = 0.7071;
m_y_0 = 0.0;
m_z_0 = 0.0;


% vec0 = vec_0_from_state_prep;
vec0 = [m_x_0,m_y_0,m_z_0];


%% ODE
for k_num = k_list
    "k = "+num2str(k_num)
    t_end = 50;
    tSpan = linspace(0, t_end, fs*t_end); 
    [T, Y] = ode15s(@(t,vec) order1(t,vec,k_num,w_num), tSpan, vec0);
    Y=real(Y);
    data = horzcat(k_num*T,Y);
    filename = sprintf('data/N00N_theta_%0.2f_k_%0.2f_N_inf_matlab_mf.csv',theta,k_num);
    writematrix(data, filename);
end

function diff = order1(t,vec,k_num,w_num)
    diff = [...
        k_num*vec(1)*vec(3); ...
        k_num*vec(2)*vec(3) - w_num*vec(3);...
        -k_num*vec(1)^2 - k_num*vec(2)^2 + w_num*vec(2);...
        ];
end
