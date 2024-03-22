%% Set up

k = 0.4; % Set kappa as per master equation. Omega is set to 1
N = 60; % Number of spins
num_eigs = 80; % Number of eigenvalues to terminate to
operator = ["jz"]; % Operators producted (eg: use ["jx","jy","jz"] if you want Sx*Sy^2)
theta=pi/4;

%% Read Liovillians and compute eigs

operator_overlap_vec = zeros(1,num_eigs);
rho_overlap_vec = zeros(1,num_eigs);
    
L_filename = sprintf('data/liouvillian/lv_k_%0.2f_N_%03d.csv',k,N);
L = readmatrix(L_filename,'Whitespace','()');


op=1;
for opname = operator
    op_filename = sprintf('data/liouvillian/%s_N_%03d.csv',opname,N);
    oper = readmatrix(op_filename,'Whitespace','()');
    op = op*oper./(N/2);
end

op=op';
op=op(:);
state=zeros(N+1,1);
state(1)=1;
state = (expm(1j*theta*oper)*state+expm(1j*(pi-theta)*oper)*state)/sqrt(2);
rho=state*state';
rho=rho./trace(rho);
rho=rho';
rho=rho(:);


[V,D,W] = eig(L);
[d,sort_ind] = sort(diag(D),'descend','ComparisonMethod','real');
d = round(d,7);
r = V(:,sort_ind);
l = W(:,sort_ind);
temp = diag(D);
eigs = temp(sort_ind);
eigs = round(eigs(1:num_eigs),4);
r = r(:,1:num_eigs);
l = l(:,1:num_eigs);
operator_overlap_vec(1,:) = op'*r;
rho_overlap_vec(1,: ) = l'*rho;


%% Rounding
operator_overlap_vec = round(operator_overlap_vec,4);
rho_overlap_vec = round(rho_overlap_vec,4);
overlap_vec = round(operator_overlap_vec.*rho_overlap_vec,4);

%% Save overlaps
op_savefilename = sprintf('data/liouvillian/overlap_%s_numeigs_%02d_N_%03d.csv',strjoin(operator,""),num_eigs,N);
writematrix(operator_overlap_vec, op_savefilename);
rho_savefilename = sprintf('data/liouvillian/overlap_State_numeigs_%02d_N_%03d.csv',num_eigs,N);
writematrix(rho_overlap_vec, rho_savefilename);
overlap_savefilename = sprintf('data/liouvillian/overlap_%s_N00N_theta_%0.2f_numeigs_%02d_N_%03d.csv',strjoin(operator,""),theta,num_eigs,N);
writematrix(horzcat(abs(imag(eigs)),abs(overlap_vec')), overlap_savefilename);
