function [] = Problem_1b()

Set_Default_Plot_Properties();

% Solution domain
Nx = 101;
x0 = 0;
xf = 1;
x = linspace(x0, xf, Nx)';
dx = x(2) - x(1);

% Solution boundary conditions
u0 = 0;
uf = 0;

% Karhunen-Loeve expansion (KLE) options
sigma = 2.0;    % Standard deviation
ell = 2.0;      % Correlation length
a = 1/2;        % Support of eigenproblem
d = 2;          % Number of terms

% K_i PCE
pk = 14; % Total order
index_pc = nD_polynomial_array(d, pk);

% Calculate the PCE expansion for K_i(x)
[Ki, Pk] = Compute_Ki(pk, sigma, ell, a, d, x);

%%%
% Compute the least-squares PC coefficients of the solution
%%%

% Solution PCE
p = 3; % Total order
P = factorial(p+d) / ( factorial(p) * factorial(d) ) - 1;
index_pc_soln = nD_polynomial_array(d, p);

% Loop over number of samples
n_samples = round([0.5, 0.9, 1.1, 4.0]*(P+1));
n_samples = [5:100];
mean = nan(1,length(n_samples));
var  = nan(1,length(n_samples));
for i = 1:length(n_samples)
    N = n_samples(i);
    Psi_mat = nan(N,P+1);
    u_vec = nan(N,1);
    for j = 1:N
        y = randn(1,2);
        Phi_i_soln = piset_hermite(y, index_pc_soln);
        u = Heat_Solution(dx, u0, uf, y(1), y(2), Ki, Pk, x, index_pc);
        Psi_mat(j,:) = Phi_i_soln;
        u_vec(j) = u(floor(length(u)/2));
    end
    c_vec = (Psi_mat' * Psi_mat) \ (Psi_mat' * u_vec);
    mean(i) = c_vec(1);
    var(i)  = sum(c_vec.^2) - mean(i)^2;
end

figure();

subplot(1,2,1);
loglog(n_samples, mean, 'DisplayName', 'mean');
% xlim([min(n_samples),max(n_samples)]);
xlabel('N');
legend('show');

subplot(1,2,2);
loglog(n_samples, var, 'DisplayName', 'var');
% xlim([min(n_samples),max(n_samples)]);
xlabel('N');
legend('show');

end





















