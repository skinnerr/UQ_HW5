function [] = Problem_1c()

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
% Compute the "true" solution based on sparse sampling
%%%

[m_true, v_true, ~] = Sample_Sparse(dx, u0, uf, Ki, Pk, x, index_pc);

% figure();
% subplot(2,1,1)
% plot(x, [nan, m_true, nan])
% xlabel('x');
% ylabel('mean');
% subplot(2,1,2)
% plot(x, [nan, v_true, nan])
% xlabel('x');
% ylabel('variance');

%%%
% Compute the least-squares PC coefficients of the solution
%%%

% Solution PCE
p_list = 1:4; % Total order
n_samples = nan(1,length(p_list));

% Loop over number of samples
mean = zeros(1,length(p_list));
var  = zeros(1,length(p_list));
for i = 1:length(p_list)
    p = p_list(i);
    P = factorial(p+d) / ( factorial(p) * factorial(d) ) - 1;
    index_pc_soln = nD_polynomial_array(d, p);
    N = round(4*(P+1));
    n_samples(i) = N;
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
var_relerr = abs(var - v_true(floor(Nx/2))) / v_true(floor(Nx/2));
semilogy(p_list, var_relerr, 'o-', 'DisplayName', 'var. rel. err.');
% xlim([min(n_samples),max(n_samples)]);
xlabel('p');
legend('show');

end





















