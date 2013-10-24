%% Estimando o Modelo de Heston

load OptionData

% initial guess, lower bound and upper bound for the parameters.
x0 = [6.5482 2.3012 -0.4176 0.0731 0.1838];
lb = [0 0 -1 0 0];
ub = [20 5 1 1 1];

HD_handler = HestonDifferences(@CallHeston, OptionData, true);

options = optimset('MaxFunEvals',20000);
HestonParameters = lsqnonlin(HD_handler, x0, lb, ub, options);


%% Comparando realizado x esperado

lambda = HestonParameters(1);
eta = HestonParameters(2);
rho = HestonParameters(3);
vbar = HestonParameters(4);
v0 = HestonParameters(5);

C_h = zeros(size(OptionData, 1), 1);
for i = 1:size(OptionData, 1)
    C_h(i) = CallHeston(OptionData(i, 1), OptionData(i, 3), OptionData(i, 2), 0, [lambda, eta, rho, vbar, v0]);
end

T = OptionData(:, 2);
K = OptionData(:, 3);
C_m = OptionData(:, 4);

plot3(K, T, C_h, 'ok');
hold on;
plot3(K, T, C_m, '*r');
grid on

figure(2);
clf;
stem3(K, T, C_m - C_h);

