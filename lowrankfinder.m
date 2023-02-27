
% obtain matrix of hicdata as 'V' from python straw library
pyrunfile("obtainHiCdata.py")
load('matlabfile.mat')


    clc;
    clear;
    close all;

    rank = 5;

    %% perform factroization
    options.verbose = 1;

    % MU
    options.alg = 'mu';
    [w_mu, infos_mu] = fro_mu_nmf(V, rank, options);

    % Hierarchical ALS
    options.alg = 'hals';
    [w_hals, infos_hals] = als_nmf(V, rank, options);

    % Accelerated Hierarchical ALS
    options.alg = 'acc_hals';
    [w_acchals, infos_acchals] = als_nmf(V, rank, options);


    %% plot
    display_graph('epoch','cost', {'Fro-MU', 'HALS', 'Acc-HALS'}, {w_mu, w_hals, w_acchals}, {infos_mu, infos_hals, infos_acchals});
    display_graph('time','cost', {'Fro-MU', 'HALS', 'Acc-HALS'}, {w_mu, w_hals, w_acchals}, {infos_mu, infos_hals, infos_acchals});



n = input('please enter an integer : ');

for i = 1:n
     display(i);
end
