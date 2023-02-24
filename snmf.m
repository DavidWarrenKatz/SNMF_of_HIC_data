function snmf()
%
% SNMF_of_HIC_data
% This code inputs a .hic file and an positive integer n. It calculates the optimal symmetric
% non-negative rank r approximation of HiC contact
% matrix for all r less than or equal to n. The results are displayed as a r(x-axis)-by-error of
% optimal-rank-r approximation (y-axis) graph.

    clc;
    clear;
    close all;

% hard-coded .hic file to practice with. This command returns a hicstraw.HiCFile object. This file contains hic data with many resolutions.
    hic = py.hicstraw.HiCFile("https://www.encodeproject.org/files/ENCFF718AWL/@@download/ENCFF718AWL.hic")
    
    matrix_object_chr4 = py.hic.getMatrixZoomData('4', '4', "observed", "NONE", "BP", 5000)n = py.int(input("Input an integer n: "))

    %% Initialize of rank to be factorized
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
    
end
