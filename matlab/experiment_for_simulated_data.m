%%
% This file produces results for Fig. 2 in paper "SEK: Sparsity exploiting
% k -mer-based estimation of bacterial community composition"

close all; clc; clear;

%%
%%%% Necessary files %%%%

% Sequence data
data = '../data/simulated/';
READ_seq = fastaread([data, 'grinder-1000-reads.fa']);
REF_seq = fastaread([data, 'trainset7_112011.fa'],'TRIMHEADERS', false);

%%
% Varibale definition
min_order = 49;
max_order = 49;

L = 400;
X = -1: 2/(L-1):1;
%Legendre matrix
Legendre_mat = generate_Legendre_matrix(max_order, X);

start_time_str = datestr(now,30);
log_dir = './log';
outputfig_dir = './outputfig';
mkdir([outputfig_dir,'/',start_time_str]);
fid = fopen([log_dir,'/',start_time_str,'.log'],'w');

% rank = Phylum, Class, Order, Family, Genus
rank = 'Genus';
fprintf(fid, 'Data source: %s\n', data);
fprintf(fid, 'L = %d, rank = %s, diversity = %d\n', L, rank, 1000);
fprintf(fid, 'min_order = %d, max_order = %d\n\n', min_order, max_order);
fprintf(fid, 'Order\t Elapsed_READ\t Elapsed_REF\t Elapsed_algo\t VD_Legendre\t It\n');

generate_reads_ground_truth(READ_seq,REF_seq,rank);
% loading ground truth for comparison
load('LargeGroundTruth.mat');
for order = min_order:max_order
    
    close all; clc;
    
    disp('current order = ');order
    %%
    % SEK implementation (using OMP^{+,1} algorithm
    %[result_SEK elapsedtime_READ elapsedtime_REF elapsedtime_algo] = SEK_implementation(k_SEK,READ_seq,REF_seq,NoOfSpecies);
    [result_Legendre elapsedtime_READ_Legendre elapsedtime_REF_Legendre elapsedtime_algo_Legendre, it] = HAMDLE_implementation(order,Legendre_mat(1:order+1,:), READ_seq,REF_seq, NoOfSpecies);


    %%
    % Comparison



    true_solution=sol_species';  % sol_species comes from Ground Truth (GroundTruth.mat)
    true_solution=true_solution/sum(true_solution);
    disp('some size');
    size(true_solution)
    size(result_Legendre)
    VD_Legendre = 0.5 * norm((true_solution - result_Legendre),1);

    disp('Variational distance performance:'); VD_Legendre
    
    fprintf(fid, '%5g\t%10g\t%10g\t%10g\t%10g\t%5g\t\n',order,elapsedtime_READ_Legendre, ...
                                        elapsedtime_REF_Legendre, ...
                                        elapsedtime_algo_Legendre, ...
                                        VD_Legendre,it);
    
    figure
    % Showing bar chart
    BarY = [true_solution; result_Legendre];
    BarX = 1:NoOfSpecies;
    bar(BarX,BarY',1)
    legend('Ground Truth', 'Legendre:OMP^{+,1}','Location','North');
    axis([0 NoOfSpecies+1 0 0.3]);
    xlabel(rank); ylabel('proportion');
    set(gca,'FontSize',18);
    
 %   log_dir = './log';
%outputfig_dir = './outputfig';
    saveas(gcf,[outputfig_dir,'/',start_time_str,'/order_',num2str(order)],'jpg');
end
fclose(fid);
