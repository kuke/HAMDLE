%%
% This file produces results for Fig. 2 in paper "SEK: Sparsity exploiting
% k -mer-based estimation of bacterial community composition"

close all; clc; clear;
addpath('common');
%%
%%%% Necessary files %%%%

% Sequence data
data = '../data/simulated/';
READ_seq = fastaread([data, 'grinder-1000-reads.fa']);
REF_seq = fastaread([data, 'trainset7_112011.fa'],'TRIMHEADERS', false);

%%
% Varibale definition
min_order = 5;
max_order = 5;

L = 400;
Lp = 100;
X = -1: 2/(L-1):1;
%Legendre matrix
Legendre_mat = generate_Legendre_matrix(max_order, X);

start_time_str = datestr(now,30);
log_dir = '../log';
outputfig_dir = '../outputfig';
mkdir([outputfig_dir,'/',start_time_str]);
fid = fopen([log_dir,'/',start_time_str,'.log'],'w');

% rank = Phylum, Class, Order, Family, Genus
rank = 'Phylum';
fprintf(fid, 'Data source: %s\n', data);
fprintf(fid, 'L = %d, rank = %s, diversity = %d\n', L, rank, 1000);
fprintf(fid, 'min_order = %d, max_order = %d\n\n', min_order, max_order);
fprintf(fid, 'Order\t Elapsed_READ\t Elapsed_REF\t Elapsed_algo\t VD_HAMDLE\t It\n');

if ~exist([data, rank, 'GroundTruth.mat'], 'file')
    generate_reads_ground_truth(READ_seq,REF_seq,rank, data); 
end
% loading ground truth 
load([data, rank, 'GroundTruth.mat']);
true_solution=sol_species';  % sol_species comes from Ground Truth (GroundTruth.mat)
true_solution=true_solution/sum(true_solution);
disp('some size');
size(true_solution)
    
for order = min_order:max_order
    close all; clc;
    
    fprintf('current order = %d\n', order);
    %%
    % SEK implementation (using OMP^{+,1} algorithm
    %[result_SEK elapsedtime_READ elapsedtime_REF elapsedtime_algo] = SEK_implementation(k_SEK,READ_seq,REF_seq,NoOfSpecies);
    [result_Legendre elapsedtime_READ_Legendre elapsedtime_REF_Legendre elapsedtime_algo_Legendre, it] = HAMDLE_implementation(order,Legendre_mat(1:order+1,:), ...
                                                                 READ_seq,REF_seq, NoOfSpecies, seq2species, Lp);


    %%
    VD_HAMDLE = 0.5 * norm((true_solution - result_Legendre),1);

    fprintf('Variational distance performance: %f\n',  VD_HAMDLE);
    
    fprintf(fid, '%5g\t%10g\t%10g\t%10g\t%10g\t%5g\t\n',order,elapsedtime_READ_Legendre, ...
                                        elapsedtime_REF_Legendre, ...
                                        elapsedtime_algo_Legendre, ...
                                        VD_HAMDLE,it);
    
    figure
    % Showing bar chart
    axes1 = axes('Position',[0.0738 0.126 0.894 0.84]);
    box(axes1,'on');
    hold(axes1,'on');
    BarY = [true_solution; result_Legendre];
    BarX = 1:NoOfSpecies;
    bar(BarX,BarY',1)
    legend('Ground Truth', 'HAMDLE','Location','North');
    axis([0 NoOfSpecies+1 0 0.35]);
    xlabel(rank); ylabel('proportion');
    set(gca,'FontSize',15);
    set (gcf,'Position',[360,278,840,420], 'color','w');

    frame = getframe(gcf);
    imwrite(frame.cdata, [outputfig_dir,'/',start_time_str,'/Order_',num2str(order),'.jpg']);
end
fclose(fid);
