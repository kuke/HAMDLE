%%
% This file produces results for Fig. 2 in paper "SEK: Sparsity exploiting
% k -mer-based estimation of bacterial community composition"

close all; clc; clear;

%%
%%%% Necessary files %%%%

% Sequence data
data = '../data/mock_community/';
READ_seq = fastaread([data, 'Measurement.fasta']);
REF_seq = fastaread([data, 'Reference.fasta']);

%%
% Varibale definition
NoOfSpecies = 21; % Number of species in the reference database
min_order = 49;
max_order = 49;

L = 450;
X = -1: 2/(L-1):1;
%Legendre matrix
Legendre_mat = generate_Legendre_matrix(max_order, X);

start_time_str = datestr(now,30);
log_dir = '../log';
outputfig_dir = '../outputfig';
mkdir([outputfig_dir,'/',start_time_str]);
fid = fopen([log_dir,'/',start_time_str,'.log'],'w');

fprintf(fid, 'Data source: %s\n', data);
fprintf(fid, 'L = %d\n', L);
fprintf(fid, 'Min_order = %d, Max_order = %d\n\n', min_order, max_order);
fprintf(fid, 'Order\t Elapsed_READ\t Elapsed_REF\t Elapsed_ALGO\t VD_HAMDLE\t It\n');

% loading ground truth for comparison
load([data, 'GroundTruth.mat']);
true_solution=sol_species';  % sol_species comes from Ground Truth (GroundTruth.mat)
true_solution=true_solution/sum(true_solution);

for order = min_order:max_order
    close all; clc;
    fprintf('current order = %d\n', order);
    %%
    % SEK implementation (using OMP^{+,1} algorithm
    %[result_SEK elapsedtime_READ elapsedtime_REF elapsedtime_algo] = SEK_implementation(k_SEK,READ_seq,REF_seq,NoOfSpecies);
    [result_Legendre elapsedtime_READ_Legendre elapsedtime_REF_Legendre elapsedtime_algo_Legendre, it] = Legendre_implementation(order,Legendre_mat(1:order+1,:), READ_seq,REF_seq, NoOfSpecies, seq2species);

    VD_HAMDLE = 0.5 * norm((true_solution - result_Legendre),1);

    fprintf('Variational distance performance: VD_HAMDLE = %f\n', VD_HAMDLE); 
    
    fprintf(fid, '%5g\t%10g\t%10g\t%10g\t%10g\t%5g\t\n',order,elapsedtime_READ_Legendre, ...
                                        elapsedtime_REF_Legendre, ...
                                        elapsedtime_algo_Legendre, ...
                                        VD_HAMDLE,it);
    
    figure
    % Showing bar chart
    result_SEK = [ 0.0828, 0.0144, 0.0076, 0.0899, 0.0497, 0.1989, 0.0121, 0.0095, 0.0041, 0.0062, 0.0210, ...
    0.0022, 0.0305, 0.0420, 0.0267, 0.0039, 0.0252, 0.2700, 0.0500, 0.0241, 0.0293];
    BarY = [true_solution; result_SEK; result_Legendre];

    BarX = 1:21;
    bar(BarX,BarY', 1)
    legend('Ground Truth', 'SEK', 'HAMDLE','Location','North');
    axis([0 22 0 0.3]);
    xlabel('species'); ylabel('proportion');
    set(gca,'FontSize',18);
 %   log_dir = './log';
%outputfig_dir = './outputfig';
    saveas(gcf,[outputfig_dir,'/',start_time_str,'/Order_',num2str(order)],'jpg');
end
fclose(fid);
