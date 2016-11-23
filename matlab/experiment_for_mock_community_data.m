% Basic experiment for the mock community data, HAMDLE
% 
close all; clc; clear;
addpath('common');

%% Import sequence data
data = '../data/mock_community/';
READ_seq = fastaread([data, 'Measurement.fasta']);
REF_seq = fastaread([data, 'Reference.fasta']);

%% Parameters 
NoOfSpecies = 21;   % Number of species in the reference database
min_order = 5;     % Minimun order of Legendre expansion
max_order = 5;     % Maximun order of Legendre expansion
L = 450;            % Length of sliding window

%% Generate Legendre matrix
X = -1: 2/(L-1):1;
Legendre_mat = generate_Legendre_matrix(max_order, X);

%% Output setup
start_time_str = datestr(now,30);
log_dir = '../log';
outputfig_dir = '../outputfig';
mkdir([outputfig_dir,'/',start_time_str]);
fid = fopen([log_dir,'/',start_time_str,'.log'],'w');

fprintf(fid, 'Data source: %s\n', data);
fprintf(fid, 'Length of sliding window = %d\n', L);
fprintf(fid, 'Min_order = %d, Max_order = %d\n\n', min_order, max_order);
fprintf(fid, 'Order\t Elapsed_READ\t Elapsed_REF\t Elapsed_ALGO\t VD_HAMDLE\t It\n');

%% Load ground truth of reads data
load([data, 'GroundTruth.mat']);
true_solution=sol_species';  % sol_species comes from Ground Truth (GroundTruth.mat)
true_solution=true_solution/sum(true_solution);

%% Execute 
for order = min_order:max_order
    close all; clc;
    fprintf('current order = %d\n', order);
    % HAMDLE implementation
    [result_Legendre elapsedtime_READ_Legendre elapsedtime_REF_Legendre elapsedtime_algo_Legendre, it] = Legendre_implementation(order,Legendre_mat(1:order+1,:), READ_seq,REF_seq, NoOfSpecies, seq2species);

    VD_HAMDLE = 0.5 * norm((true_solution - result_Legendre), 1);

    fprintf('Variational distance performance: VD_HAMDLE = %f\n', VD_HAMDLE); 
    
    fprintf(fid, '%5g\t%10g\t%10g\t%10g\t%10g\t%5g\t\n',order,elapsedtime_READ_Legendre, ...
                                        elapsedtime_REF_Legendre, ...
                                        elapsedtime_algo_Legendre, ...
                                        VD_HAMDLE,it);
    
    figure
    
    result_SEK = [ 0.0828, 0.0144, 0.0076, 0.0899, 0.0497, 0.1989, 0.0121, 0.0095, 0.0041, 0.0062, 0.0210, ...
    0.0022, 0.0305, 0.0420, 0.0267, 0.0039, 0.0252, 0.2700, 0.0500, 0.0241, 0.0293];
    species = {'A.baumannii', 'A.odontolyticus', 'B.cereus', 'B.vulgatus', 'C.beijerinckii', ...
            'D.radiodurans', 'E.coli', 'E.faecalis', 'H.pylori', 'L.gasseri', ...
            'L.monocytogenes', 'M.smithii', 'N.meningitidis', 'P.acnes', 'P.aeruginosa', ...
            'R.sphaeroides', 'S.agalactiae', 'S.aureus', 'S.epidermidis', 'S.mutans', 'S.pneumoniae'};
    
    axes1 = axes('Position',[0.0738 0.198 0.894 0.771]);
    box(axes1,'on');
    hold(axes1,'on');
    % Showing bar chart
    BarY = [true_solution; result_SEK; result_Legendre];
    BarX = 1:NoOfSpecies;
    bar(BarX,BarY', 1)
    legend('Ground Truth', 'SEK', 'HAMDLE','Location','North');

    axis([0 NoOfSpecies+1 0 0.3]);
    ylabel('proportion');
    set(gca,'FontSize',15);
    set (gcf,'Position',[360,278,840,420], 'color','w');
    
    % label species in X axis
    set(gca, 'XTick', 1:length(species), 'XTickLabel', species);
    xtb = get(gca, 'XTickLabel');
    xt = get(gca, 'XTick');
    yt = get(gca, 'YTick');
    xtextp = xt;
    ytextp = yt(1)*ones(1, length(xt));
    text(xtextp, ytextp-0.01, xtb, 'HorizontalAlignment', 'right', 'rotation', 45, 'fontsize', 12);
    set(gca,'xticklabel','');
    frame = getframe(gcf);
    imwrite(frame.cdata, [outputfig_dir,'/',start_time_str,'/Order_',num2str(order),'.jpg']);
end
fclose(fid);
