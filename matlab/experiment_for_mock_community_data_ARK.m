%% ARK (Aggregation of Reads by K-means) framework algorithm

% ARK framework will be used for HAMDLE
% Adapted from ARK paper---"ARK: Aggregation of Reads by K-Means for
% Estimation of Bacterial Community Composition"

%% ARK method 

% For experimental evaluation, Saikat Chatterjee used the mock communities databse of 16S rRNA, 
% as used in BeBAC and Legendre papers (both papers are in Bioinformatics journal)

clear all; clc; clear;
addpath('common');

%% variable definitions

data = '../data/mock_community/';
%k = 4;  % k for k-mers
NoOfSpecies = 21; % the database contains 21 species for references and reads

% Variables for Legendre based on OMP^{+,1} (please see Algorithm 1 in Legendre paper)
I=1000; % Maximum allowable iteration in Legendre:OMP^{+,1}    
nu=0.001; % allowable tolerance in l1 norm computation for OMP^{+,1}

% Variable for Quikr
lambda = 10000; 

% Variable for clustering (LBG algorithm based K-means clustering) in ARK algorithm
eta=0.0005; % User choice: Tolerance till convergence of squared error distance (for all algorithms: Legendre, Quikr and Taxy).
%MaxNoOfClusters = 32; % Q_{max}    User choice: Maximum allowable number of clusters 
MaxNoOfClusters = 2;

%% Loading GroundTruth variable which contains ground truth and reference species identification
load([data, 'GroundTruth.mat']);
true_solution=sol_species'; 
true_solution=true_solution/sum(true_solution);

%% Set parameters
L = 450;
X = -1:2/(L-1):1;
order = 49;
Legendre_mat = generate_Legendre_matrix(order, X);

%% System Matrix Generation (by processing reference data followed by computing k-mers)
% Note that this is off-line computation

% System Matrix Generation for Legendre
% Legendre requires more time to compute system matrix. Hence, if the directory contains system matrix then skip computation (automatically skips this step if necessary .mat files are present).

REF_seq = fastaread([data, 'Reference.fasta']);
generate_ref_coeff_vector_Legendre(REF_seq, Legendre_mat, order); % saved as kmer_vector_allref_Legendre.mat
allref_vector_compile_Legendre('Legendre_coeff_vector_allref.mat', seq2species); % saved as kmer_vector_allref_trans_Legendre.mat
load('Legendre_coeff_vector_allref_trans_Legendre.mat');
X = all_species_Legendre_coeff_trans; % X is the system matrix for Legendre


%% Legendre vector from reads (Online computation)
% Note that all three methods (Legendre, Quikr and Taxy) use same k-mers vectors. Hence one program function is sufficient
% Further, if k-mers from reads are already pre-computed and exist in the directory then skip computation to save time (automatically skips this step if necessary .mat files are present).

READ_seq = fastaread([data,'Measurement.fasta']);
generate_read_coeff_vector_Legendre(READ_seq, Legendre_mat, order); % saved as kmer_vector_reads_Legendre.mat
load('Legendre_coeff_vector_reads.mat'); 


%% ARK HAMDLE simulation (ARK_HAMDLE)

% Initialization
Composition_ARK_HAMDLE  = [];
ChangeInComposition_ARK_HAMDLE = 1;
NoOfClusters_Legendre = 0;
data_ARK_HAMDLE = [];
tstart_ARK_HAMDLE = tic;

while (ChangeInComposition_ARK_HAMDLE  > eta) && (NoOfClusters_Legendre < MaxNoOfClusters)  % (stopping criteria for LBG based clustering)
    % ---- This is LBG based K-means clustering ----
    % The variable C_ARK_HAMDLE contains mean vectors of clusters
    if NoOfClusters_Legendre == 0
        C_ARK_HAMDLE = mean(Legendre_coeff_vector_reads); ClusterProbability = [1];
    else
        [C_ARK_HAMDLE, ClusterProbability] = LBG(Legendre_coeff_vector_reads, C_ARK_HAMDLE, ClusterProbability);  % LBG algorithm increases the number of clusters as output from the input no of clusters by one 
    end
    
    NoOfClusters_Legendre = length(ClusterProbability);

    % After clustering, HAMDLE:OMP^{+,1} is used for each cluster
    Mu_ARK_HAMDLE = C_ARK_HAMDLE';  % Cluster mean vectors)
    tmp = 0; gamma = 0;
    result_ARK_HAMDLE = zeros(1,NoOfSpecies);
    
    for i=1:NoOfClusters_Legendre
        [tmp, ~]= OMP_plus_1_for_HAMDLE(X,Mu_ARK_HAMDLE(:,i),nu,I);
        tmp_ARK_HAMDLE = zeros(1,NoOfSpecies);
        for j=1:length(tmp)
            if tmp(j) ~=0
                tmp_ARK_HAMDLE(fragment2species(j)) = tmp_ARK_HAMDLE(fragment2species(j)) + tmp(j);
            end
        end
        
        % result_ARK_HAMDLE contains final estimate of community composition by linear addition       
        result_ARK_HAMDLE = result_ARK_HAMDLE + ClusterProbability(i)*tmp_ARK_HAMDLE;
    end

    VD_ARK_HAMDLE = 0.5 * norm((true_solution - result_ARK_HAMDLE), 1);
    fprintf('VD_ARK_HAMDLE: %f\n', VD_ARK_HAMDLE );
    
    if NoOfClusters_Legendre > 1, ChangeInComposition_ARK_HAMDLE  = norm((Composition_ARK_HAMDLE(end,:) - result_ARK_HAMDLE), 1);    end
    
    fprintf('VD with in current and previous iterations at %u clusters: %d\n', NoOfClusters_Legendre, ChangeInComposition_ARK_HAMDLE)
    
    data_ARK_HAMDLE = [data_ARK_HAMDLE; NoOfClusters_Legendre ChangeInComposition_ARK_HAMDLE VD_ARK_HAMDLE];
    
    Composition_ARK_HAMDLE  = [Composition_ARK_HAMDLE; result_ARK_HAMDLE];
     
    pause(1);
end

elapsedtime_ARK_HAMDLE=toc(tstart_ARK_HAMDLE);

fprintf('Elapsed time till convergence: %f\n',  elapsedtime_ARK_HAMDLE);
fprintf('Number of clusters at convergence: %f\n', NoOfClusters_Legendre); 

%% Saving outputs

save('Convergencedata.mat','elapsedtime_ARK_HAMDLE','data_ARK_HAMDLE','Composition_ARK_HAMDLE');

disp ('---------------------------------');
fprintf('Elapsed time till convergence: %f\n', elapsedtime_ARK_HAMDLE);
figure;
plot(data_ARK_HAMDLE(:,1), data_ARK_HAMDLE(:,3), 'b', 'linewidth', 1.5); xlabel('Number of clusters'); ylabel('VD'); 
set(gca,'FontSize',15);
frame = getframe(gcf);
imwrite(frame.cdata, 'VD_vs_Q_result.jpg');
savefig('VD_vs_Q_result.fig');

%% Showing bar chart
species = {'A.baumannii', 'A.odontolyticus', 'B.cereus', 'B.vulgatus', 'C.beijerinckii', ...
            'D.radiodurans', 'E.coli', 'E.faecalis', 'H.pylori', 'L.gasseri', ...
            'L.monocytogenes', 'M.smithii', 'N.meningitidis', 'P.acnes', 'P.aeruginosa', ...
            'R.sphaeroides', 'S.agalactiae', 'S.aureus', 'S.epidermidis', 'S.mutans', 'S.pneumoniae'};
        
[VD_min, index_min]=min(data_ARK_HAMDLE(:,3))
BarYY = [true_solution; Composition_ARK_HAMDLE(1,:); Composition_ARK_HAMDLE(index_min,:)];
BarX = 1 : NoOfSpecies;
figure;
axes1 = axes('Position',[0.0738 0.198 0.894 0.771]);
box(axes1,'on');
hold(axes1,'on');

bar(BarX,BarYY',1);
legend('Ground Truth', 'HAMDLE', 'ARK HAMDLE', 'Location','North');
axis([0 NoOfSpecies+1 0 0.3]);
ylabel('proportion');
set(gca,'FontSize',15);
set(gcf,'Position',[360,278,840,420], 'color','w');

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
imwrite(frame.cdata, 'ARK_HAMDLE_result.jpg');
savefig('ARK_HAMDLE_result.fig');


