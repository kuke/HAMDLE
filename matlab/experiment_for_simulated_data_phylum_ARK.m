%% ARK (Aggregation of Reads by K-means) framework algorithm

% ARK framework will be used for HAMDLE for the simulated data
% Adapted from ARK paper---"ARK: Aggregation of Reads by K-Means for
% Estimation of Bacterial Community Composition"

%% ARK method 

% For experimental evaluation, Saikat Chatterjee used the mock communities databse of 16S rRNA, 
% as used in BeBAC and Legendre papers (both papers are in Bioinformatics journal)

clear all; clc; clear;
addpath('common');

%% variable definitions
data = '../data/simulated/';
need_reprocess_data = 0; % =1 when it is first time to run or parameters changed, otherwise = 0
NoOfSpecies = 39; % the database contains 39 phylum for references and reads

% Variables for Legendre based on OMP^{+,1} (please see Algorithm 1 in HAMDLE paper)
I=1000; % Maximum allowable iteration in Legendre:OMP^{+,1}    
nu=0.00001; % allowable tolerance in l1 norm computation for OMP^{+,1}

% Variable for clustering (LBG algorithm based K-means clustering) in ARK algorithm
eta=0.0005; % User choice: Tolerance till convergence of squared error distance (for all algorithms: Legendre, Quikr and Taxy).
%MaxNoOfClusters = 32; % Q_{max}    User choice: Maximum allowable number of clusters 
MaxNoOfClusters = 2;
rank = 'Phylum';

%% Loading GroundTruth variable which contains ground truth and reference species identification
load([data, rank, 'GroundTruth.mat']);
true_solution=sol_species'; 
true_solution=true_solution/sum(true_solution);

%% Set parameters
L = 400;
Lp = 100;
X = -1:2/(L-1):1;
order = 49;
Legendre_mat = generate_Legendre_matrix(order, X);

%% System Matrix Generation (by processing reference data followed by computing k-mers)
% Note that this is off-line computation

% System Matrix Generation for Legendre
% Legendre requires more time to compute system matrix. Hence, if the directory contains system matrix then skip computation (automatically skips this step if necessary .mat files are present).

if need_reprocess_data == 1
    REF_seq = fastaread([data, 'trainset7_112011.fa'],'TRIMHEADERS', false);
    generate_ref_coeff_vector_Legendre(REF_seq, Legendre_mat, order, Lp); % saved as kmer_vector_allref_Legendre.mat
    allref_vector_compile_Legendre('Legendre_coeff_vector_allref.mat', seq2species); % saved as kmer_vector_allref_trans_Legendre.mat
end
load('Legendre_coeff_vector_allref_trans_Legendre.mat');
X = all_species_Legendre_coeff_trans; % X is the system matrix for Legendre


%% Legendre vector from reads (Online computation)
if need_reprocess_data == 1
    READ_seq = fastaread([data, 'grinder-1000-reads.fa']);
    generate_read_coeff_vector_Legendre(READ_seq, Legendre_mat, order); % saved as kmer_vector_reads_Legendre.mat
end
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
xlim([1 MaxNoOfClusters]);
frame = getframe(gcf);
imwrite(frame.cdata, 'VD_vs_Q_simulated_phylum_result.jpg');
savefig('VD_vs_Q_simulated_phylum_result.jpg.fig');

%% Showing bar chart
result_SEK = load([data, 'phylum_result_SEK.txt']);       
[VD_min, index_min]=min(data_ARK_HAMDLE(:,3));
fprintf('VD_min = %f, index_min = %d\n', VD_min, index_min);
BarYY = [true_solution; result_SEK; Composition_ARK_HAMDLE(1,:); Composition_ARK_HAMDLE(index_min,:)];
BarX = 1 : NoOfSpecies;
figure;
axes1 = axes('Position',[0.074 0.117 0.89 0.85]);
box(axes1,'on');
hold(axes1,'on');

bar(BarX,BarYY',1);
legend('Ground Truth', 'SEK', 'HAMDLE', 'ARK HAMDLE', 'Location','North');
axis([0 NoOfSpecies+1 0 0.35]);
xlabel('Phylum'); ylabel('proportion');
set(gca,'FontSize',15);
set(gcf,'Position',[360,278,840,420], 'color','w');

% label species in X axis
frame = getframe(gcf);
imwrite(frame.cdata, ['ARK_HAMDLE_simulated_phylum_result', '_49', '.jpg']);
savefig(['ARK_HAMDLE_simulated_phylum_result','_49','.fig']);


