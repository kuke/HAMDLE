%% compute within- and between- class max cosine similarity
addpath('common');
load('Legendre_coeff_vector_allref_trans_Legendre.mat');
[len_feat, num_fragments] = size(all_species_Legendre_coeff_trans);
num_species = 21;
index_species = find((fragment2species-[0,fragment2species(1:end-1)])==1);
offset_species = zeros(1, num_species);
for i = 1:num_species-1
    offset_species(i) = index_species(i+1)-index_species(i)-1;
end

offset_species(num_species) = size(fragment2species,2)-index_species(num_species);
max_cos_similarity = zeros(num_species, num_species);
sum(offset_species)
for i = 1:num_species
    max_cos_similarity(i, i) = 1;
end

for i = 1:num_species
    for j = i+1:num_species
        max_cos_similarity(i, j) = compute_between_class_max_cosine_similarity( ...
                    all_species_Legendre_coeff_trans(:, index_species(i):index_species(i)+offset_species(i)), ...
                    all_species_Legendre_coeff_trans(:, index_species(j):index_species(j)+offset_species(j)));
        max_cos_similarity(j, i) = max_cos_similarity(i, j);
    end
end

disp('similarity = ')
max_cos_similarity
colormap('hot');
imagesc(max_cos_similarity)
colorbar
xlabel('species'); ylabel('species');
set(gca,'FontSize',18);
set(gca,'XTick',0:5:21);
set(gca,'YTick',0:5:21);


