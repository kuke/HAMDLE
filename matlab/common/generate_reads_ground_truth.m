%%should be called when test dataset is changed
function generate_reads_ground_truth(READ_seq, REF_seq,rank)
% rank = 'Phylum', 'Class', 'Order', 'Family', 'Genus'
    key_index = 2;
    switch rank
        case 'Phylum'
            key_index = 2;
        case 'Class'
            key_index = 3;
        case 'Order'
            key_index = 4;
        case 'Family'
            key_index = 5;
        case 'Genus'
            key_index = 6;
    end
    rank_map = containers.Map();
    num_reference = size(REF_seq, 1);
    seq2species = zeros(num_reference, 1);
    key_id = 1;
    for n = 1:num_reference
        semicolon = strfind(REF_seq(n).Header,';');
        if(key_index > size(semicolon,2)) 
            continue;
        else if (key_index == size(semicolon,2))
            key_name = REF_seq(n).Header(semicolon(key_index)+1:end);
            else
                key_name = REF_seq(n).Header(semicolon(key_index)+1:semicolon(key_index+1)-1);
            end
        end
        
        quota = strfind(key_name, '"');
        if size(quota, 2) == 2
            key_name = key_name(quota(1)+1:quota(2)-1);
        end
        if ~isKey(rank_map, key_name)
            rank_map(key_name) = key_id;
            key_id = key_id+1;
        end
        seq2species(n) = rank_map(key_name);
    end
    
    disp('rank_map size');
    size(rank_map)
    num_fragments = size(READ_seq, 1);
    ref_seqs = zeros(num_fragments, 1);
    rank_keys = keys(rank_map);
    for n = 1:num_fragments
        header = READ_seq(n).Header;
        for i = 1:size(rank_keys, 2)
            if size(strfind(header,rank_keys{i}),1)>0
                ref_seqs(n) = rank_map(rank_keys{i});
                break;
            end
        end
%      min_reads = min(min_reads, size(READ_seq(n).Sequence));
%      start_ = strfind(header, 'reference=')+10;
%       end_ = strfind(header, ' position')-1;
%      header = header(start_:end_);
    % ref_seqs(n) = ref_id(header); 
    end
 
    
    NoOfSpecies = size(rank_map,1);
    sol_species = zeros(NoOfSpecies, 1);
    for n = 1:num_fragments
        header = READ_seq(n).Header;
        for i = 1:size(rank_keys,2)
            if size(strfind(header,rank_keys{i}),1)>0
                sol_species(rank_map(rank_keys{i})) = sol_species(rank_map(rank_keys{i}))+1;
                break;
            end
        end
    end
    
    fragments_ref_map = containers.Map();
    for n = 1:num_fragments
        header = READ_seq(n).Header;
        start_ = strfind(header, 'reference=')+10;
        end_ = strfind(header, ' position')-1;
        header = header(start_:end_);
        if ~isKey(fragments_ref_map, header)
            fragments_ref_map(header) = 1;
        else
            fragments_ref_map(header) = fragments_ref_map(header)+1;
        end
    end
 
    sol_seqs = zeros(1, num_reference);
    fragments_keys = keys(fragments_ref_map);
    for n =1:num_reference
        header = REF_seq(n).Header;
        for i = 1:size(fragments_keys,2)
            if size(strfind(header, fragments_keys{i}),1)>0
                sol_seqs(n) = fragments_ref_map(fragments_keys{i});
                break;
            end
        end
    end

 filename = strcat('LargeGroundTruth.mat');
 save(filename,'rank_map','ref_seqs','seq2species','sol_seqs','sol_species','NoOfSpecies');

end