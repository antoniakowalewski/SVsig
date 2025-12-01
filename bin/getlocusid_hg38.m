function locusid=getlocusid_hg38(chr,pos,refgene_lookup,pad_high,pad_low)

chr_col   = cell2mat(refgene_lookup(:,1));  % chromosome
start_col = cell2mat(refgene_lookup(:,2));  % gene start
end_col   = cell2mat(refgene_lookup(:,3));  % gene end
loc=find(chr_col == chr & start_col < pos & end_col > pos);
if ~isempty(loc)
    locusid=refgene_lookup{loc(1),4};
else
    loc_neg = find(chr_col == chr & end_col < pos & end_col > pos - pad_low);
    loc_pos = find(chr_col == chr & start_col > pos & start_col < pos + pad_high);
    if isempty(loc_neg) && ~isempty(loc_pos)
        locusid = refgene_lookup{loc_pos(1), 4};
    elseif ~isempty(loc_neg) && isempty(loc_pos)
        locusid = refgene_lookup{loc_neg(1), 4};
    elseif ~isempty(loc_neg) && ~isempty(loc_pos)
        dist_neg = abs(end_col(loc_neg(1)) - pos);
        dist_pos = abs(start_col(loc_pos(1)) - pos);
        if dist_neg > dist_pos
            locusid = refgene_lookup{loc_neg(1),4};
        else
            locusid = refgene_lookup{loc_pos(1),4};
        end
    else
        locusid = -1e10;
    end
end