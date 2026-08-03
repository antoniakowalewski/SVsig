%%%Function to renormalize tiles following event ratios. 
%Used at the end of break_invasion_copy.m and double_break_join_model.m}

function [normalized_mat] = renormalize_tiles(mat_ratios, mat, events, bins, CHR)

%ic ratio is ratio of inter:intra events
%intra ratio is ratio of short:long events 

annot_tiles=tiles_annot_copy('length',events,bins,CHR);

% Verify dimensions before logical indexing
if size(mat_ratios,1) ~= size(annot_tiles,1) || ...
   size(mat_ratios,2) ~= size(annot_tiles,2) || ...
   size(mat,1) ~= size(annot_tiles,1) || ...
   size(mat,2) ~= size(annot_tiles,2)

    error('renormalize_tiles:SizeMismatch', ...
        ['mat_ratios is %dx%d, mat is %dx%d, and annot_tiles is %dx%d. ' ...
         'All must have matching dimensions.'], ...
        size(mat_ratios,1), size(mat_ratios,2), ...
        size(mat,1), size(mat,2), ...
        size(annot_tiles,1), size(annot_tiles,2));
end

ratio_total = full(sum(mat_ratios(:)));
mat_total = full(sum(mat(:)));

if ratio_total <= 0 || ~isfinite(ratio_total)
    error('SVsig:InsufficientObservedSignal', ...
        ['Unable to estimate the background model because the observed event matrix ' ...
         'has zero or non-finite total mass (%g). This typically means that too few ' ...
         'structural variants remain after filtering to estimate the required event ' ...
         'categories. Consider using a precomputed background model or relaxing filters.'], ...
        ratio_total);
end

if mat_total <= 0 || ~isfinite(mat_total)
    error('SVsig:InsufficientModelSignal', ...
        ['Unable to estimate the background model because the expected probability ' ...
         'matrix has zero or non-finite total mass (%g). This usually indicates that ' ...
         'one or more event categories cannot be estimated from the remaining cohort data. ' ...
         'Consider using a precomputed background model or relaxing filters.'], ...
        mat_total);
end

% Desired ratios from observed events
diag_short_ratio = full(sum(mat_ratios(annot_tiles(:,:,1)))) / ratio_total;
short_ratio      = full(sum(mat_ratios(annot_tiles(:,:,2)))) / ratio_total;
long_ratio       = full(sum(mat_ratios(annot_tiles(:,:,3)))) / ratio_total;
inter_ratio      = full(sum(mat_ratios(annot_tiles(:,:,4)))) / ratio_total;

% Current ratios in model
diag_short_annot = full(sum(mat(annot_tiles(:,:,1)))) / mat_total;
short_annot      = full(sum(mat(annot_tiles(:,:,2)))) / mat_total;
long_annot       = full(sum(mat(annot_tiles(:,:,3)))) / mat_total;
inter_annot      = full(sum(mat(annot_tiles(:,:,4)))) / mat_total;

desired_ratios = [diag_short_ratio, short_ratio, long_ratio, inter_ratio];
current_ratios = [diag_short_annot, short_annot, long_annot, inter_annot];

normalized_mat = zeros(size(mat));

for a = 1:4
    mask = annot_tiles(:,:,a);

    if desired_ratios(a) == 0
        % No observed events in this category, so assign zero mass.
        normalized_mat(mask) = 0;

    elseif current_ratios(a) > 0 && isfinite(current_ratios(a))
        normalized_mat(mask) = ...
            (desired_ratios(a) / current_ratios(a)) .* mat(mask);

    else
        % Observed events exist, but model assigns no mass to this category.
        error('renormalize_tiles:UnsupportedCategory', ...
            ['Category %d has desired ratio %g, but the model ratio is %g. ' ...
             'The model cannot represent observed events in this category.'], ...
            a, desired_ratios(a), current_ratios(a));
    end
end

normalized_total = sum(normalized_mat(:));

if normalized_total <= 0 || ~isfinite(normalized_total)
    error('renormalize_tiles:InvalidNormalizedTotal', ...
        'Normalized matrix has invalid total mass: %g.', normalized_total);
end

new_diag_short = sum(normalized_mat(annot_tiles(:,:,1))) / normalized_total;
new_short      = sum(normalized_mat(annot_tiles(:,:,2))) / normalized_total;
new_long       = sum(normalized_mat(annot_tiles(:,:,3))) / normalized_total;
new_inter      = sum(normalized_mat(annot_tiles(:,:,4))) / normalized_total;

fprintf('target diag-short, short, long, inter ratios = %g, %g, %g, %g\n', ...
    diag_short_ratio, short_ratio, long_ratio, inter_ratio);

fprintf('old diag-short, short, long, inter ratios = %g, %g, %g, %g\n', ...
    diag_short_annot, short_annot, long_annot, inter_annot);

fprintf('new diag-short, short, long, inter ratios = %g, %g, %g, %g\n', ...
    new_diag_short, new_short, new_long, new_inter);

end