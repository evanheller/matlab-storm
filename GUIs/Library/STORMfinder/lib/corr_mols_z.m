function [matched,unmatched] = corr_mols_z(set1_pos,set2_pos,tform_start,match_radius,max_zdiff)

% This function takes two input sets of positions and matches the points to
% each other based on the initial transform tform_start and the radial
% tolerance match_radius. set1 is the master set, and set2 is matched to
% it. Both input variables set1_pos and set2_pos must be structs with
% length 1. set1_pos and set2_pos should have fields .x and .y that contain
% vectors of the x and y positions. The input tform_start is the spatial
% transformation matrix (see cp2tform). The output matched and unmatched
% each have fields set1_inds and set2_inds.
%
% SYNTAX: [matched,unmatched] = corr_mols(set1_pos,set2_pos,tform_start,match_radius)
% By Josh Vaugan 

% with additions by Alistair Boettiger

% initialize variables
num_mols_set1 = length(set1_pos.x);
matching_set2_inds = zeros(1,num_mols_set1);


if isa(set1_pos.x,'single') 
    set1_pos.x = double(set1_pos.x);
    set1_pos.y = double(set1_pos.y);
    set1_pos.z = double(set1_pos.z);
    set2_pos.x = double(set2_pos.x);
    set2_pos.y = double(set2_pos.y);
    set2_pos.z = double(set2_pos.z);
end

% overlap set2 onto set1 (shift, or also warp, depending on tform_start)
[registered_set2_pos.x,registered_set2_pos.y] = tforminv(tform_start,set2_pos.x',set2_pos.y');
registered_set2_pos.z=double(set2_pos.z');

% correlate molecules in the two sets

for n = 1:num_mols_set1 % loop over each set1 index
    xdiff = set1_pos.x(n) - registered_set2_pos.x;
    ydiff = set1_pos.y(n) - registered_set2_pos.y;
    zdiff = set1_pos.z(n) - registered_set2_pos.z;
    matching_set2_ind = find(sqrt(xdiff.^2 + ydiff.^2)<match_radius & abs(zdiff) < max_zdiff );

    if length(matching_set2_ind)==1 
       matching_set2_inds(n) = matching_set2_ind;
    end
end

% get only the matching molecules
matched.set1_inds = find(matching_set2_inds>0);
matched.set2_inds = matching_set2_inds(matched.set1_inds);
unmatched.set1_inds = find(matching_set2_inds==0);
unmatched.set2_inds = setdiff(1:length(set2_pos.x),matched.set2_inds);

