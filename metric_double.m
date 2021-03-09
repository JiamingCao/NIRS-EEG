function [metricA, metricB] = metric_double(center, data, brainmesh, thresh)
data(brainmesh.pos(:,3)<120)=0; % consider only the dorsal side of the brain, to avoid artefacts from spikes from unrelated regions
data = data/max(data);  % normalize
useful_data = data>=thresh;
centerA = center(1, :);
centerB = center(2, :);
distance_to_A2 = (brainmesh.pos(:,1) - centerA(1)).^2 + (brainmesh.pos(:,2) - centerA(2)).^2 + (brainmesh.pos(:,3) - centerA(3)).^2;
distance_to_B2 = (brainmesh.pos(:,1) - centerB(1)).^2 + (brainmesh.pos(:,2) - centerB(2)).^2 + (brainmesh.pos(:,3) - centerB(3)).^2;
A_recon = useful_data & (distance_to_A2 < distance_to_B2);
B_recon = useful_data & (distance_to_A2 >= distance_to_B2);
if any(A_recon) && any(B_recon)
    distances = pdist2(brainmesh.pos(A_recon,:), brainmesh.pos(B_recon,:));
else
    distances = 0;
end

if ~any(A_recon) && any(B_recon)
    metricA = nan;
    centerB_recon = mean(brainmesh.pos(useful_data, :), 1);
    bias = norm(centerB_recon - centerB);
    variance = mean((brainmesh.pos(useful_data,1) - centerB_recon(1)).^2 + (brainmesh.pos(useful_data,2) - centerB_recon(2)).^2 + (brainmesh.pos(useful_data,3) - centerB_recon(3)).^2);
    metricB = sqrt(bias^2 + variance);
elseif ~any(B_recon) && any(A_recon)
    centerA_recon = mean(brainmesh.pos(useful_data, :), 1);
    bias = norm(centerA_recon - centerA);
    variance = mean((brainmesh.pos(useful_data,1) - centerA_recon(1)).^2 + (brainmesh.pos(useful_data,2) - centerA_recon(2)).^2 + (brainmesh.pos(useful_data,3) - centerA_recon(3)).^2);
    metricA = sqrt(bias^2 + variance);
    metricB = nan;
elseif min(distances(:))<5
    center_recon = mean(brainmesh.pos(useful_data, :), 1);
    biasA = norm(center_recon - centerA);
    variance = mean((brainmesh.pos(useful_data,1) - center_recon(1)).^2 + (brainmesh.pos(useful_data,2) - center_recon(2)).^2 + (brainmesh.pos(useful_data,3) - center_recon(3)).^2);
    metricA = sqrt(biasA^2 + variance);
    biasB = norm(center_recon - centerB);
    variance = mean((brainmesh.pos(useful_data,1) - center_recon(1)).^2 + (brainmesh.pos(useful_data,2) - center_recon(2)).^2 + (brainmesh.pos(useful_data,3) - center_recon(3)).^2);
    metricB = sqrt(biasB^2 + variance);
else
    bias = norm(centerA - mean(brainmesh.pos(A_recon, :), 1));
    variance = mean(distance_to_A2(A_recon));
    metricA = sqrt(bias^2 + variance);
    bias = norm(centerB - mean(brainmesh.pos(B_recon, :), 1));
    variance = mean(distance_to_B2(B_recon));
    metricB = sqrt(bias^2 + variance);
end
