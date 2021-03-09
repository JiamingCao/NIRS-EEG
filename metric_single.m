function metric = metric_single(center, data, brainmesh, thresh)
data(brainmesh.pos(:,3)<120)=0;     % consider only the dorsal side of the brain, to avoid artefacts from spikes from unrelated regions
data = data/max(data);  % normalize
useful_data = data>=thresh;
center_recon = mean(brainmesh.pos(useful_data, :), 1);
bias = norm(center_recon - center);
variance = mean((brainmesh.pos(useful_data,1) - center_recon(1)).^2 + (brainmesh.pos(useful_data,2) - center_recon(2)).^2 + (brainmesh.pos(useful_data,3) - center_recon(3)).^2);
metric = sqrt(bias^2 + variance);
