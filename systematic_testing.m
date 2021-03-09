addpath(genpath('~/Documents/MATLAB/fieldtrip-20180426'));
addpath(genpath('/home/jiaming/Documents/MATLAB/NIRFAST-9.0'));
% Load the pre-generated models (icbm152)
clear;
load icbm152_09c/eegmodel.mat
% load icbm152_09c/eegmodel64.mat
load icbm152_09c/nirsmodel.mat
% load icbm152_09c/nirsmodel_hd.mat
load icbm152_09c/shared.mat

clear J_*  jacobian_full

% starts here
xyrange = [110, 145, 90, 130];
z_thresh = 120;
N_trials = 500;

idx = brainmesh.pos(:,1)>xyrange(1) & brainmesh.pos(:,1)<xyrange(2) & brainmesh.pos(:,2)>xyrange(3) & brainmesh.pos(:,2)<xyrange(4) & brainmesh.pos(:,3)>z_thresh;
pool_loc = brainmesh.pos(idx, :);
n_loc = sum(idx);
dist_mat = pdist2(pool_loc, pool_loc);

Fs = 0.2;   % kHz
T_vec = 0:1/Fs:30000;    % 60 sec
n_T = length(T_vec);
idx_t(1) = find(T_vec == 50);
idx_t(2) = find(T_vec == 75);
idx_t(3) = find(T_vec == 100);
metric_thresh = 0.5; % for metric calculation
improved = zeros(N_trials, 1);

metrics_tot = zeros(4,2,N_trials);
centers_tot = zeros(2,3,N_trials);

for N = 1:N_trials
    fprintf('Trial %d\n', N);
    metrics = zeros(4, 2); % [50ms EEG, 50ms Joint; 75ms EEG A, 75ms Joint A; 75ms EEG B, 75mm Joint B; 100ms EEG, 100ms Joint]
    
    center_idx1 = randperm(n_loc, 1);
    center1 = pool_loc(center_idx1, :);
    tmp_idx = find(dist_mat(center_idx1, :)>20 & dist_mat(center_idx1, :)<30);
    center2 = pool_loc(tmp_idx(randperm(length(tmp_idx), 1)), :);
    center = [center1; center2];
    centers_tot(:,:,N) = center;
    
    active_idx = cell(size(center, 1), 1);
    for i = 1:size(center, 1)
        tmp = sort(sum((brainmesh.pos - center(i,:)).^2, 2));
        active_idx{i} = find(sum((brainmesh.pos - center(i,:)).^2, 2) < tmp(11));
    %     [~, active_idx(i)] = min(sum((brainmesh.pos - center(i, :)).^2, 2));
    end
    n_dipole = length(idx_dipole);
    
    % Design stim
    onset = 100:200:20e3;  % ms
    peak_shape = normpdf(linspace(-3,3,100*Fs)); % 50ms peaks
    peak_shape = 2*peak_shape/max(peak_shape);
    onset_idx = round(onset*Fs);
    stim_mat = zeros(n_dipole, n_T);
    stim_mat(active_idx{1}, onset_idx) = 1;
    stim_mat(active_idx{2}, onset_idx+round(50*Fs)) = 1;

    % Forward
    activity_mat = 0.1 * randn(size(stim_mat)) + (filter(peak_shape, 1, stim_mat'))';
    hrf = hrf0(20, Fs*1000, 0, 1);
    hemodynamics = (filter(hrf, 1, stim_mat'))';
    hbo_source = 12*randn(size(stim_mat))*max(hemodynamics(:)) + 3*hemodynamics;
    hbd_source = 2*randn(size(stim_mat))*max(hemodynamics(:)) - hemodynamics;
    eeg = L * activity_mat*100;
    eeg = eeg + sqrt(0.1 * (eeg.^2)) .* randn(size(eeg));
    elec = [];
    elec.pnt = headmodel.elec.chanpos;
    elec.pnt = elec.pnt - mean(elec.pnt);
    elec.label = headmodel.elec.label;

    % DOT recordings
    n_nirs = size(jacobian_nirs,1)/2;
    dOD = jacobian_nirs * [hbo_source; hbd_source];
    dOD = dOD + [randn(n_nirs,n_T).*sqrt(0.1*(dOD(1:n_nirs,:).^2));randn(n_nirs,n_T).*sqrt(0.05*(dOD(n_nirs+1:end,:).^2))];
    
    % Inverse
    % EEG only
    Qn_eeg = {speye(size(eeg, 1))};
    Qp_eeg = {speye(size(L, 2))};
    eeg_trial = zeros(size(L,1), 200*Fs);
    for i=1:length(onset)
        eeg_trial = eeg_trial + eeg(:, onset_idx(i):onset_idx(i)+200*Fs-1);
    end
    eeg_trial = eeg_trial/length(onset);
    maxeeg=size(eeg_trial,2);
    Beta_eeg = zeros(size(L, 2), maxeeg);
    parfor tstep = 1:maxeeg
        [~, Beta_eeg(:,tstep), ~] = REML(eeg_trial(:, tstep), L, [],Qn_eeg, Qp_eeg, 500);
    end
    
    % average NIRS
    Qn_nirs = {blkdiag(speye(n_nirs), sparse(n_nirs, n_nirs)), blkdiag(sparse(n_nirs, n_nirs), speye(n_nirs))};
    Qp_nirs = {blkdiag(speye(size(jacobian_nirs, 2)/2), sparse(size(jacobian_nirs, 2)/2, size(jacobian_nirs, 2)/2)), blkdiag(sparse(size(jacobian_nirs, 2)/2, size(jacobian_nirs, 2)/2), speye(size(jacobian_nirs, 2)/2))};
    % Qp_nirs = {blkdiag(speye(n_dipole), sparse(n_dipole, n_dipole)), blkdiag(sparse(n_dipole, n_dipole), speye(n_dipole)),...
    %     [sparse(n_dipole, n_dipole), -speye(n_dipole);-speye(n_dipole), sparse(n_dipole, n_dipole)]};

    [b,a]=butter(2, 1/(Fs*1000/2));
    dOD = filtfilt(b, a, dOD')';
    dOD_avg = mean(dOD,2);
    [~, Beta_nirs, ~] = REML(dOD_avg, jacobian_nirs, [],Qn_nirs, Qp_nirs, 500);
    % Joint reconstruction
    nirspow=abs(Beta_nirs(1:length(Beta_nirs)/2));
    thresh = 0.1;
    nirspow_norm = nirspow.*(nirspow > thresh) / max(nirspow);
    Qp_eeg2 = {spdiags(1-exp(-(nirspow_norm + 0.1)/1), 0, size(L, 2), size(L, 2))};
    Beta_eeg2 = zeros(size(L, 2), maxeeg);
    parfor tstep = 1:maxeeg
        [~, Beta_eeg2(:,tstep), ~] = REML(eeg_trial(:, tstep), L, [],Qn_eeg, Qp_eeg2, 500);
    end
    
    % metrics
    metrics(1,1) = metric_single(center(1,:), Beta_eeg(:, idx_t(1)), brainmesh, metric_thresh);
    metrics(1,2) = metric_single(center(1,:), Beta_eeg2(:, idx_t(1)), brainmesh, metric_thresh);
    [metrics(2,1), metrics(3,1)] = metric_double(center, Beta_eeg(:, idx_t(2)), brainmesh, metric_thresh);
    [metrics(2,2), metrics(3,2)] = metric_double(center, Beta_eeg2(:, idx_t(2)), brainmesh, metric_thresh);
    metrics(4,1) = metric_single(center(2,:), Beta_eeg(:, idx_t(3)), brainmesh, metric_thresh);
    metrics(4,2) = metric_single(center(2,:), Beta_eeg2(:, idx_t(3)), brainmesh, metric_thresh);
    
    metrics_tot(:,:,N) = metrics;
    
    improved(N) = all(metrics(:,1) > metrics(:,2));
end

save('systematic_smaller_roi.mat','improved','metrics_tot', 'centers_tot', 'xyrange')
