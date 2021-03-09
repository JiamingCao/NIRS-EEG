[A,B] = meshgrid(0.04:0.02:0.2, 0.4:0.2:5);
params = [A(:), B(:)];
idx_t(1) = find(T_vec == 50);
idx_t(2) = find(T_vec == 75);
idx_t(3) = find(T_vec == 100);
metric_thresh = 0.5; 
ratios = zeros(4, length(params));
parfor i=1:length(params)
    fprintf('Trial %d\n', i);
    metrics = zeros(4, 2);
    a = params(i,1);
    b = params(i,2);
    Qp_eeg2 = {spdiags(1-exp(-(nirspow_norm + a)/b), 0, size(L, 2), size(L, 2))};
    Beta_eeg2 = zeros(size(L, 2), maxeeg);
    for tstep = 1:maxeeg
        [~, Beta_eeg2(:,tstep), ~] = REML(eeg_trial(:, tstep), L, [],Qn_eeg, Qp_eeg2, 500);
    end
    metrics(1,1) = metric_single(center(1,:), Beta_eeg(:, idx_t(1)), brainmesh, metric_thresh);
    metrics(1,2) = metric_single(center(1,:), Beta_eeg2(:, idx_t(1)), brainmesh, metric_thresh);
    [metrics(2,1), metrics(3,1)] = metric_double(center, Beta_eeg(:, idx_t(2)), brainmesh, metric_thresh);
    [metrics(2,2), metrics(3,2)] = metric_double(center, Beta_eeg2(:, idx_t(2)), brainmesh, metric_thresh);
    metrics(4,1) = metric_single(center(2,:), Beta_eeg(:, idx_t(3)), brainmesh, metric_thresh);
    metrics(4,2) = metric_single(center(2,:), Beta_eeg2(:, idx_t(3)), brainmesh, metric_thresh);
    ratios(:,i) = metrics(:,2)./metrics(:,1);
end