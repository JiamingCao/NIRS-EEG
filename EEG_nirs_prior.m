addpath(genpath('~/Documents/MATLAB/fieldtrip-20180426'));
addpath(genpath('/home/jiaming/Documents/MATLAB/NIRFAST-9.0'));
addpath('topoplot')     % topoplot functions from EEGLAB
addpath('inhull/')
% addpath('/home/jiaming/Documents/CMUDrive/Research/Codes/Multi/FilterM')
% addpath('brainsuite_matlab')    % to read .dfs files
% Load the pre-generated models (icbm152)
clear;
% load icbm152_09c/eegmodel.mat
load icbm152_09c/eegmodel64.mat
load icbm152_09c/nirsmodel.mat
load icbm152_09c/shared.mat

clear J_*  jacobian_full
%% Generate some simple forward data
% Add some activation
% center = [146, 110, 140; 132, 110, 150];      % one under source
% center = [130, 123, 140; 119, 110, 150];      % neither under elec, both diff nirs (explodes)
% center = [146, 110, 140; 154, 126, 127];      % same elec same nirs (1)
% center = [122, 160, 142; 137, 160, 131];      % same elec no nirs (2)
center = [148, 127, 140; 146, 105, 140];      % same elec, diff nirs (3)
% center = [146, 128, 134; 154, 111, 131];      % same elec, one nirs (4)
% center = [136, 111, 147; 136, 122, 147];      % no elec same nirs (5)
% center = [82, 91, 160; 121, 96, 154];         % no elec no nirs (6)
% center = [146, 100, 140; 132, 120, 150];      % neither under elec, both diff nirs (7)
% center = [113, 79, 155; 129, 110, 153];       % no elec one nirs (8)
% center = [161, 116, 120; 131, 163, 132];      % diff elec no nirs (10)
% center = [146, 110, 140; 125, 86, 152];       % both elec both nirs (11); second eeg not showing up as supposed
    
% center = [146, 110, 140; 139, 145, 138];      % both elec one nirs (12)
% center = [131, 163, 132; 145, 160, 121];      % one elec no nirs (13)
% center = [146, 110, 140; 132, 120, 150];      % one elec,both diff nirs (14)
% center = [146, 110, 140; 162, 105, 125];      % one elec, also nirs (15)
% center = [161, 116, 120; 133, 110, 150];      % one elec, other nirs (16)

active_idx = cell(size(center, 1), 1);
for i = 1:size(center, 1)
    tmp = sort(sum((brainmesh.pos - center(i,:)).^2, 2));
    active_idx{i} = find(sum((brainmesh.pos - center(i,:)).^2, 2) < tmp(11));
%     [~, active_idx(i)] = min(sum((brainmesh.pos - center(i, :)).^2, 2));
end
n_dipole = length(idx_dipole);
Fs = 0.2;   % kHz
T_vec = 0:1/Fs:30000;    % 60 sec
n_T = length(T_vec);

%% Design stim
onset = 100:200:20e3;  % ms
peak_shape = normpdf(linspace(-3,3,100*Fs)); % 50ms peaks
peak_shape = 2*peak_shape/max(peak_shape);
onset_idx = round(onset*Fs);
stim_mat = zeros(n_dipole, n_T);
stim_mat(active_idx{1}, onset_idx) = 1;
stim_mat(active_idx{2}, onset_idx+round(50*Fs)) = 1;

%% Forward
activity_mat = 0.1 * randn(size(stim_mat)) + (filter(peak_shape, 1, stim_mat'))';

% time1 = 100;    %ms, first stim
% time2 = 125;
% time3 = 300;
% time4 = 325;
% activity_mat(active_idx{1}, find(T_vec>=time1, 1):find(T_vec>=time1, 1)+length(peak_shape)-1) = repmat(peak_shape, 10, 1);
% activity_mat(active_idx{2}, find(T_vec>=time2, 1):find(T_vec>=time2, 1)+length(peak_shape)-1) = repmat(peak_shape, 10, 1);
% activity_mat(active_idx{1}, find(T_vec>=time3, 1):find(T_vec>=time3, 1)+length(peak_shape)-1) = repmat(peak_shape, 10, 1);
% activity_mat(active_idx{2}, find(T_vec>=time4, 1):find(T_vec>=time4, 1)+length(peak_shape)-1) = repmat(peak_shape, 10, 1);

hrf = hrf0(20, Fs*1000, 0, 1);
hemodynamics = (filter(hrf, 1, stim_mat'))';
hbo_source = 0.1*randn(size(stim_mat)) + 30*hemodynamics;
hbd_source = 0.1*randn(size(stim_mat)) - 10*hemodynamics;

% hbo_source = 3*filter(hrf, 1, (activity_mat').^2);
% hbo_source = hbo_source';
% hbd_source = -hbo_source/3;
% hbo_source = 0.1 * randn(n_dipole,n_T);
% hbd_source = 0.1 * randn(n_dipole,n_T);
% hbo_source(active_idx{1}, find(T_vec>=time1, 1):find(T_vec>=time1, 1)+length(hrf)-1) = hbo_source(active_idx{1}, find(T_vec>=time1, 1):find(T_vec>=time1, 1)+length(hrf)-1)+30*repmat(hrf, 10, 1);
% hbo_source(active_idx{2}, find(T_vec>=time2, 1):find(T_vec>=time2, 1)+length(hrf)-1) = hbo_source(active_idx{2}, find(T_vec>=time2, 1):find(T_vec>=time2, 1)+length(hrf)-1)+30*repmat(hrf, 10, 1);
% hbo_source(active_idx{1}, find(T_vec>=time3, 1):find(T_vec>=time3, 1)+length(hrf)-1) = hbo_source(active_idx{1}, find(T_vec>=time3, 1):find(T_vec>=time3, 1)+length(hrf)-1)+30*repmat(hrf, 10, 1);
% hbo_source(active_idx{2}, find(T_vec>=time4, 1):find(T_vec>=time4, 1)+length(hrf)-1) = hbo_source(active_idx{2}, find(T_vec>=time4, 1):find(T_vec>=time4, 1)+length(hrf)-1)+30*repmat(hrf, 10, 1);
% hbd_source(active_idx{1}, find(T_vec>=time1, 1):find(T_vec>=time1, 1)+length(hrf)-1) = hbd_source(active_idx{1}, find(T_vec>=time1, 1):find(T_vec>=time1, 1)+length(hrf)-1)-10*repmat(hrf, 10, 1);
% hbd_source(active_idx{2}, find(T_vec>=time2, 1):find(T_vec>=time2, 1)+length(hrf)-1) = hbd_source(active_idx{2}, find(T_vec>=time2, 1):find(T_vec>=time2, 1)+length(hrf)-1)-10*repmat(hrf, 10, 1);
% hbd_source(active_idx{1}, find(T_vec>=time3, 1):find(T_vec>=time3, 1)+length(hrf)-1) = hbd_source(active_idx{1}, find(T_vec>=time3, 1):find(T_vec>=time3, 1)+length(hrf)-1)-10*repmat(hrf, 10, 1);
% hbd_source(active_idx{2}, find(T_vec>=time4, 1):find(T_vec>=time4, 1)+length(hrf)-1) = hbd_source(active_idx{2}, find(T_vec>=time4, 1):find(T_vec>=time4, 1)+length(hrf)-1)-10*repmat(hrf, 10, 1);

% take a look
% colors = nan(length(brainmesh.pos), 3);
% colors(idx_dipole, :) = vals2colormap(activity, 'jet', [-max(abs(activity)), max(abs(activity))]);
% maxeeg = find(T_vec>=time4+100,1);
% figure;
% for tstep = 1:find(T_vec>=time4+100,1)
%     cla;
%     colors =  vals2colormap(activity_mat(:,tstep), 'jet', [-max(abs(activity_mat(:))), max(abs(activity_mat(:)))]);
%     ft_plot_mesh(brainmesh, 'vertexcolor', colors);
%     colorbar, colormap('jet'), caxis([-max(abs(activity_mat(:))), max(abs(activity_mat(:)))]);
%     title(num2str(T_vec(tstep)));
%     hold on
%     for i = 1:length(nirs_mesh.link)
%         line([nirs_mesh.source.coord(nirs_mesh.link(i, 1), 1), nirs_mesh.meas.coord(nirs_mesh.link(i, 2), 1)], ...
%         [nirs_mesh.source.coord(nirs_mesh.link(i, 1), 2), nirs_mesh.meas.coord(nirs_mesh.link(i, 2), 2)], ...
%         [nirs_mesh.source.coord(nirs_mesh.link(i, 1), 3), nirs_mesh.meas.coord(nirs_mesh.link(i, 2), 3)], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 4);
%     end
%     scatter3(nirs_mesh.source.coord(:,1), nirs_mesh.source.coord(:,2), nirs_mesh.source.coord(:,3), 100, 'ro', 'filled');
%     scatter3(nirs_mesh.meas.coord(:,1), nirs_mesh.meas.coord(:,2), nirs_mesh.meas.coord(:,3), 100, 'bo', 'filled');
%     pause(0.2);
% end

% Scalp eeg recordings
eeg = L * activity_mat*100;
eeg = eeg + sqrt(0.01 * (eeg.^2)) .* randn(size(eeg));
elec = [];
elec.pnt = headmodel.elec.chanpos;
elec.pnt = elec.pnt - mean(elec.pnt);
elec.label = headmodel.elec.label;
% figure;
% for tstep = 1:find(T_vec>=time4+100,1)
%     cla;
%     topoplot(eeg(:,tstep),elec,'colormap','jet','nosedir','+Y', 'electrodes','numbers');
%     caxis([-0.06, 0.06]);
%     title(num2str(T_vec(tstep)));
%     pause(0.2);
% end

% NIRS recordings
n_nirs = size(jacobian_nirs,1)/2;
dOD = jacobian_nirs * [hbo_source; hbd_source];
dOD = dOD + [randn(n_nirs,n_T).*sqrt(0.1*(dOD(1:n_nirs,:).^2));randn(n_nirs,n_T).*sqrt(0.05*(dOD(n_nirs+1:end,:).^2))];

% plot_nirs_data(dOD, nirs_mesh, {'750', '850'})

% Beer-Lambert
hemoglobin = zeros(2*n_nirs, n_T);
src = nirs_mesh.source.coord;
det = nirs_mesh.meas.coord;
link = nirs_mesh.link;
DPF750 = 6.59;
DPF850 = 5.82;  % General equation for the differential pathlength factor of the frontal human head depending on wavelength and age; Scholkmann & Wolf
excoeff = importdata('excoef.txt');
excoeff = excoeff.data;
epsilon = [excoeff(find(excoeff(:,1)==750), 2:3); excoeff(find(excoeff(:,1)==850), 2:3)];

for tstep = 1:n_T
    dOD2 = [dOD(1:n_nirs, tstep), dOD(n_nirs+1:end, tstep)];
    tmp = zeros(n_nirs,2);
    for ch = 1:n_nirs
        dist = norm(src(link(ch, 1),:) - det(link(ch, 2),:));
        A = epsilon .* [DPF750, DPF850; DPF750, DPF850] * dist;
        tmp(ch, :) = -pinv(A)*dOD2(ch, :)';
    end
    hemoglobin(:, tstep) = tmp(:);
end

% plot_nirs_data(hemoglobin(:, 26), nirs_mesh, {'HbO', 'HbD'});
    

%% Inverse
% EEG only
Qn_eeg = {speye(size(eeg, 1))};
Qp_eeg = {speye(size(L, 2))};
% tmp = sparse(diag(0.1*rand(n_dipole, 1)));
% tmp(active_idx, active_idx) = 1;
% Qp_eeg = {tmp};
% eeg_trial = (eeg(:,T_vec>time1-50&T_vec<time2+100)+eeg(:,T_vec>time3-50&T_vec<time4+100))/2;
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
% figure;
% for tstep = 1:maxeeg
%     cla;
%     colors = vals2colormap(Beta_eeg(:,tstep), 'jet', [-20, 20]);
%     ft_plot_mesh(brainmesh, 'vertexcolor', colors);
%     colormap('jet'), caxis([-20, 20]);
%     title(['eeg only ',num2str(T_vec(tstep)), 'ms']);
%     pause(0.2)
% end

%% average NIRS
Qn_nirs = {blkdiag(speye(n_nirs), sparse(n_nirs, n_nirs)), blkdiag(sparse(n_nirs, n_nirs), speye(n_nirs))};
% Qp_nirs = {blkdiag(speye(size(jacobian_nirs, 2)/2), sparse(size(jacobian_nirs, 2)/2, size(jacobian_nirs, 2)/2)), blkdiag(sparse(size(jacobian_nirs, 2)/2, size(jacobian_nirs, 2)/2), speye(size(jacobian_nirs, 2)/2))};
Qp_nirs = {blkdiag(speye(n_dipole), sparse(n_dipole, n_dipole)), blkdiag(sparse(n_dipole, n_dipole), speye(n_dipole)),...
    [sparse(n_dipole, n_dipole), -speye(n_dipole);-speye(n_dipole), sparse(n_dipole, n_dipole)]};
dOD_avg = mean(dOD,2);
[~, Beta_nirs, ~] = REML(dOD_avg, jacobian_nirs, [],Qn_nirs, Qp_nirs, 500);
% colors = vals2colormap(Beta_nirs(1:length(Beta_nirs)/2), 'jet', [-max(abs(Beta_nirs(1:length(Beta_nirs)/2))), max(abs(Beta_nirs(1:length(Beta_nirs)/2)))]);
% figure, ft_plot_mesh(brainmesh, 'vertexcolor', colors);
% % colorbar;
% colormap('jet'), caxis([-max(abs(Beta_nirs(1:length(Beta_nirs)/2))), max(abs(Beta_nirs(1:length(Beta_nirs)/2)))]);
% title('HbO')
% 
% colors = vals2colormap(Beta_nirs(length(Beta_nirs)/2+1:end), 'jet', [-max(abs(Beta_nirs(length(Beta_nirs)/2+1:end))), max(abs(Beta_nirs(length(Beta_nirs)/2+1:end)))]);
% figure, ft_plot_mesh(brainmesh, 'vertexcolor', colors);
% % colorbar;
% colormap('jet'), caxis([-max(abs(Beta_nirs(length(Beta_nirs)/2+1:end))), max(abs(Beta_nirs(length(Beta_nirs)/2+1:end)))]);
% title('HbD')
% 
% dOD2 = [dOD_avg(1:length(dOD_avg)/2), dOD_avg(length(dOD_avg)/2+1:end)];
% for ch = 1:length(link)
%     dist = norm(src(link(ch, 1),:) - det(link(ch, 2),:));
%     A = epsilon .* [DPF750, DPF850; DPF750, DPF850] * dist;
%     hemoglobin2(ch, :) = -pinv(A)*dOD2(ch, :)';
% end
% plot_nirs_data(hemoglobin2(:), nirs_mesh, {'HbO', 'HbD'});
%% NIRS projection and prior using NIRS projection
Beta_HbO_proj = Proj_to_cortex(nirs_mesh, hemoglobin, brainmesh);
% take a look
colors =  vals2colormap(Beta_HbO_proj, 'jet', [-1, 1]);
figure, ft_plot_mesh(brainmesh, 'vertexcolor', colors);
colormap('jet'), caxis([-1, 1]);
title('HbO Projection')

Qp_eeg_proj = {spdiags(1-exp(-(Beta_HbO_proj + 0.25)/1), 0, size(L, 2), size(L, 2))};
Beta_eeg_proj = zeros(size(L, 2), maxeeg);
parfor tstep = 1:maxeeg
    [~, Beta_eeg_proj(:,tstep), ~] = REML(eeg_trial(:, tstep), L, [],Qn_eeg, Qp_eeg_proj, 500);
end
% figure;
% for tstep = 1:maxeeg
%     cla;
% %     colors = vals2colormap(Beta_eeg_proj(:,tstep), 'jet', [-20, 20]);
%     ft_plot_mesh(brainmesh, 'vertexcolor', Beta_eeg_proj(:,tstep), 'colormap', brain_cmap, 'clim', [-20,20]);
%     camlight headlight
% %     colormap('jet'), caxis([-20, 20]);
%     title(['eeg nirs\_proj ',num2str(T_vec(tstep)), 'ms']);
%     pause(0.2)
% end

%% EEG with NIRS prior
% Qn_eeg2 = {speye(length(eeg))};
nirspow=abs(Beta_nirs(1:length(Beta_nirs)/2));
thresh = 0.1;
nirspow_norm = nirspow.*(nirspow > thresh) / max(nirspow);
% thresh = 0.25;
% tmp = speye(n_dipole);
% tmp = tmp+spdiags(nirspow>thresh,0,n_dipole,n_dipole)*9;
% Qp_eeg2 = {tmp};
% Qp_eeg2 = {speye(size(L, 2)), spdiags(1-exp(-(nirspow_norm+0.2)/2), 0, size(L, 2), size(L, 2))};
Qp_eeg2 = {spdiags(1-exp(-(nirspow_norm + 0.25)/1), 0, size(L, 2), size(L, 2))};
Beta_eeg2 = zeros(size(L, 2), maxeeg);
parfor tstep = 1:maxeeg
    [~, Beta_eeg2(:,tstep), ~] = REML(eeg_trial(:, tstep), L, [],Qn_eeg, Qp_eeg2, 500);
end
% figure;
% for tstep = 1:maxeeg
%     cla;
% %     colors = vals2colormap(Beta_eeg2(:,tstep), 'jet', [-20, 20]);
%     ft_plot_mesh(brainmesh, 'vertexcolor', Beta_eeg2(:,tstep), 'colormap', brain_cmap, 'clim', [-20,20]);
%     camlight headlight
% %     colormap('jet'), caxis([-20, 20]);
%     title(['eeg nirs ',num2str(T_vec(tstep)), 'ms']);
%     pause(0.2)
% end
%% some plotting
condition = 'no_elec_diff_nirs';
if ~exist(['Fig_JournalPaper/', condition], 'dir')
    mkdir(['Fig_JournalPaper/', condition]);
end
save_dir = ['Fig_JournalPaper/', condition, '/'];
% compare of eeg
figure;
subplot(1,2,1),subplot(1,2,2);
for tstep = 1:maxeeg
    subplot(1,2,1),cla;
%     colors = vals2colormap(Beta_eeg(:,tstep), 'jet', [-20, 20]);
    ft_plot_mesh(brainmesh, 'vertexcolor', Beta_eeg(:,tstep), 'colormap', redblue, 'clim', [-20,20]);
    camlight headlight
%     colormap('jet'), caxis([-20, 20]);
    title(['EEG only ',num2str(T_vec(tstep)), 'ms']);
    subplot(1,2,2),cla;
%     colors = vals2colormap(Beta_eeg2(:,tstep), 'jet', [-20, 20]);
    ft_plot_mesh(brainmesh, 'vertexcolor', Beta_eeg2(:,tstep), 'colormap', redblue, 'clim', [-20,20]);
    camlight headlight
%     colormap('jet'), caxis([-20, 20]);
%     colorbar;
    title(['With NIRS prior ',num2str(T_vec(tstep)), 'ms']);
    
    frame = getframe(gcf);
    img = frame2im(frame);
    [img, cmap] = rgb2ind(img,256);
    if tstep==1
        imwrite(img, cmap, [save_dir,'compare_',condition, '64.gif'],'gif','LoopCount', Inf,'Delaytime', 0.2);
    else
        imwrite(img, cmap, [save_dir,'compare_',condition,'64.gif'],'gif','WriteMode', 'append','Delaytime', 0.2);
    end
%     pause(0.2)
end

idx_t(1) = find(T_vec == 50);
idx_t(2) = find(T_vec == 75);
idx_t(3) = find(T_vec == 100);
figure;
for i=1:3
    cla;
%     subplot(1,3,i);
%     colors = vals2colormap(Beta_eeg2(:,idx_t(i)), 'jet', [-20, 20]);
    ft_plot_mesh(brainmesh, 'vertexcolor', Beta_eeg2(:,idx_t(i)), 'colormap', redblue, 'clim', [-max(Beta_eeg2(:,idx_t(i))),max(Beta_eeg2(:,idx_t(i)))]);
    camlight headlight
    xlim([96,170]),ylim([20,220]),zlim([0,160])
%     colormap('jet'), caxis([-20, 20]);
%     title([num2str(T_vec(idx_t(i))), 'ms']);
    savefig([save_dir,num2str(T_vec((idx_t(i)))),'ms_', condition,'64']);
    saveas(gcf, [save_dir,num2str(T_vec((idx_t(i)))),'ms_', condition, '64.png']);
end
% suptitle('EEG with NIRS prior');

% ground truth
true_trial = zeros(n_dipole, maxeeg);
for i=1:length(onset)
    true_trial = true_trial+activity_mat(:,onset_idx(i):onset_idx(i)+200*Fs-1);
end
true_trial = true_trial/length(onset);

figure;
for tstep = 1:maxeeg
    cla;
%     colors = vals2colormap(true_trial(:,tstep), 'jet', [-2, 2]);
    ft_plot_mesh(brainmesh, 'vertexcolor', true_trial(:,tstep), 'colormap', redblue, 'clim', [-2,2]);
    camlight headlight
    xlim([96,170]),ylim([20,220]),zlim([0,160])
%     colormap('jet'), caxis([-2, 2]);
    title(['Ground truth ',num2str(T_vec(tstep)), 'ms']);
    
    frame = getframe(gcf);
    img = frame2im(frame);
    [img, cmap] = rgb2ind(img,256);
    if tstep==1
        imwrite(img, cmap, [save_dir,'truth_',condition,'64.gif'],'gif','LoopCount', Inf,'Delaytime', 0.2);
    else
        imwrite(img, cmap,[save_dir,'truth_',condition,'64.gif'],'gif','WriteMode', 'append','Delaytime', 0.2);
    end
%     pause(0.2)
end

figure;
for i=1:3
    cla;
%     subplot(1,3,i);
%     colors = vals2colormap(true_trial(:,idx_t(i)), 'jet', [-2, 2]);
    ft_plot_mesh(brainmesh, 'vertexcolor', true_trial(:,idx_t(i)), 'colormap', redblue, 'clim', [-2,2]);
    camlight headlight
    xlim([96,170]),ylim([20,220]),zlim([0,160])
%     title([num2str(T_vec(idx_t(i))), 'ms']);
    savefig([save_dir,num2str(T_vec((idx_t(i)))),'ms_truth_', condition,'64']);
    saveas(gcf, [save_dir,num2str(T_vec((idx_t(i)))),'ms_truth_', condition, '64.png']);
end
% suptitle('Ground truth');

figure;
% subplot(1,2,1);
% colors = vals2colormap(Beta_nirs(1:length(Beta_nirs)/2), 'jet', [-max(abs(Beta_nirs(1:length(Beta_nirs)/2))), max(abs(Beta_nirs(1:length(Beta_nirs)/2)))]);
ft_plot_mesh(brainmesh, 'vertexcolor', Beta_nirs(1:length(Beta_nirs)/2), 'colormap', redblue, 'clim', [-max(abs(Beta_nirs(1:length(Beta_nirs)/2))), max(abs(Beta_nirs(1:length(Beta_nirs)/2)))]);
camlight headlight
% colormap('jet'), caxis([-max(abs(Beta_nirs(1:length(Beta_nirs)/2))), max(abs(Beta_nirs(1:length(Beta_nirs)/2)))]);
% title('HbO')

% subplot(1,2,2);
% colors = vals2colormap(Beta_nirs(length(Beta_nirs)/2+1:end), 'jet', [-max(abs(Beta_nirs(length(Beta_nirs)/2+1:end))), max(abs(Beta_nirs(length(Beta_nirs)/2+1:end)))]);
% ft_plot_mesh(brainmesh, 'vertexcolor', colors);
% colormap('jet'), caxis([-max(abs(Beta_nirs(length(Beta_nirs)/2+1:end))), max(abs(Beta_nirs(length(Beta_nirs)/2+1:end)))]);
% title('HbD')

savefig([save_dir,'nirs_recon_',condition,'64']);
saveas(gcf, [save_dir,'nirs_recon_',condition, '64.png']);
%% EEG with NIRS projection
% compare of eeg
figure;
subplot(1,2,1),subplot(1,2,2);
for tstep = 1:maxeeg
    subplot(1,2,1),cla;
%     colors = vals2colormap(Beta_eeg(:,tstep), 'jet', [-20, 20]);
    ft_plot_mesh(brainmesh, 'vertexcolor', Beta_eeg(:,tstep), 'colormap', redblue, 'clim', [-20,20]);
    camlight headlight
%     colormap('jet'), caxis([-20, 20]);
    title(['EEG only ',num2str(T_vec(tstep)), 'ms']);
    subplot(1,2,2),cla;
%     colors = vals2colormap(Beta_eeg_proj(:,tstep), 'jet', [-20, 20]);
    ft_plot_mesh(brainmesh, 'vertexcolor', Beta_eeg_proj(:,tstep), 'colormap', redblue, 'clim', [-20,20]);
    camlight headlight
%     colormap('jet'), caxis([-20, 20]);
    title(['With NIRS projection prior ',num2str(T_vec(tstep)), 'ms']);
    
    frame = getframe(gcf);
    img = frame2im(frame);
    [img, cmap] = rgb2ind(img,256);
    if tstep==1
        imwrite(img, cmap, [save_dir,'compare_',condition,'_nirsproj64.gif'],'gif','LoopCount', Inf,'Delaytime', 0.2);
    else
        imwrite(img, cmap, [save_dir,'compare_',condition,'_nirsproj64.gif'],'gif','WriteMode', 'append','Delaytime', 0.2);
    end
end
 
idx_t(1) = find(T_vec == 50);
idx_t(2) = find(T_vec == 75);
idx_t(3) = find(T_vec == 100);
figure;
for i=1:3
    cla;
%     subplot(1,3,i);
%     colors = vals2colormap(Beta_eeg2(:,idx_t(i)), 'jet', [-20, 20]);
    ft_plot_mesh(brainmesh, 'vertexcolor', Beta_eeg_proj(:,idx_t(i)), 'colormap', redblue, 'clim', [-max(Beta_eeg_proj(:,idx_t(i))),max(Beta_eeg_proj(:,idx_t(i)))]);
    camlight headlight
    xlim([96,170]),ylim([20,220]),zlim([0,160])
%     colormap('jet'), caxis([-20, 20]);
%     title([num2str(T_vec(idx_t(i))), 'ms']);
    savefig([save_dir,num2str(T_vec((idx_t(i)))),'ms_', condition, '_nirsproj64']);
    saveas(gcf, [save_dir,num2str(T_vec((idx_t(i)))),'ms_', condition, '_nirsproj64.png']);
end

% colors =  vals2colormap(Beta_HbO_proj, 'jet', [-1, 1]);
figure, ft_plot_mesh(brainmesh, 'vertexcolor', Beta_HbO_proj, 'colormap', redblue, 'clim', [-1,1]);
camlight headlight
% colormap('jet'), caxis([-1, 1]);
% title('HbO Projection')
savefig([save_dir,'nirs_proj_',condition,'64']);
saveas(gcf, [save_dir,'nirs_proj_',condition, '64.png']);

%% eeg only
figure;
for i=1:3
    cla;
%     subplot(1,3,i);
%     colors = vals2colormap(Beta_eeg(:,idx_t(i)), 'jet', [-20, 20]);
    ft_plot_mesh(brainmesh, 'vertexcolor', Beta_eeg(:,idx_t(i)), 'colormap', redblue, 'clim', [-max(Beta_eeg(:,idx_t(i))),max(Beta_eeg(:,idx_t(i)))]);
    camlight headlight
    xlim([96,170]),ylim([20,220]),zlim([0,160])
%     colormap('jet'), caxis([-20, 20]);
%     title([num2str(T_vec(idx_t(i))), 'ms']);
    savefig([save_dir,num2str(T_vec((idx_t(i)))),'ms_eeg_', condition,'64']);
    saveas(gcf, [save_dir,num2str(T_vec((idx_t(i)))),'ms_eeg_', condition, '64.png']);
end
% suptitle('EEG only');

%% sensor locations
activity = zeros(length(idx_dipole), 1);
for i=1:length(active_idx)
    activity(active_idx{i}) = 1;
end
% colors =  vals2colormap(activity, 'jet', [-max(abs(activity)), max(abs(activity))]);
figure, ft_plot_mesh(brainmesh, 'vertexcolor', activity, 'colormap', redblue, 'clim', [-max(abs(activity)), max(abs(activity))]);
camlight headlight

% colormap('jet'), caxis([-max(abs(activity)), max(abs(activity))]);
% Plot the optodes and electrodes
hold on
for i = 1:length(nirs_mesh.link)
    line([nirs_mesh.source.coord(nirs_mesh.link(i, 1), 1), nirs_mesh.meas.coord(nirs_mesh.link(i, 2), 1)], ...
    [nirs_mesh.source.coord(nirs_mesh.link(i, 1), 2), nirs_mesh.meas.coord(nirs_mesh.link(i, 2), 2)], ...
    [nirs_mesh.source.coord(nirs_mesh.link(i, 1), 3), nirs_mesh.meas.coord(nirs_mesh.link(i, 2), 3)], 'Color', [0,0.5,0], 'LineStyle', '--', 'LineWidth', 3);
end
scatter3(nirs_mesh.source.coord(:,1), nirs_mesh.source.coord(:,2), nirs_mesh.source.coord(:,3), 100, 'ro','LineWidth', 2, 'MarkerEdgeAlpha', 0.6);
scatter3(nirs_mesh.meas.coord(:,1), nirs_mesh.meas.coord(:,2), nirs_mesh.meas.coord(:,3), 100, 'bo','LineWidth', 2, 'MarkerEdgeAlpha', 0.6);
scatter3(headmodel.elec.chanpos(:,1), headmodel.elec.chanpos(:,2), headmodel.elec.chanpos(:,3), 100, 'ko', 'filled');
savefig([save_dir,'location_probe_',condition,'64']);
saveas(gcf, [save_dir,'location_probe_',condition, '64.png']);

