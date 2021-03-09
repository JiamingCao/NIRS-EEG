% This is the right one
clear;
addpath(genpath('~/Documents/MATLAB/fieldtrip-20180426'));
addpath(genpath('~/Documents/MATLAB/iso2mesh-master'));
addpath(genpath('~/Documents/MATLAB/NIRFAST-9.0'));
% Use the presegmented icbm152_09c brain
% Use a simpler layered model, as in BEM
tmp = ft_read_mri('mni_icbm152_nlin_asym_09c_nifti/mni_icbm152_nlin_asym_09c/mni_icbm152_t1_tal_nlin_asym_09c.skull.label.nii.gz');
fprintf('Meshing...\n')
opt = [];   % set up the structure for iso2mesh, and create tetrahedral mesh
opt.maxnode = 1e6;
opt.radbound = 0.75;
iso_thresh = 15.5:1:18.5;
[node, elem, face] = v2m(uint8(tmp.anatomy), iso_thresh, opt, 5, 'cgalmesh');
node = node(:, 1:3);
% [pial_n, pial_f]=v2s(uint8(tmp.anatomy>18.519),0.5,2,'cgalmesh');
fprintf('Done Meshing.\n')

% Take a look. 
figure; hold on;    % this is iso2mesh function titled plotmesh, but renamed to not be confused with the nirfast plotmesh
iso2mesh_plotmesh(node,elem(elem(:,5)==4,:),'z<92','FaceColor',[0.35 0.35 0.35],'EdgeAlpha',0.6) %%brain
iso2mesh_plotmesh(node,elem(elem(:,5)==3,:),'z<92','FaceColor',[0.2 0.6 1],'EdgeAlpha',0.6) %%csf
iso2mesh_plotmesh(node,elem(elem(:,5)==2,:),'z<92','FaceColor',[1 1 0.9],'EdgeAlpha',0.6) %%bone
iso2mesh_plotmesh(node,elem(elem(:,5)==1,:),'z<92','FaceColor',[1 0.8 0.7],'EdgeAlpha',0.6)%% scalp

% Remove dangling elements and repair
elem = elem(elem(:, 5)>0, :);
[node2, elem2] = meshcheckrepair(node, elem, 'isolated');
elem3 = meshreorient(node2, elem2(:, 1:4)); % fix the orientation, such that simbio doesn't crash
elem3(:, 5) = elem2(:, 5);
% Save in tet format
savetetgenele(elem3,'icbm152_09c/mesh_b2m_dense.ele');
savetetgennode(node2,'icbm152_09c/mesh_b2m_dense.node');
% save('icbm152_09c/pial', 'pial_n', 'pial_f');

% Load in fieldtrip
clear;
eegmesh = ft_read_headshape('icbm152_09c/mesh_b2m_dense.ele', 'format', 'tetgen_ele');
eegmesh.tet = eegmesh.tet + 1;  % fix indexing issue
% Construct a headmodel using simbio
fprintf('Constructing a simbio FEM model...\n')
cfg = [];
cfg.method = 'simbio';  % use "simbio" for FEM forward modeling
cfg.conductivity = [1, 1/80, 5, 1]; % resistivity for each layer. As long as the ratios are right
eegvol = ft_prepare_headmodel(cfg, eegmesh);

% Place the electrodes
fprintf('Placing electrodes...\n')
elec = ft_read_sens('standard_1020.elc');   % load the default electrode positions
label32 = {'fp1', 'fp2', 'f7', 'f3', 'fz', 'f4', 'f8', 'ft9', 'fc5', 'fc1', 'fc2', 'fc6', ...
    'ft10', 't7', 'c3', 'cz', 'c4', 't8', 'tp9', 'cp5', 'cp1', 'cp2', 'cp6', 'tp10', 'p7', 'p3', 'pz', 'p4', 'p8', 'o1', 'oz', 'o2'};
[C, idx, ~] = intersect(lower(elec.label), lower(label32));
if length(C) ~= length(label32)
    error('Label mismatch!');
end
idx = sort(idx); % select the desired subset of the electrodes
elec.chanpos = elec.chanpos((idx), :);
elec.elecpos = elec.elecpos((idx), :);
elec.chantype = elec.chantype((idx), :);
elec.chanunit = elec.chanunit((idx), :);
elec.label = elec.label((idx), :);

homogenous = [1,0,0,97;0,1.05,0,131.25;0,0,1,85;0,0,0,1];   % Manually determined values for icbm152 brain
elec_aligned = ft_transform_sens(homogenous, elec); % manual registration of electrodes to the scalp, to roughly move them close
elec_aligned.homogeneous = homogenous;
cfg = [];
cfg.method = 'project'; % project the electrode locations on the scalp
cfg.headshape = eegmesh;
cfg.elec = elec_aligned;
elec_aligned = ft_electroderealign(cfg);
% Take a look
figure, hold on;
ft_plot_mesh(eegmesh);
camlight;
ft_plot_sens(elec_aligned, 'style', 'r');
fprintf('Finished placing electrodes.\n')

% Generate leadfield matrix
headmodel = [];
headmodel.elec = elec_aligned;
headmodel.vol = eegvol;

% Save the headmodel
save('icbm152_09c/eegmodel_dense', 'headmodel', '-v7.3')
clear

% Now load to NIRFAST (modified from auto generated code)
fprintf('Loading mesh to NIRFAST...\n')
solidmesh2nirfast('icbm152_09c/mesh_b2m_dense.ele','/home/jiaming/Dropbox (CMU Biophotonics Lab)/CMU Biophotonics/Users/Jiaming/Research/Codes/Multi/icbm152_09c/nirs_mesh_dense','stnd');
nirs_mesh = load_mesh('icbm152_09c/nirs_mesh_dense');
% Place some optodes
% nirs_mesh.link =[ 1 1 1; 1 2 1; 1 6 1; 1 7 1; 2 2 1; 2 3 1; 2 7 1; 2 8 1; 3 3 1; 3 4 1; 3 8 1; 3 9 1; 4 4 1; 4 5 1; 4 9 1; 4 10 1;];
% nirs_mesh.source.coord =[59.0552 115.0496 173.7791;81.1611 116.4031 181.8012;106.4249 117.7565 182.4198;129.4331 118.2076 174.5855];    % not carefully designed, should update later
% nirs_mesh.link = [1 1 1; 1 2 1; 2 1 1; 2 3 1; 3 2 1; 3 4 1; 4 3 1; 4 4 1; 5 1 1; 5 2 1; 5 3 1; 5 4 1; ];
% nirs_mesh.source.coord = [128.0188 130.8028 172.0269; 101.6362 106.0456 185.0513; 150.3915 105.6647 164.1516; 127.9184 88.5242 176.8405; 126.7757 104.522 178.3064];
% nirs_mesh.source.coord = [101.2777 130.9334 180.3975; 148.8679 131.5659 158.8562; 100.4935 89.6669 184.5726; 151.1533 88.5242 162.877; 126.7757 104.522 178.3064];
nirs_mesh.link = [1 1;1 2; 2 2; 2 3; 3 3; 3 4; 4 1; 4 2; 4 5; 4 6;5 2; 5 3; 5 6; 5 7; 6 3; 6 4; 6 7; 6 8;7 5; 7 6; 8 6; 8 7; 9 7; 9 8];
nirs_mesh.link = [nirs_mesh.link, ones(length(nirs_mesh.link), 1)];
nirs_mesh.source.coord = [99.6536 139.4106 178.1236;
130.4069 141.517 167.1399;
155.2622 140.2532 149.7199;
100.0749 108.6574 185.0202;
132.9345 110.7638 174.585;
156.9473 110.7638 158.3452;
100.0749 81.6957 182.5371;
137.9898 82.5382 171.014;
162.0026 82.5382 150.0803];
nirs_mesh.source.num = (1:size(nirs_mesh.source.coord,1))';
nirs_mesh.source.fwhm = zeros(size(nirs_mesh.source.coord,1),1);
nirs_mesh.source.fixed =0;
nirs_mesh.source.distributed =0;
% nirs_mesh.meas.coord =[47.3256 125.877 163.1881;67.1757 126.3281 174.7211;92.4396 127.6816 181.5611;117.7034 129.9373 176.6471;142.9673 131.2907 162.8498;49.5813 98.8086 170.1871;70.3337 103.32 179.6286;93.3419 105.5757 185.0013;120.8614 105.1246 179.9774;144.3207 104.2223 168.537];
% nirs_mesh.meas.coord = [101.2777 130.9334 180.3975; 148.8679 131.5659 158.8562; 100.4935 89.6669 184.5726; 151.1533 88.5242 162.877;];
nirs_mesh.meas.coord = [86.1728 126.7723 181.0293;
117.7686 127.1936 176.8518;
145.5728 128.0362 162.8086;
167.058 126.7723 141.6222;
85.7515 94.7553 184.8744;
121.56 97.7042 179.9752;
149.3643 97.2829 165.3878;
169.5856 97.7042 144.7947];
nirs_mesh.meas.num = (1:size(nirs_mesh.meas.coord,1))';
nirs_mesh.meas.fixed =0;
nirs_mesh.ri = 1.4 * ones(size(nirs_mesh.ri));  % refractive index
save_mesh(nirs_mesh,'/home/jiaming/Dropbox (CMU Biophotonics Lab)/CMU Biophotonics/Users/Jiaming/Research/Codes/Multi/icbm152_09c/nirs_mesh_dense');
clear nirs_mesh
nirs_mesh = load_mesh('/home/jiaming/Dropbox (CMU Biophotonics Lab)/CMU Biophotonics/Users/Jiaming/Research/Codes/Multi/icbm152_09c/nirs_mesh_dense');   % Need to load again to move optodes to correct positions

% Get standard Jacobian first
% 750nm, data from Eggebrecht et al. 2012
fprintf('Calculating Jacobian for 750nm...\n')
nirs_mesh.mua(nirs_mesh.region == 1) = 0.0170;      % mm^-1
nirs_mesh.mua(nirs_mesh.region == 2) = 0.0116;
nirs_mesh.mua(nirs_mesh.region == 3) = 0.004;
nirs_mesh.mua(nirs_mesh.region == 4) = 0.0180;      % assume all gray matter
% nirs_mesh.mua(nirs_mesh.region == 5) = 0.0167;
nirs_mesh.mus(nirs_mesh.region == 1) = 0.74;        % mm^-1; mus_prime
nirs_mesh.mus(nirs_mesh.region == 2) = 0.94;
nirs_mesh.mus(nirs_mesh.region == 3) = 0.3;
nirs_mesh.mus(nirs_mesh.region == 4) = 0.8359;
% nirs_mesh.mus(nirs_mesh.region == 5) = 1.1908;
J_750 = jacobian(nirs_mesh, 0);

% 850nm, data from Eggevrecht et al. 2012
fprintf('Calculating Jacobian for 850nm...\n')
nirs_mesh.mua(nirs_mesh.region == 1) = 0.0190;      % mm^-1
nirs_mesh.mua(nirs_mesh.region == 2) = 0.0139;
nirs_mesh.mua(nirs_mesh.region == 3) = 0.004;
nirs_mesh.mua(nirs_mesh.region == 4) = 0.0192;
% nirs_mesh.mua(nirs_mesh.region == 5) = 0.0208;
nirs_mesh.mus(nirs_mesh.region == 1) = 0.64;        % mm^-1; mus_prime
nirs_mesh.mus(nirs_mesh.region == 2) = 0.84;
nirs_mesh.mus(nirs_mesh.region == 3) = 0.3;
nirs_mesh.mus(nirs_mesh.region == 4) = 0.6726;
% nirs_mesh.mus(nirs_mesh.region == 5) = 1.0107;
J_850 = jacobian(nirs_mesh, 0);

% Now get the spectral Jacobian; origional is "del OD/del mua", convert to "del OD/del HbO (HbD)"
% [HbO_lambda1 HbD_lambda1]
% [HbO_lambda2 HbD_lambda2]
excoeff = importdata('excoef.txt');
excoeff = excoeff.data;
epsilon = [excoeff(find(excoeff(:,1)==750), 2:3); excoeff(find(excoeff(:,1)==850), 2:3)];
jacobian_full = [J_750.complete*epsilon(1,1), J_750.complete*epsilon(1,2); J_850.complete*epsilon(2,1), J_850.complete*epsilon(2,2)];
% clear J_750 J_850 excoeff epsilon

% Find nodes corresponding to GM and use them as sources
% In this part, I'm limiting only the nodes on the outer surface of the brain
% layer to be the neural source voxels. This is optional, depending on the
% usage of the forward modeling
fprintf('Trying to locate surface nodes...')
% load('icbm152_09c/pial');
% idx_dipole = zeros(length(pial_n), 1);
% for i = 1:length(pial_n)
%     [~, idx_dipole(i)] = min(sum((nirs_mesh.nodes - pial_n(i, :)).^2, 2));
% end
idx_brain = find(nirs_mesh.region == 4);
in = zeros(size(nirs_mesh.elements));
for i = 1:numel(nirs_mesh.elements) % mark if the nodes of the tetrahedrons are within the brain layer
    in(i) = ~all(idx_brain - nirs_mesh.elements(i));
end
inside = find(sum(in, 2) == 4); % a tetrahedron is considered to be in the brain layer only if all for nodes are in the brain layer
elem_brain = nirs_mesh.elements(inside, :);
openface = volface(elem_brain); % get the outer surface of the brain layer
idx_dipole = unique(openface(:)); % define the surface points as the neural sources
jacobian_nirs = jacobian_full(:, [idx_dipole; idx_dipole+length(nirs_mesh.bndvtx)]);
% jacobian_nirs = jacobian_full(:, [idx_brain; idx_brain+length(nirs_mesh.bndvtx)]);

% Back to eeg, calculate the leadfield matrix
fprintf('Calculating leadfield matrix...\n')
load('icbm152_09c/eegmodel_dense');

% define where the source voxels are. Here using the same as in NIRS, but
% can be something else. This definition of sources is not optional
headmodel.grid.pos = headmodel.vol.pos(idx_dipole,:);
headmodel.grid.inside = 1:length(idx_dipole);
% headmodel.grid.pos = pial_n;
% headmodel.grid.inside = 1:length(pial_n);

lf = ft_prepare_leadfield(headmodel);
% Calculate normals for each node in the mesh
fprintf('Calculating normals...\n');
brainmesh = [];
tmp = openface(:);
for i=1:length(tmp)
    tmp(i) = find(idx_dipole == tmp(i), 1);
end
openface = reshape(tmp, size(openface));
brainmesh.pos = nirs_mesh.nodes(idx_dipole, :);
brainmesh.tri = openface;

triangulated = triangulation(brainmesh.tri, brainmesh.pos);
normals = vertexNormal(triangulated);
% Orientation-constrained leadfield matrix
% Coming out of FieldTrip are dipoles with x,y,z components. Now project
% them onto the direction of the normals calculated in the previous step
fprintf('Applying orientation constraints to dipoles...\n')
L = zeros(length(lf.label), length(lf.leadfield));
for i = 1:size(L, 2)
    L(:,i) =  lf.leadfield{i} * normals(i, :)';
end
if any(isnan(L(:)))
	warning('NaN values found in the leadfield matrix were set to zero');
	L(isnan(L)) = 0;
end
fprintf('Done!\n');

% Generate some random data and take a look
data = 5 * randn(length(idx_dipole), 1);
% colors = nan(length(brainmesh.pos), 3);
colors = vals2colormap(data, 'jet');
figure, ft_plot_mesh(brainmesh, 'vertexcolor', colors);

% A fancy way to visualize (projecting onto a standard brain)
% cfg = [];     % Might be incorrect
% cfg.method = 'surface';
% cfg.funparameter = 'pow';
% 
% tmp = [];
% tmp.pos = visualize_mesh.pos;
% tmp.pow = data;
% ft_sourceplot(cfg, tmp);

fprintf('Finished generating visualization mesh.\n')

% save
save('icbm152_09c/eegmodel_dense', 'headmodel', 'lf', 'L', '-v7.3')
save('icbm152_09c/nirsmodel_dense.mat', 'nirs_mesh', 'J_750', 'J_850', 'jacobian_nirs', 'jacobian_full', '-v7.3');
save('icbm152_09c/shared_dense.mat', 'idx_dipole', 'brainmesh');
