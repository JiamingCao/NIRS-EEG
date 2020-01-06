% addpath('inhull/')
function Beta_HbO_proj = Proj_to_cortex(nirs_mesh, hemoglobin, brainmesh)
proj_on_cortex = zeros(size(nirs_mesh.link, 1), 3);
% fill the convex hull with points
xgrid = floor(min(brainmesh.pos(:,1))):ceil(max(brainmesh.pos(:,1)));
ygrid = floor(min(brainmesh.pos(:,2))):ceil(max(brainmesh.pos(:,2)));
zgrid = floor(min(brainmesh.pos(:,3))):ceil(max(brainmesh.pos(:,3)));
[X,Y,Z] = meshgrid(xgrid, ygrid, zgrid);
points_all = [X(:),Y(:),Z(:)];
points_in_hull = points_all(inhull(points_all, brainmesh.pos), :);

for ch = 1:size(nirs_mesh.link, 1)
    ch_center = 0.5*(nirs_mesh.source.coord(nirs_mesh.link(ch, 1), :) + nirs_mesh.meas.coord(nirs_mesh.link(ch, 2), :));
    % intersection with the hull
    [~, idx_hull] = sort(vecnorm(points_in_hull-ch_center, 2, 2));
    inter_hull(ch,:) = points_in_hull(idx_hull(1), :);
    % now find the intersection with the cortical surface
    line_vec = ch_center - inter_hull(ch,:);
    proj_len_on_line = (brainmesh.pos - inter_hull(ch, :)) * line_vec' / norm(line_vec);
    dist_sq = sum((brainmesh.pos - inter_hull(ch, :)).^2, 2) - proj_len_on_line.^2;
    [~, idx_cortex] = sort(dist_sq);
    for i = 1:length(idx_cortex)
        if norm(brainmesh.pos(idx_cortex(i), :) - ch_center) < 50
            proj_on_cortex(ch, :) = brainmesh.pos(idx_cortex(i), :);
            break;
        end
    end
end

HbO_norm = mean(hemoglobin(1:size(hemoglobin, 1)/2, :), 2);
HbO_norm = HbO_norm/max(HbO_norm);
% interpolation
Beta_HbO_proj = zeros(length(idx_cortex), 1);
for voxel = 1:length(idx_cortex)
    w = 1./sum((proj_on_cortex - brainmesh.pos(voxel, :)).^2, 2);
    w(w<1/400) = 0;
    if ~any(w)
        continue;
    end
    if any(isinf(w))
        tmp_idx = find(isinf(w));
        Beta_HbO_proj(voxel) = HbO_norm(tmp_idx);
    else
        Beta_HbO_proj(voxel) = dot(HbO_norm, w)/sum(w);
    end
end

% take a look
% colors =  vals2colormap(Beta_HbO_proj, 'jet', [-1, 1]);
% figure, ft_plot_mesh(brainmesh, 'vertexcolor', colors);
% colormap('jet'), caxis([-1, 1]);

% figure, ft_plot_mesh(brainmesh);
% hold on;
% scatter3(proj_on_cortex(:,1), proj_on_cortex(:,2), proj_on_cortex(:,3), 'filled');
    