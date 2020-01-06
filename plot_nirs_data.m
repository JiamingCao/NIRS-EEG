function plot_nirs_data(oxy_deoxy, nirs_mesh, titles, dim, brainmesh)
% Plots oxy and deoxy data, by projecting optodes onto a 2D plane
% oxy_deoxy should be oxy and deoxy concatenated
% nirs_mesh is in NIRFAST format, requiring:
% nirs_mesn.source.coord
% nirs_mesh.meas.coord
% nirs_mesh.link
% titles is cell array to control the titles of the figures

if nargin < 5
    dim = 2;
end

optodes = [nirs_mesh.source.coord; nirs_mesh.meas.coord];
HbO = oxy_deoxy(1:length(oxy_deoxy)/2);
HbD = oxy_deoxy(length(oxy_deoxy)/2+1:end);
link = nirs_mesh.link;
if dim ==2
    % Recenter
    src = nirs_mesh.source.coord - mean(optodes);
    det = nirs_mesh.meas.coord - mean(optodes);
    % Get spherical azimuthal and radius, and rotate back to same plane
    [th, ~, r] = cart2sph(src(:,1), src(:,2), src(:,3));
    [src_x, src_y] = pol2cart(th, r);
    [th, ~, r] = cart2sph(det(:,1), det(:,2), det(:,3));
    [det_x, det_y] = pol2cart(th, r);
else
    src = nirs_mesh.source.coord;
    det = nirs_mesh.meas.coord;
end


% Plot oxy
figure, hold on, title(titles{1})
colors_hbo = vals2colormap(HbO, 'jet', [-max(abs(HbO)), max(abs(HbO))]);
for i=1:size(link, 1)
    if dim == 2
        line([src_x(link(i, 1)), det_x(link(i, 2))], [src_y(link(i, 1)), det_y(link(i, 2))], 'Color', colors_hbo(i, :), 'LineWidth', 5);
    else
        line([src(link(i, 1),1), det(link(i, 2),1)], [src(link(i, 1),2), det(link(i, 2),2)], [src(link(i, 1),3), det(link(i, 2),3)], 'Color', colors_hbo(i, :), 'LineWidth', 5);
    end
end
axis off
colorbar;
colormap('jet');
caxis([-max(abs(HbO)), max(abs(HbO))]);
if dim ==2
    scatter(src_x, src_y, 100, 'ro', 'filled');
    scatter(det_x, det_y, 100, 'bo', 'filled');
else
    scatter3(src(:,1), src(:,2), src(:,3), 100, 'ro', 'filled');
    scatter3(det(:,1), det(:,2), det(:,3), 100, 'bo', 'filled');
    ft_plot_mesh(brainmesh);
end

% Plot deoxy
figure, hold on, title(titles{2})
colors_hbd = vals2colormap(HbD, 'jet', [-max(abs(HbD)), max(abs(HbD))]);
for i=1:size(link, 1)
    if dim == 2
        line([src_x(link(i, 1)), det_x(link(i, 2))], [src_y(link(i, 1)), det_y(link(i, 2))], 'Color', colors_hbd(i, :), 'LineWidth', 5);
    else
        line([src(link(i, 1),1), det(link(i, 2),1)], [src(link(i, 1),2), det(link(i, 2),2)], [src(link(i, 1),3), det(link(i, 2),3)], 'Color', colors_hbd(i, :), 'LineWidth', 5);
    end
end
axis off
colorbar;
colormap('jet');
caxis([-max(abs(HbD)), max(abs(HbD))]);
if dim ==2
    scatter(src_x, src_y, 100, 'ro', 'filled');
    scatter(det_x, det_y, 100, 'bo', 'filled');
else
    scatter3(src(:,1), src(:,2), src(:,3), 100, 'ro', 'filled');
    scatter3(det(:,1), det(:,2), det(:,3), 100, 'bo', 'filled');
    ft_plot_mesh(brainmesh);
end


function rgb = vals2colormap(vals, colormap, crange)
% Take in a vector of N values and return and return a Nx3 matrix of RGB
% values associated with a given colormap
%
% rgb = AFQ_vals2colormap(vals, [colormap = 'jet'], [crange])
%
% Inputs:
% vals     = A vector of values to map to a colormap or a cell array of
%            vectors of values
% colormap = A matlab colormap. Examples: colormap = 'autumn';
%            colormap = 'jet'; colormap = 'hot';
% crange   = The values to map to the minimum and maximum of the colormap.
%            Defualts to the full range of values in vals.
%
% Outputs:
% rgb      = Nx3 matrix of rgb values mapping each value in vals to the
%            corresponding rgb colors.  If vals is a cell array then rgb
%            will be a cell array of the same length
%
% Example:
% vals = rand(1,100);
% rgb = AFQ_vals2colormap(vals, 'hot');
%
% Copyright Jason D. Yeatman, June 2012

if ~exist('colormap','var') || isempty(colormap)
    colormap = 'jet';
end
% Generate the colormap
if ischar(colormap)
    cmap = eval([colormap '(256)']);
else
    % If the colors were provided calculate the number of unique colors
    ncolors = size(colormap,1);
    rp = floor(256./ncolors);
    rm = 256 - rp.*ncolors;
    cmap = [];
    % Repeat each color an equal number of times
    for ii = 1:ncolors
       cmap = vertcat(cmap,repmat(colormap(ii,:),rp,1));
    end
    cmap = vertcat(repmat(colormap(1,:),round(rm./2),1),cmap,repmat(colormap(end,:),round(rm./2),1));
end
%
if ~iscell(vals)
    if ~exist('crange','var') || isempty(crange)
        crange = [min(vals) max(vals)];
    end
    % Normalize the values to be between 1 and 256
    vals(vals < crange(1)) = crange(1);
    vals(vals > crange(2)) = crange(2);
    valsN = round(((vals - crange(1)) ./ diff(crange)) .* 255)+1;
    % Convert any nans to ones
    valsN(isnan(valsN)) = 1;
    % Convert the normalized values to the RGB values of the colormap
    rgb = cmap(valsN, :);
elseif iscell(vals)
    if ~exist('crange','var') || isempty(crange)
        if size(vals{1},1) > size(vals{1},2)
            crange = [min(vertcat(vals{:})) max(vertcat(vals{:}))];
        else
            crange = [min(horzcat(vals{:})) max(horzcat(vals{:}))];
        end
    end
    for ii = 1:length(vals)
        % Normalize the values to be between 1 and 256 for cell ii
        valsN = vals{ii};
        valsN(valsN < crange(1)) = crange(1);
        valsN(valsN > crange(2)) = crange(2);
        valsN = round(((valsN - crange(1)) ./ diff(crange)) .* 255)+1;
        % Convert any nans to ones
        valsN(isnan(valsN)) = 1;
        % Convert the normalized values to the RGB values of the colormap
        rgb{ii} = cmap(valsN, :);
    end
end
return
