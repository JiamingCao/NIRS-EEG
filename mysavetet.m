function mysavetet(filename, nodes, tri)
% Save BEM mesh in Tetgen format. Hopefully it's loadable in NIRFAST
% Modified from surf_to_tetgen.m (Fieldtrip)


%Write the node file, Copied from original
node_filename = [filename '.node'];
precision_str = '%u   %.40e  %.40e  %.40e';

N_nodes = size(nodes,1);
if(nargin > 4)
    N_attributes = size(attributes_nodes,2);
    A = [nodes attributes_nodes];
    for i=1:N_attributes
        precision_str = [precision_str '  %.5e'];
    end
    precision_str = [precision_str '\n'];
else 
    N_attributes = 0;   
    A = nodes;
    precision_str = [precision_str '\n'];
end

fid = fopen(node_filename, 'w');
fprintf(fid, [num2str(N_nodes) ' 3 ' num2str(N_attributes) ' 0\n']);
for i=1:(N_nodes)
    fprintf(fid,precision_str,i,A(i,:));
end
fclose(fid);

% Write the ele file
ele_filename = [filename '.ele'];
num_tets = size(tri, 1);

fid = fopen(ele_filename, 'w');
fprintf(fid, [num2str(num_tets),'  3   0\n']);
for i=1:num_tets
    fprintf(fid, ['\t', num2str(i), '\t%u\t%u\t%u\n'], tri(i, :));
end

