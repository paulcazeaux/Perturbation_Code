fid = fopen('hBN/wannier90_00001.xsf');
[~] = textscan(fid, 'CRYSTAL\nPRIMVEC\n', 'CommentStyle','#');
tmp = textscan(fid, '%n %n %n\n', 3);
L = [tmp{1} tmp{2} tmp{3}]';
[~] = textscan(fid, 'CONVVEC\n%n %n %n\n%n %n %n\n%n %n %n\n', 'CommentStyle','#');
tmp = textscan(fid, 'PRIMCOORD\n%n %n', 'CommentStyle','#');
num_atoms = tmp{1}; assert(tmp{2} == 1);
tmp = textscan(fid, '%s %n %n %n\n', num_atoms);
atom_names = tmp{1};
atom_positions = [tmp{2}';tmp{3}';tmp{4}'];
tmp = textscan(fid, 'BEGIN_BLOCK_DATAGRID_3D\n3D_field\nBEGIN_DATAGRID_3D_UNKNOWN\n%n %n %n\n', 'CommentStyle','#');
ngx = tmp{1}; ngy = tmp{2}; ngz = tmp{3};
tmp = textscan(fid, '%f %f %f\n', 1);
grid_origin = [tmp{1};tmp{2};tmp{3}];
tmp = textscan(fid, '%f %f %f\n', 1);
a1 = [tmp{1};tmp{2};tmp{3}];
tmp = textscan(fid, '%f %f %f\n', 1);
a2 = [tmp{1};tmp{2};tmp{3}];
tmp = textscan(fid, '%f %f %f\n', 1);
a3 = [tmp{1};tmp{2};tmp{3}];

datagrid = textscan(fid, '%f');