classdef Wannier < Basis
    %LCAO Basis class - implementation of LCAO data structure and routines
    
    properties
        size
        centers
        
        Quad
        
        grid
        values
        dvol
        width
        wannier_centres
    end
    
    methods
        result = Evaluation(obj, Grid, R, ind);        
        value = W(obj, V, patch)
    end
    
    methods % Constructor
        function obj = Wannier(folder)
            if nargin > 0
                [obj.size, obj.wannier_centres, ~, ~] = read_centres_file(sprintf('%s/wannier90_centres.xyz', folder));
                file_list = cell(obj.size, 1);
                for i=1:obj.size
                    file_list{i} = sprintf('%s/wannier90_%05i.xsf', folder,i);
                end
                obj.centers = zeros(3,5);

                % We load and process the Wannier functions
                [L, grid_origin, a1_vec, a2_vec, a3_vec, N_xyz] ...
                    = read_wannier_grid_file(file_list{1});
                obj.dvol = abs(det([a1_vec a2_vec a3_vec]))/(prod(N_xyz-1)*abs(det(L)));
                obj.grid=zeros([3 N_xyz]);

                for ind1=1:N_xyz(1)
                    for ind2=1:N_xyz(2)
                        for ind3=1:N_xyz(3)
                            obj.grid(:,ind1,ind2,ind3) = a1_vec*(ind1-1)/(N_xyz(1)-1)+a2_vec*(ind2-1)/(N_xyz(2)-1)+a3_vec*(ind3-1)/(N_xyz(3)-1)+grid_origin;
                        end
                    end
                end
                % We keep only one half of the grid in the z direction to
                % suppress a spurious image
                ngz = N_xyz(3)/2;
                zrange=floor(ngz/2)+2:ceil(3*ngz/2)+2;
                obj.grid = obj.grid(:,:,:,zrange);
                
                for i = 1:obj.size
                    % We load and process the Wannier functions
                    [~, origin, a1, a2, a3, N, ~, ~, datagrid] ...
                        = read_wannier_grid_file(file_list{i});
                    
                    assert(all(N_xyz==N));
                    assert(all(grid_origin==origin));
                    assert(all(a1_vec==a1));
                    assert(all(a2_vec==a2));
                    assert(all(a3_vec==a3));
                    obj.values(:, :, :, i) = datagrid(:,:,zrange);             
                end
                
                obj.Quad.z = squeeze(obj.grid(3,1,1,:));
                obj.Quad.w = 1;
                
                % Finally, we find the largest circle inscribed in the grid
                obj.width = sqrt(min(sum(obj.grid(1:2,1,:,1).^2,1)));
            end
        end
        
        function result = Overlap(obj, R1, R2)
            if all(R1 == R2)
                result = eye(obj.size);
            else
                result = zeros(obj.size);
            end
        end
    end
end

function [m, wannier_centres, atom_positions, atom_names] = read_centres_file(filename)
    fid = fopen(filename);
    info = textscan(fid, '%n\n Wannier centres, written by Wannier90 on%sat%s');
    centres = textscan(fid, 'X %f %f %f');
    wannier_centres = [centres{1}';centres{2}';centres{3}'];
    atoms = textscan(fid, '%s %f %f %f');
    atom_names = atoms{1};
    atom_positions = [atoms{2}';atoms{3}';atoms{4}'];
    m = size(wannier_centres, 2);
    assert(m + size(atom_positions, 2) == info{1});
    fclose(fid);
end

function [L, grid_origin, a1, a2, a3, N_xyz, atom_names, atom_positions, datagrid] = read_wannier_grid_file(filename)
    fid = fopen(filename);
    
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
    N_xyz = [tmp{1} tmp{2} tmp{3}];
    tmp = textscan(fid, '%f %f %f\n', 1);
    grid_origin = [tmp{1};tmp{2};tmp{3}];
    tmp = textscan(fid, '%f %f %f\n', 1);
    a1 = [tmp{1};tmp{2};tmp{3}];
    tmp = textscan(fid, '%f %f %f\n', 1);
    a2 = [tmp{1};tmp{2};tmp{3}];
    tmp = textscan(fid, '%f %f %f\n', 1);
    a3 = [tmp{1};tmp{2};tmp{3}];
    
    if nargout > 8
        datagrid = textscan(fid, '%f');
        datagrid = reshape(datagrid{1}, N_xyz);
    end
    fclose(fid);
end

