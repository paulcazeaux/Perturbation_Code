function res = symmetrize_5bands(patch,h,symmetrize_sublattices)
    res = h;
    
    R1 = [-1 -1 ; 1 0];
    R2 = [0 1 ; -1 -1];

        % 3-way symmetry
        % A<->A and B<->B sublattice terms
    m_list = 1:patch.size;
    while ~isempty(m_list);
        m = m_list(1);
        R = patch.positions_on_lattice(:,m);

        [t1,m1] = ismember((R1*R)', patch.positions_on_lattice', 'rows');
        [t2,m2] = ismember((R2*R)', patch.positions_on_lattice', 'rows');
        [ts,ms] = ismember(-R', patch.positions_on_lattice', 'rows');
        [ts1,ms1] = ismember((-R1*R)', patch.positions_on_lattice', 'rows');
        [ts2,ms2] = ismember((-R2*R)', patch.positions_on_lattice', 'rows');

        b = zeros(4); c = 0;
        if (t1 && t2 && ts && ts1 && ts2) % All the rotated bonds are in the patch
            b = (h(1:4,1:4,m) + h([2 3 1 4],[2 3 1 4],m1) + h([3 1 2 4],[3 1 2 4],m2)...
                    + h(1:4,1:4,ms)' + h([2 3 1 4],[2 3 1 4],ms1)' + h([3 1 2 4],[3 1 2 4],ms2)')/6;
            c = (h(5,5,m) + h(5,5,m1) + h(5,5,m2) ...
                    + h(5,5,ms)'+h(5,5,ms1)' + h(5,5,ms2)')/6;
        end

            res(1:4,1:4,m) = b;
            res(5,5,m) = c;
        if t1
            res(1:4,1:4,m1) = b([3 1 2 4],[3 1 2 4]);
            res(5,5,m1) = c;
        end
        if t2
            res(1:4,1:4,m2) = b([2 3 1 4],[2 3 1 4]);
            res(5,5,m2) = c;
        end
        if ts
            res(1:4,1:4,ms) = b';
            res(5,5,ms) = c';
        end
        if ts1
            res(1:4,1:4,ms1) = b([3 1 2 4],[3 1 2 4])';
            res(5,5,ms1) = c';
        end
        if ts2
            res(1:4,1:4,ms2) = b([2 3 1 4],[2 3 1 4])';
            res(5,5,ms2) = c';
        end


        m_list = m_list(m_list~=m & m_list~=m1 & m_list~=m2 ...
                    & m_list~=ms & m_list~=ms1 & m_list~=ms2);
    end

    % A->B sublattice terms
    m_list = 1:patch.size;
    while ~isempty(m_list);
        m = m_list(1);
        R = patch.positions_on_lattice(:,m);
        [t1,m1] = ismember((R1*R + [1;0])', patch.positions_on_lattice', 'rows');
        [t2,m2] = ismember((R2*R + [0;1])', patch.positions_on_lattice', 'rows');
        [ts,ms] = ismember(-R', patch.positions_on_lattice', 'rows');
        [ts1,ms1] = ismember((-R1*R + [-1;0])', patch.positions_on_lattice', 'rows');
        [ts2,ms2] = ismember((-R2*R + [0;-1])', patch.positions_on_lattice', 'rows');
        b = zeros(4,1);
        if (t1 && t2 && ts && ts1 && ts2) % One of the rotated bonds is not in the patch
            b = (h(1:4,5,m) + h([2 3 1 4],5,m1) + h([3 1 2 4],5,m2)...
                    + h(5,1:4,ms)' + h(5,[2 3 1 4],ms1)' + h(5,[3 1 2 4],ms2)')/6;
        end
            res(1:4,5,m) = b; 
        if t1
            res(1:4,5,m1) = b([3 1 2 4]);
        end
        if t2
            res(1:4,5,m2) = b([2 3 1 4]);
        end
        if ts
            res(5,1:4,ms) = b';
        end
        if ts1
            res(5,1:4,ms1) = b([3 1 2 4])';
        end
        if ts2
            res(5,1:4,ms2) = b([2 3 1 4])';
        end
        m_list = m_list(m_list~=m & m_list~=m1 & m_list~=m2);
    end

    % B->A sublattice terms
    m_list = 1:patch.size;
    while ~isempty(m_list);
        m = m_list(1);
        R = patch.positions_on_lattice(:,m);
        [t1,m1] = ismember((R1*R + [-1;0])', patch.positions_on_lattice', 'rows');
        [t2,m2] = ismember((R2*R + [0;-1])', patch.positions_on_lattice', 'rows');
        [ts,ms] = ismember(-R', patch.positions_on_lattice', 'rows');
        [ts1,ms1] = ismember((-R1*R + [1;0])', patch.positions_on_lattice', 'rows');
        [ts2,ms2] = ismember((-R2*R + [0;1])', patch.positions_on_lattice', 'rows');
        b = zeros(4,1);
        if (t1 && t2 && ts && ts1 && ts2) % One of the rotated bonds is not in the patch
            b = (h(5,1:4,m) + h(5,[2 3 1 4],m1) + h(5,[3 1 2 4],m2)...
                    + h(1:4,5,ms)' + h([2 3 1 4],5,ms1)' + h([3 1 2 4],5,ms2)')/6;
        end

            res(5,1:4,m) = b;
        if t1
            res(5,1:4,m1) = b([3 1 2 4]);
        end
        if t2
            res(5,1:4,m2) = b([2 3 1 4]);
        end
        if ts
            res(1:4,5,ms) = b';
        end
        if ts1
            res(1:4,5,ms1) = b([3 1 2 4])';
        end
        if ts2
            res(1:4,5,ms2) = b([2 3 1 4])';
        end
        m_list = m_list(m_list~=m & m_list~=m1 & m_list~=m2);
    end
    
    % x-axis symmetry
    S = [0 1; 1 0];
    m_list = 1:patch.size;
    while ~isempty(m_list);
        m = m_list(1);
        R = patch.positions_on_lattice(:,m);
        [ts,ms] = ismember((S*R)', patch.positions_on_lattice', 'rows');
        b = zeros(5);
        if ts
            b = (res(:,:,m) + res([2 1 3 4 5],[2 1 3 4 5],ms))/2;
        end
            res(:,:,m) = b;
        if ts
            res([2 1 3 4 5],[2 1 3 4 5],ms) = b;
        end
        m_list = m_list(m_list~=m & m_list~=ms);
    end
    
    
    % Symmetrize A / B pz sublattices ( graphene only !)
    if symmetrize_sublattices
        res_sym = res;
        m_list = 1:patch.size;
        while ~isempty(m_list);
            m = m_list(1);
            R = patch.positions_on_lattice(:,m);
            
            [ts,ms] = ismember(-R', patch.positions_on_lattice', 'rows');
            
            [t1,m1] = ismember((-R + [1;0])', patch.positions_on_lattice', 'rows');
            if t1
                res_sym(1,[3 4 5],m)  = (res(1,[3 4 5],m)+res(1,[3 5 4],m1))/2;
            else
                res_sym(1,[3 4 5],m)  = 0;
            end            
            
            [t2,m2] = ismember((-R + [0;1])', patch.positions_on_lattice', 'rows');
            if t2
                res_sym(2,[3 4 5],m) = (res(2,[3 4 5],m)+res(2,[3 5 4],m2))/2;
            else
                res_sym(2,[3 4 5],m)  = 0;
            end
            
            [t1,m1] = ismember((-R - [1;0])', patch.positions_on_lattice', 'rows');
            if t1
                res_sym([3 4 5],1,m)  = (res([3 4 5],1,m)+res([3 5 4],1,m1))/2;
            else
                res_sym([3 4 5],1,m)  = 0;
            end
            
            [t2,m2] = ismember((-R - [0;1])', patch.positions_on_lattice', 'rows');
            if t2
                res_sym([3 4 5],2,m) = (res([3 4 5],2,m)+res([3 5 4],2,m2))/2;
            else
                res_sym(2,[3 4 5],m)  = 0;
            end
            
            [t12,m12] = ismember((-R + [1;-1])', patch.positions_on_lattice', 'rows');
            if t12
                res_sym(1,2,m) = (res(1,2,m) + res(1,2,m12))/2;
            else
                res_sym(1,2,m) = 0;
            end
          
            [t21,m21] = ismember((-R + [-1;1])', patch.positions_on_lattice', 'rows');
            if t21
                res_sym(2,1,m) = (res(2,1,m) + res(2,1,m21))/2;
            else
                res_sym(2,1,m) = 0;
            end
            
            if ts
                res_sym(1,1,m) = (res(1,1,m) + res(1,1,ms))/2;
                res_sym(2,2,m) = (res(2,2,m) + res(2,2,ms))/2;
                res_sym([3 4 5],[3 4 5],m) = (res([3 4 5],[3 4 5],m)+res([3 5 4],[3 5 4],ms))/2;
            else
                res_sym(1,1,m) = 0;
                res_sym(2,2,m) = 0;
                res_sym([3 4 5],[3 4 5],m) = 0;
            end
            
            res_sym(:,:,ms) = res_sym(:,:,m)';
            m_list = m_list(m_list~=m & m_list~=ms);
        end
        res = res_sym;
    end
end
