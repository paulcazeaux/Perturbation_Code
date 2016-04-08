function result = W(obj, V, patch)
% calculates (40) from Perturbation Notes
% V is a PerturbingPotential object

    result = zeros(obj.size, obj.size, size(V.K_list,2), patch.size);
    Tr = V.lattice.Lr;
    N = obj.dvol/sqrt(abs(det(Tr)));
    
    % Parameters of the grid
    ngx = size(obj.grid, 2); ngy = size(obj.grid, 3);
	P = [obj.grid(:,2,1,1) obj.grid(:,1,2,1) obj.grid(:,1,1,2)] - repmat(obj.grid(:,1,1,1),1,3);
    
    Z = 1+P(3,3)\(obj.Quad.z-obj.grid(3,1,1,1));
    
    % (x,y) plane Fourier grid shift
   fourier_shift = floor([ngx+1 ngy+1]/2);
        
    for ind1 = 1:obj.size
        for ind2 = 1:obj.size
            R2 = obj.centers(:,ind2);
            R0 = obj.grid(1:2,1,1,1) + obj.centers(1:2, ind2);
            interp_values_2 = interp3(obj.values(:,:,:,ind2), 1:ngx, 1:ngy, Z-R2(3), 'splines', 0);
            for r = 1:patch.size
                R1 = [patch.positions(:,r);0] + obj.centers(:,ind1);
                if norm(R1-R2) > 2*obj.width % We assume no overlap between the Wannier functions
                    result(ind1, ind2, :, r) = 0;
                else
                    R1 = P\R1;
                    R2 = P\R2;
                    product = interp3(obj.values(:,:,:,ind1), (1:ngx)-R1(1)+R2(1), (1:ngy)-R1(2)+R2(2), Z-R1(3), 'splines', 0)...
                            .*interp_values_2;
                    product_fourier2d = circshift(fft2(product),fourier_shift);
                    K_points = -[ngx/(2*pi) 0; 0 ngy/(2*pi)]*P(1:2,1:2)'*V.K_list;
                    K_points(1,:) = K_points(1,:) + fourier_shift(1) + 1;
                    K_points(2,:) = K_points(2,:) + fourier_shift(2) + 1;
                    Fxy = zeros(length(Z),size(K_points,2));
                    for nz = 1:length(Z)
                        Fxy(nz,:) = interp2(product_fourier2d(:,:,nz), K_points(1,:), K_points(2,:), 'linear', 0);
                    end
                    for k = 1:size(V.K_list,2)
%                         for i = 1:ngx
%                             for j = 1:ngy
%                                 Fxy(:,k) = Fxy(:,k) + exp(2*1i*pi*(K_points(1,k)*(i-1)/ngx + K_points(2,k)*(j-1)/ngy))*squeeze(product(i,j,:));
%                             end
%                         end
                        result(ind1, ind2, k,r) = N*exp(1i*R0'*V.K_list(:,k))*sum(V.Evaluate(obj.Quad.z,k).*Fxy(:,k).*obj.Quad.w);
                    end
                end
            end
        end
    end
    
    % Ensure symmetry
    for r1 = 1:patch.size
        [~,r2] = ismember(-patch.positions_on_lattice(:,r1)', patch.positions_on_lattice', 'rows');
        for k = 1:size(V.K_list,2)
            result(:,:,k,r1) = .5*(result(:,:,k,r1) + exp(1i*V.K_list(:,k)'*patch.positions(:,r1))*result(:,:,k,r2).');
            result(:,:,k,r2) = exp(-1i*V.K_list(:,k)'*patch.positions(:,r1))*result(:,:,k,r1).';
        end
    end
end