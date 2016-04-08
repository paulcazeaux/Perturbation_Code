%clc;
fprintf('Importing Chebfun')
addpath('chebfun')
fprintf('   Done\n')

fprintf('==============\nSetup...')
% Set up the Wannier basis
folder = 'graphene';

basis = Wannier(folder);


% Then the lattice and cutoff patch
L = Lattice([sqrt(3)/2 sqrt(3)/2; -1/2 1/2]);
patch = Patch(L, 1.5);

% Now the perturbation lattice and potential
p=3; q=5; n=7;
T = Lattice([2*p+q -sqrt(3)*q; sqrt(3)*q 2*p+q]/(2*n)*L.L);
h = 3.3;

K  = importdata('hBN/BN-total-pot/g2d.dat')';
z       = importdata('hBN/BN-total-pot/z.dat');
VG_z    = importdata('hBN/BN-total-pot/VG_z.dat');
tol = 1e-0;

Values = VG_z(:,1:2:end) + 1i*VG_z(:,2:2:end);
I = sum(conj(Values).*Values) > tol^2;
theta = pi/6;
Rot = [cos(theta) -sin(theta); sin(theta) cos(theta)];

K = T.Lr*round(L.L'*Rot*K(:,I));
Values = Values(:,I);
size_K = size(K,2);
FourierCoefficients = cell(size_K,1);
for k = 1:size_K
	FourierCoefficients{k} = chebfun(Values(:,k), [-10+h 10+h], 'trig');
end

V = Potential(T, K, FourierCoefficients); % Perturbing potential Fourier transformed in the (x,y) plane 

% Here we construct the superlattice
S = SuperLattice(L, [p -q; q p+q]);
fprintf('   Done\n')

% Compute s, h
fprintf('==============\nComputing tight-binding matrices...')
h = zeros(basis.size,basis.size,patch.size);
s = zeros(basis.size,basis.size,patch.size);
load('graphene/graphene_wan_5band.mat');
O = [0;0;0];

for m = 1:patch.size
    R = [patch.positions_on_lattice(:,m); 0];
    [~,ind] = ismember(R', ham_r, 'rows');
    if ind ~=0
        s(:,:,m) = basis.Overlap(O, R);
        h(:,:,m) = ham_real(:,:,ind).' + 1i * ham_imag(:,:,ind).';
    end
end

% Enforce the symmetries for a nice Dirac cone
h =  symmetrize_5bands(patch,h,true);

fprintf('   Done\n')

% Compute D0
fprintf('==============\nComputing D0...')
gap_position = 4; Nq = 99;
D0 = PeriodicDensityMatrix(L, basis, patch, h, s, gap_position, Nq);
fprintf('   Done\n')

clf;
surf(squeeze(D0.Qgrid(1,:,:)), squeeze(D0.Qgrid(2,:,:)), squeeze(D0.eigenvalues_q(1,:,:)))
%axis equal
hold on
for j=2:basis.size
    surf(squeeze(D0.Qgrid(1,:,:)), squeeze(D0.Qgrid(2,:,:)), squeeze(D0.eigenvalues_q(j,:,:)))
end

Gap = squeeze(D0.eigenvalues_q(5,:,:) -  D0.eigenvalues_q(4,:,:));
[gap, I] = min(Gap(:));
[I,J] = ind2sub(size(Gap), I);
gap
D0.Qgrid(:,I,J)/(2*pi)
%%
% Compute D1
fprintf('==============\nComputing D1...')
D1 = PerturbationDensityMatrix(D0, patch, V, Nq);
fprintf('   Done\n')

%% Nq convergence
fprintf('==============\nComputing D0...')
gap_position = 4; Nq = 1;
D0 = PeriodicDensityMatrix(L, basis, patch, h, s, gap_position, Nq);
fprintf('   Done\n')

Nq_list = [1:50 101];
l = length(Nq_list);

O = [0;0];
B = zeros(basis.size, basis.size, l);
err = zeros(l,1);
for i=l:-1:1;
    Nq = Nq_list(i)
    D1 = D1.Change_Nq(Nq);
    B(:,:,i) = D1.Block(O,O);
    err(i) = norm(B(:,:,i) - B(:,:,end))/norm(B(:,:,end));
    err(i)
end
%%
figure(2); clf
loglog(Nq_list, err);
