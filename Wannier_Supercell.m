clc;
fprintf('Importing Chebfun')
addpath('chebfun')
fprintf('   Done\n')

fprintf('==============\nSetup...')
% Set up the Wannier basis
folder = 'hBN';

basis = Wannier(folder);


% Then the lattice and cutoff patch
L = Lattice([2.1755  2.1755;  -1.0877  1.0877]);
patch = Patch(L, 3);

% Now the perturbation lattice and potential
p=3; q=5; n=7;
T = Lattice([2*p+q -sqrt(3)*q; sqrt(3)*q 2*p+q]/(2*n)*L.L);
h = 3.3;

K  = importdata('hBN/BN-total-pot/g2d.dat')';
z       = importdata('hBN/BN-total-pot/z.dat');
VG_z    = importdata('hBN/BN-total-pot/VG_z.dat');
tol = 1e-2;

Values = VG_z(:,1:2:end) + 1i*VG_z(:,2:2:end);
I = sum(conj(Values).*Values) > tol^2;
K = K(:,I);
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
Atomic_Cutoff = 2;
h = zeros(basis.size,basis.size,patch.size);
s = zeros(basis.size,basis.size,patch.size);
load('hBN/hBN_wan_5band.mat');
O = [0;0;0];
patch_O = Patch(patch.lattice, Atomic_Cutoff, O);

for m = 1:patch.size
    R = [patch.positions_on_lattice(:,m); 0];
    [~,ind] = ismember(R', ham_r, 'rows');
    s(:,:,m) = basis.Overlap(O, R);
    h(:,:,m) = ham_real(:,:,ind) + 1i * ham_imag(:,:,ind);
end
fprintf('   Done\n')

% Compute D0
fprintf('==============\nComputing D0...')
gap_position = 4; Nq = 30;
D0 = PeriodicDensityMatrix(L, basis, patch, h, s, gap_position, Nq);
fprintf('   Done\n')
% %%
% clf;
% surf(squeeze(D0.Qgrid(1,:,:)), squeeze(D0.Qgrid(2,:,:)), squeeze(D0.eigenvalues_q(1,:,:)))
% axis equal
% hold on
% for j=2:basis.size
%     surf(squeeze(D0.Qgrid(1,:,:)), squeeze(D0.Qgrid(2,:,:)), squeeze(D0.eigenvalues_q(j,:,:)))
% end


% Compute D1
fprintf('==============\nComputing D1...')
D1 = PerturbationDensityMatrix(D0, patch, V, Nq);
fprintf('   Done\n')
%%

% Construct superlattice matrices
[super_basis,X] = S.Basis_Setup(basis);
[super_s, super_patch] = S.Change(basis, patch, s);
super_h = S.Change(basis, patch, h);
% Assuming that V is LS-Periodic
super_w = squeeze(sum(S.Change_W(basis, patch, V.K_list, basis.W(V, patch)),3));
%super_w = squeeze(sum(super_basis.W(V, super_patch),3));

% Now the perturbation error tests
fprintf('==============\nComputing the perturbation error for different potential scalings...')
N = 10;
t = 100.^((0:N-1)/(N-1)-1);
O = [0;0];
B0 = zeros(super_basis.size, super_basis.size, super_patch.size);
B1 = zeros(super_basis.size, super_basis.size, super_patch.size);
%
parfor r = 1:super_patch.size
    R = super_patch.positions(:,r);
    B0(:,:,r) = S.Supercell_Block(D0, R, O);
    B1(:,:,r) = S.Supercell_Block(D1, R, O);
end

err0 = zeros(1,N); err1 = zeros(1,N);
parfor i = 1:N
    D = PeriodicDensityMatrix(S, super_basis, super_patch, super_h+t(i)*super_w, super_s, S.N*gap_position, Nq);
    for r = 1:super_patch.size
        R = super_patch.positions(:,r);
        B = D.Block(R,O);
        err1(i) = max(err1(i), norm(B - B0(:,:,r) - t(i)*B1(:,:,r)));
        err0(i) = max(err0(i), norm(B - B0(:,:,r)));
    end
end

fprintf('   Done\n')
figure(1); clf;
loglog(t, err0, 'b-'); hold on;
loglog(t, err1, 'r-');
