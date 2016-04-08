K_list  = importdata('hBN/BN-total-pot/g2d.dat')';
z       = importdata('hBN/BN-total-pot/z.dat');
VG_z    = importdata('hBN/BN-total-pot/VG_z.dat');
tol = 0;

Values = VG_z(:,1:2:end) + 1i*VG_z(:,2:2:end);

I = sum(conj(Values).*Values) > tol^2;
K_list = K_list(:,I);
Values = Values(:,I);
size_K = size(K_list,2);
FourierCoefficients = cell(size_K,1);
for k = 1:size_K
	FourierCoefficients{k} = chebfun(Values(:,k), [-10 10], 'equi');
end


figure(1); clf;

P = @(theta) [cos(theta) -sin(theta); sin(theta) cos(theta)];

V = chebfun2(0, [-2 2 -2 2]);
subplot(1,3,3)
hold on
for k = 1:size_K
    V = V + chebfun2(@(x,y) reshape(exp(2*1i*pi*([x(:) y(:)]*K_list(:,k)))*FourierCoefficients{k}(0), size(x)), [-2 2 -2 2]);
    scatter(K_list(1,k), K_list(2,k), '+')
end
V = real(V);
axis equal
subplot(1,3,[1 2])
plot(V);
axis equal
%NB_Potential = Potential(lattice,K_list,FourierCoefficients,Nz);
%%
figure(2); clf;
C = zeros(1,size_K);
for k = 2:size_K
    C(k) =  norm(FourierCoefficients{k});
end
scatter(K_list(1,:), K_list(2,:), [],  C)
axis equal
hold on

scatter(K_list(1,250), K_list(2,250), 'r+')
scatter(K_list(1,251), K_list(2,251), 'r+')
scatter(K_list(1,273), K_list(2,273), 'r+')
scatter(R(1,:), R(2,:))
