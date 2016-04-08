
fid = fopen('hBN/Wannier90_Output/UNK00001.1');
hr1 = fread(fid, 1, 'int32');
ngx = fread(fid, 1, 'ulong');
ngy = fread(fid, 1, 'ulong');
ngz = fread(fid, 1, 'ulong');
k = fread(fid, 1, 'ulong');
nbnd = fread(fid, 1, 'ulong');
assert(hr1 == fread(fid, 1, 'int32'));
fclose(fid);
bloch_values = zeros(ngx*ngy*ngz, nbnd, 361);

for k=1:361
    fid = fopen(sprintf('hBN/Wannier90_Output/UNK%05d.1', k));
    hr1 = fread(fid, 1, 'int32');
    assert(ngx == fread(fid, 1, 'ulong'));
    assert(ngy == fread(fid, 1, 'ulong'));
    assert(ngz == fread(fid, 1, 'ulong'));
    assert(k == fread(fid, 1, 'ulong'));
    assert(nbnd == fread(fid, 1, 'ulong'));
    assert(hr1 == fread(fid, 1, 'int32'));
    for band=1:nbnd
        hr2 = fread(fid, 1, 'int32');
        assert(hr2 == ngx*ngy*ngz*2*8);
        bloch_values(:,band, k) = [1 1i]*fread(fid,[2 ngx*ngy*ngz],'float64');
        assert(hr2 == fread(fid, 1, 'int32'));
    end
    fclose(fid);
end

%%
N = zeros(nbnd,361);
for k=1:361
    for band=1:20
        N(band,k) = sqrt(sum(abs(bloch_values(:,band,k)).^2)/(ngx*ngy*ngz));
    end
end
errorbar(mean(N,2), std(N,0,2))

%%
figure(2)
data = reshape(abs(bloch_values(:,1,1)), [ngx ngy ngz]);
surf([squeeze(data(:,1,:)) squeeze(data(:,1,:))]);
