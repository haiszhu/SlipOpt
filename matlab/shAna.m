function shc = shAna(f)
% shAna-compatible wrapper backed by Fortran sht_ana_rotgrid_mex.
% Input f is np-by-nc on GL-uniform grid with np=(p+1)*(2*p).

[d1, d2] = size(f);
p = (sqrt(2*d1 + 1) - 1) / 2;
if abs(p - round(p)) > 1e-12
    error('shAna_mex:bad_size', ...
        'Input rows must satisfy np=(p+1)*(2*p).');
end

p = round(p);
np = d1;
nc = d2;
nshc = (p + 1) * (p + 1);

shc = complex(zeros(nshc, nc));
shc = sht_ana_rotgrid_mex(p, np, nc, f, nshc, shc);
end
