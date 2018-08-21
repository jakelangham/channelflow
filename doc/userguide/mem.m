N =    [8   16  32  48 64 96  128];
rk3 =  [1.5 2.6 10  28 64 202 468];
cnab = [1.5 2.1 6.8 19 40 129 299];


% minimum data has 14 fields
% un1 vn1 wn1 pn1,  un vn wn pn,  fn, fgn, hn, fn_1, fgn_1, hn_1

% cnab taussolver saves 6 additional fields: 20 fields total
% un1 vn1 wn1 pn1,  un vn wn pn,  fn, fgn, hn, fn_1, fgn_1, hn_1
% P0, v0, P+, v+, P-, v-

% cnab taussolver saves 3*6 additional fields: 32 fields total
% un1 vn1 wn1 pn1,  un vn wn pn,  fn, fgn, hn, fn_1, fgn_1, hn_1
% (P0, v0, P+, v+, P-, v-) times 3 substeps


% N^3 doubles/field
% 8   bytes/double
% 1/1048576  bytes/MB
% mindata = 14*8/1048576 N^3
mindata = (14*8/1048576)*N.^3;
cnabtau = (20*8/1048576)*N.^3;
rk3tau  = (32*8/1048576)*N.^3;

% mindata = [0.055 0.44 3.5000  12   28   95  224]

loglog(N,rk3, 'kx', N,cnab,'ko', N, rk3tau, 'k:', N, cnabtau, 'k--', ...
    N, mindata, 'k-');

axis([5 150 0.02 1000])
legend('RK3    measured', 'CNAB measured', 'RK3    estimate', 'CNAB estimate', ...
    'minimal estimate', 4)

ylabel('megabytes', 'FontSize', 14);
xlabel('N', 'FontSize', 14);

print -deps2 memory.eps
