function eta = estimate_eta(hi,lo)
lo(isnan(lo))  = 1e-6;
n = numel(hi);
k = randperm(n,round(0.5*n));
hi_k = hi(k);
hi_min = min(hi,[],'all');
k((hi_k>0.2*hi_min)|(hi_k>=0)) = [];

lo_k = lo(k);
k((lo_k<1e-9)) = [];
eta=-sum(hi(k),'all')/sum(lo(k),'all');

