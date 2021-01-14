function pointMax = lagrange_2d(S, ri, ci)
X = [(ri-1)^2  (ri-1)*(ci-1)  (ci-1)^2  ri-1  ci-1  1;
     (ri-1)^2  (ri-1)*(ci)    (ci)^2    ri-1  ci    1;
     (ri-1)^2  (ri-1)*(ci+1)  (ci+1)^2  ri-1  ci+1  1;
     (ri)^2    (ri)*(ci-1)    (ci-1)^2  ri    ci-1  1;
     (ri)^2    (ri)*(ci)      (ci)^2    ri    ci    1;
     (ri)^2    (ri)*(ci+1)    (ci+1)^2  ri    ci+1  1;
     (ri+1)^2  (ri+1)*(ci-1)  (ci-1)^2  ri+1  ci-1  1;
     (ri+1)^2  (ri+1)*(ci)    (ci)^2    ri+1  ci    1;
     (ri+1)^2  (ri+1)*(ci+1)  (ci+1)^2  ri+1  ci+1  1;
     ];
 
Y = [S(ri-1, ci-1) S(ri-1, ci) ... 
     S(ri-1, ci+1) S(ri, ci-1) ...
     S(ri, ci) S(ri, ci+1) ...
     S(ri+1, ci-1) S(ri+1, ci) ...
     S(ri+1, ci+1)]';
params = X\Y;

A = [2*params(1) params(2); params(2) 2*params(3)];
b = [-params(4); -params(5)];

pointMax = A\b;
end