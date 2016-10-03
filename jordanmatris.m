function J  = jordanmatris(A,tol)
if nargin < 2 tol = 0.2; end;
[ev, mult] = heltalsev(A,tol);


uniq_ev = unique(round(ev));
cell_index = 1;
J = zeros(sum(unique(mult)));
for i=1:length(uniq_ev)
    r = [];
    
    A_tmp1 = A-(eye(length(A)) .* uniq_ev(i));
    k = 1;
    prev_r = 0;
    while 1
        A_tmp2 = A_tmp1^k;
        r(k) = rank(A_tmp2);
        if r(k) == prev_r
            break;
        end
        prev_r = r(k);
        k = k + 1;
    end
    r = r';
    p = zeros(length(r),1)+length(A);
    p = p - r;
    
    b = zeros(length(p),1);
    
    b(1) = p(1);
    k = 2;
    while k<=length(b)
       b(k) = p(k) - p(k-1);
       k = k + 1;
    end
    
    n = zeros(length(b)+1,1);
    k = 1;
    while k<length(b)
       n(k) = b(k)-b(k+1);
       k = k + 1;
    end
    
    for j = 1:length(n)
       if n(j) ~= 0
           if j == 1
               J(cell_index:cell_index+n(j)-1,cell_index:cell_index+n(j)-1) = diag(uniq_ev(i).*ones(n(j),1));
           else
               tmp = diag(uniq_ev(i).*ones(1,j*n(j))) + diag(ones(1,j*n(j)-1),1);
               if n(j) > 1
                    tmp(j,j+1) = 0;
               end
               J(cell_index:cell_index+(j*n(j))-1,cell_index:cell_index+(j*n(j))-1) = tmp;
           end
           cell_index = cell_index + (j*n(j));
       end
    end
end

end

