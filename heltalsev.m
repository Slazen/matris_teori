function [ev , mult] = heltalsev(A,tol)
if nargin < 2, tol = 0.2; end

ev = eig(A);
mult = zeros(length(ev),1);
flag = 0;
for i=1:length(ev)
   if ev(i) - round(ev(i)) > tol
       %Break call error
       error('Eigenvalue number not a integer');
   end
end

if  ~flag
     for i=1:length(ev)
        mult(i,1) = length(find(abs(ev - ev(i)) < tol ));
     end
end
end