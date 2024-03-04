function [e] = minimize(e)

N = size(e,1);

e0 = picomax.eps.sym.initialize(N).';


% for n=N:-1:1
%     for m=n:-1:1
%         s1 = solve(e(m,n)==e0(m,n),e(m,n));
%         s2 = solve(e(m,n)==e0(m,n),e0(m,n));
%         if s1~=0 && s2~=0
%             e = subs(e,s2,s1);
%         end
%     end
% end


unique_e = unique(e);
for i=1:numel(unique_e)
    j = find(e==unique_e(i));
    e(j) = e0(j(1));
end

