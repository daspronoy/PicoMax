function [] = generate(o)

if isempty(o.num.qnum) || isempty(o.num.qpath)
    o.dat.qi = 0;
    return
end


qpoints = regexp(o.num.qpath,'([^ ,:]*)','tokens');
qpoints = vertcat(qpoints{:});
% Nq = numel(qpoints);
nq = sum(o.num.qnum)+1;
qvector = zeros(3,nq);

cnt = 1;
for i=1:numel(qpoints)-1
    q1 = o.highsym_point(qpoints{i},o.sys.crystal);
    q2 = o.highsym_point(qpoints{i+1},o.sys.crystal);
    dq = (q2-q1)/o.num.qnum(i);
    for j=1:o.num.qnum(i)
        qvector(:,cnt) = q1+dq*(j-1);
        cnt = cnt+1;
    end
end
qvector(:,cnt) = q2;
o.num.qvec = qvector;
end

