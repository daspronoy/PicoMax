function [] = generate(o)

if isempty(o.num.QNUM) || isempty(o.num.QPATH)
    o.dat.qi = 0;
    return
end


qpoints = regexp(o.num.QPATH,'([^ ,:]*)','tokens');
qpoints = vertcat(qpoints{:});
% Nq = numel(qpoints);
Nq = sum(o.num.QNUM)+1;
qvector = zeros(3,Nq);

cnt = 1;
for i=1:numel(qpoints)-1
    q1 = o.highsym_point(qpoints{i},o.sys.crystal);
    q2 = o.highsym_point(qpoints{i+1},o.sys.crystal);
    dq = (q2-q1)/o.num.QNUM(i);
    for j=1:o.num.QNUM(i)
        qvector(:,cnt) = q1+dq*(j-1);
        cnt = cnt+1;
    end
end
qvector(:,cnt) = q2;
o.num.QVEC = qvector;
end

