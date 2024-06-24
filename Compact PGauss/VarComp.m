function [Dq] = VarComp(q)
[a,b]=size(q);

%Finding mean in each dimension
v=mean(q');

%Create a matrix where each column is the mean
M=repmat(v',1,b);

%Apply variance formula
S=(q-M).^2/(b-1);
Dq=sum(S')';
end