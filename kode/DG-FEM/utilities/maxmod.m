function mfunc = maxmod(v)
% Purpose: Implement the maxmod function.
% Hesthaven 2017
m = size(v,1); mfunc = zeros(1,size(v,2));
s = sum(sign(v),1)/m;
ids = find(abs(s)==1);
if(~isempty(ids))
mfunc(ids) = s(ids).*max(abs(v(:,ids)),[],1);
end
return
