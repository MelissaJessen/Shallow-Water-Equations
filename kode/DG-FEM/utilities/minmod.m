% function psi=minmod(v)
% hesthaven
% function psi = minmod(v)
% % Purpose: Implement the midmod function visa vector
% N = size(v,1); m = size(v,2); %psi = zeros(N,1);
% s = sum(sign(v),2)/m; ids = find(abs(s)==1);
% 
% if (~isempty(ids))
% psi(ids) = s(ids) .* min(abs(v(ids,:)),[],2);
% end
% 
% return;

function mfunc = minmod(v)
% Purpose: Implement the midmod function v is a vector
m = size(v,1); mfunc = zeros(1,size(v,2));
s = sum(sign(v),1)/m;
ids = find(abs(s)==1);
if(~isempty(ids))
mfunc(ids) = s(ids).*min(abs(v(:,ids)),[],1);
end
return;