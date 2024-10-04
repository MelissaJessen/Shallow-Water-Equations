function [bx,b] = Build_slopes(x,x0,alphas)
% could add much more functionality as x0 will be empty for one pipe
% section etc.

bx = ones(size(x));
% 
% for i = 1:length(x0)
%     if i == 1
%         bx(x<=x0(i)) = bx(x<=x0(i))*alphas(1);
%     elseif i==length(x)
%         bx(x>x0(end)) = alphas(end);
%     else
%         bx(x>x0(i-1) & x<x0(i)) = alphas(i);
%     end
%         
% end
% 
% b = x.*bx;

% b(x>x0(1)) = b(x>x0(1))-

end