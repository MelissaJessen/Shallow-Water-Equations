function [bx,b] = Build_slopes(x,x0,alphas)
bx = zeros(size(x));
b = zeros(size(x));

if nargout>1
lowtmp = 0;
hightmp = 0;
end
for i = 1:length(alphas)
    if i == 1
        bx(x<=x0(i)) = alphas(1);
        if nargout>1
        b(x<=x0(i))  = x(x<=x0(i))*alphas(1);
        hightmp = b(x<=x0(i));
        end
    elseif i==length(alphas)
        bx(x>x0(end)) = alphas(end);
        if nargout>1
        b(x>x0(end))  = x(x>x0(end))*alphas(end);
        lowtmp = hightmp;
        hightmp = b(x>x0(end));
        b(x>x0(end)) = b(x>x0(end))-(hightmp(1)-lowtmp(end));
        end

    else
        bx(x>x0(i-1) & x<=x0(i)) = alphas(i);
        if nargout>1
        b(x>x0(i-1) & x<=x0(i))  = x(x>x0(i-1) & x<=x0(i))*alphas(i);
        lowtmp = hightmp;
        hightmp = b(x>x0(i-1) & x<=x0(i));
        b(x>x0(i-1) & x<=x0(i)) = b(x>x0(i-1) & x<=x0(i))-(hightmp(1)-lowtmp(end));
        hightmp = b(x>x0(i-1) & x<=x0(i));
        end

    end
        
end
end