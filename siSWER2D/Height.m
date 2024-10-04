% Nonlinear height function 
%   Input: 
%     - eta = free surface elevation 
%     - b   = bottom depth 
%   Output: 
%     - height = fluid depth
function H = Height(eta,b)
    
    H = max(0,eta-b);     

end
