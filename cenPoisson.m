function phi = cenPoisson(phi,f,h,niterMax,epsilon)

r = 1+epsilon;
n = 0;
M = size(phi,1)-2;
global Mfine
Mfine = M;

while n < niterMax && r > epsilon
    
    phi = cenMulGrid(phi,f,h);
    r = max(max(abs(centerResid(phi,f,h))));
    
    n = n+1;
    
end

end
