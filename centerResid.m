function r = centerResid(phi,f,h)

M = size(phi,1)-2; N = size(phi,2)-2;
r = zeros(M+2,N+2);

%residual calculation
for j = 2:N+1
    
    for i = 2:M+1
        
        r(i,j) = f(i,j)-((phi(i-1,j)-2*phi(i,j)+phi(i+1,j))/h^2+(phi(i,j-1)-2*phi(i,j)+phi(i,j+1))/h^2);
        
    end
    
end

end

