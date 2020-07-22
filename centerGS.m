function phi = centerGS(phi,f,h,niter)

M = size(phi,1)-2;
N = size(phi,2)-2;
%CHANGE BOUNDARY CONDITIONS FOR EACH PROJECT
phi = bcGS(phi);

for k = 1:niter
    
    %loop over j
    for j = 2:N+1
        
        %loop over i
        for i = 2:M+1
            phi(i,j) = 0.25*(phi(i-1,j)+phi(i+1,j)+phi(i,j-1)+phi(i,j+1))-.25*h^2*f(i,j);
        
        end
    
        %CHANGE BOUNDARY CONDITIONS FOR EACH PROJECT
    phi = bcGS(phi); %function updates boundary conditions    
    
    end
    
end

end

