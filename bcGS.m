function phi = bcGS(phi)

M = size(phi,1)-2;
N = size(phi,2)-2;
global Mfine

%if M == Mfine for HW 7 all BCs are d/dn = 0
    
%     for j = 1:N+2
%     
%     phi(M+2,j) = phi(M+1,j);
%     phi(1,j) = 2-phi(2,j);
%     
%     end
% 
%     for i = 1:M+2
%     
%         phi(i,1) = phi(i,2);
%         phi(i,N+2) = phi(i,N+1);
%     
%     end
%     
% else
    for j = 1:N+2
    
    phi(M+2,j) = phi(M+1,j);
    phi(1,j) = phi(2,j);
    
    end

    for i = 1:M+2
    
        phi(i,1) = phi(i,2);
        phi(i,N+2) = phi(i,N+1);
    
    end
% end
end