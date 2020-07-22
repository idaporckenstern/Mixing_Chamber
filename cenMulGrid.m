function phi = cenMulGrid(phi,f,h)

M = size(phi,1)-2; N = size(phi,2)-2;
phi = centerGS(phi,f,h,1);
   
if mod(M,2) == 0 && mod(N,2) == 0
    rh = centerResid(phi,f,h);
    r2h = cenRestrict(rh);
    e2h = zeros(M/2+2,N/2+2);
    e2h = cenMulGrid(e2h,r2h,2*h);
    eh = cenProlong(e2h);
    phi = phi+eh;
    %CHANGE BOUNDARY CONDITIONS FOR EACH PROBLEM
    phi = bcGS(phi);
    phi = centerGS(phi,f,h,1);
end

end