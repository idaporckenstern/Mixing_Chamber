function [u,v] = uBC(u,v,t)

M = size(u,1)-1;
N = size(u,2)-2;
Lx = 4;
Ly = Lx;
h = Lx/M;

%inlet locations
inStar1 = Ly/16;
inEnd1 = inStar1+Ly/8;
 
inStar2 = Lx-Lx/16-Lx/8;
inEnd2 = inStar2+Lx/8;

inStar3 = Ly-Ly/16-Ly/8;
inEnd3 = inStar3+Ly/8;

inStar4 = Lx/16;
inEnd4 = inStar4+Lx/8;

%Outlet locations
outStar1 = inStar2-Lx/4-Lx/8;
outEnd1 = outStar1+Lx/8;

outStar2 = inEnd4+Lx/4;
outEnd2 = outStar2+Lx/8;

%inlet condition constants and stuff
f = @(z,L) 6*z/L*(1-z/L);
V1 = 1/sqrt(2); V3 = V1;
V2 = 1; V4 = V2;

%outlet open/close conditions
if t < 10 || mod(t-10,6) < 3
    out1 = 1;
else
    out1 = 0;
end

if t < 10 || mod(t-10,6) >= 3
    out2 = 1;
else
    out2 = 0;
end

% u ghost boundary conditions
x = linspace(0,Lx,M+1);
y = linspace(0-h/2,Ly+h/2,N+2);

for i = 1:M+1
    
    u(i,1) = -u(i,2);
    u(i,N+2) = -u(i,N+1);
        
    if (x(i) > outStar1 && (x(i) < outEnd1)) && out1 == 1
        u(i,1) = u(i,2);
    end
    
    if (x(i) > outStar2 && (x(i) < outEnd2)) && out2 == 1
        u(i,N+2) = u(i,N+1);
    end
    
end

% v ghost boundary conditions
x = linspace(0-h/2,Lx+h/2,M+2);
y = linspace(0,Ly,N+1);

for j = 1:N+1

    v(1,j) = -v(2,j);
    v(M+2,j) = -v(M+1,j);
        
    if y(j) > inStar1 && y(j) < inEnd1
        z = y(j)-inStar1;
        v(1,j) = 2*(V1*f(z,Ly/8))-v(2,j);
    end
    
    if y(j) > inStar3 && y(j) < inEnd3
        z = y(j)-inStar3;
        v(M+2,j) = 2*(-V3*f(z,Ly/8))-v(M+1,j);
    end    
    
end

end