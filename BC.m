function [u,v,Y] = BC(u,v,Y,t)

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

% u boundary conditions
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

for j = 1:N+2

    u(1,j) = 0;
    u(M+1,j) = 0;
        
    if y(j) > inStar1 && y(j) < inEnd1
        z = y(j)-inStar1;
        u(1,j) = V1*f(z,Ly/8);
    end
    
    if y(j) > inStar3 && y(j) < inEnd3
        z = y(j)-inStar3;
        u(M+1,j) = -V3*f(z,Ly/8);
    end
    
end

% v boundary conditions
x = linspace(0-h/2,Lx+h/2,M+2);
y = linspace(0,Ly,N+1);

for i = 1:M+2
    
    v(i,1) = 0;
    v(i,N+1) = 0;
        
    if x(i) > inStar2 && x(i) < inEnd2
        z = x(i)-inStar2;
        v(i,1) = V2*f(z,Lx/8);
    end
    
    if (x(i) > outStar1 && (x(i) < outEnd1)) && out1 == 1
        v(i,1) = (4*v(i,2)-v(i,3))/3;
    end
    
    if x(i) > inStar4 && x(i) < inEnd4
        z = x(i)-inStar4;
        v(i,N+1) = -V4*f(z,Lx/8);
    end

    if (x(i) > outStar2 && (x(i) < outEnd2)) && out2 == 1
         v(i,N+1) = (4*v(i,N)-v(i,N-1))/3;
    end
       
end

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

% Y boundary conditions
x = linspace(0-h/2,Lx+h/2,M+2);
y = linspace(0-h/2,Ly+h/2,N+2);

if mod(t,10) < 5
    fY = 1;
else
    fY = 0;
end

for i = 1:M+2
    
    Y(i,1) = Y(i,2);
    Y(i,N+2) = Y(i,N+1);
        
    if x(i) >= inStar2 && x(i) <= inEnd2
        Y(i,1) = 2*fY-Y(i,2);
    end
    
    if x(i) >= inStar4 && x(i) <= inEnd4
        Y(i,N+2) = 2*fY-Y(i,N+1);
    end
    
end

for j = 1:N+2
    
    Y(1,j) = Y(2,j);
    Y(M+2,j) = Y(M+1,j);
        
    if y(j) >= inStar1 && y(j) <= inEnd1
        Y(1,j) = 2*fY-Y(2,j);
    end
    
    if y(j) >= inStar3 && y(j) <= inEnd3
        Y(M+2,j) = 2*fY-Y(M+1,j);
    end        

end

end

