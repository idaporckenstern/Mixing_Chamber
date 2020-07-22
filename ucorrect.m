function v = ucorrect(v,u,t)

Lx = 4;
Ly = Lx;
M = size(v,1)-2;
N = size(v,2)-1;
h = Lx/M;
x = linspace(0-h/2,Lx+h/2,M+2);

%inlet locations
inStar2 = Lx-Lx/16-Lx/8;

inStar4 = Lx/16;
inEnd4 = inStar4+Lx/8;

%Outlet locations
outStar1 = inStar2-Lx/4-Lx/8;
outEnd1 = outStar1+Lx/8;

outStar2 = inEnd4+Lx/4;
outEnd2 = outStar2+Lx/8;

q = sum(u(1,2:N+1)*h)-sum(u(M+1,2:N+1)*h)+sum(v(2:M+1,1)*h)-sum(v(2:M+1,N+1)*h);

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

if out1 == 1 && out2 == 1
    Lout = Lx/4;
else Lout = Lx/8;
end

ucorr = q/Lout;
    
for i = 1:M+1

    if (x(i) > outStar1 && (x(i) < outEnd1)) && out1 == 1
        v(i,1) = v(i,1)-ucorr;
    end

    if (x(i) > outStar2 && (x(i) < outEnd2)) && out2 == 1
            v(i,N+1) = v(i,N+1)+ucorr;
    end
        
    end
        
end