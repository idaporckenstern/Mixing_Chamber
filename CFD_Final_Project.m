%% Setting up problem

clear all
format long
format compact
close all
clc

tPlot = [2 5 10 15 20 24 0];
tEnd = 25;
t = 0;
CFL = .8;
M = 128;
N = M;
Lx = 4;
Ly = Lx;
h = Lx/M;
Re = 40;
Sc = 2;
count = 1;
doPlot = 0;
n = 1;

%spacial values for each vector
xu = linspace(0,Lx,M+1);
yu = linspace(0-h/2,Ly+h/2,N+2);

xv = linspace(0-h/2,Lx+h/2,M+2);
yv = linspace(0,Ly,N+1);

xY = xv;
yY = yu;

xp = xY;
yp = yY;

%% Part 2

% initial conditions
u = zeros(M+1,N+2);
v = zeros(M+2,N+1);
Y = zeros(M+2,N+2);
p = zeros(M+2,N+2);

u2 = zeros(M+1,N+2);
v2 = zeros(M+2,N+1);
Y2 = zeros(M+2,N+2);

%set boundary conditions
[u,v,Y] = BC(u,v,Y,t);

while t < tEnd
    
    if t == 0
          for j = 2:N+1

              for i = 1:M

                  uhat(i,j-1) = (u(i,j)+u(i+1,j))/2;

              end

          end

          for j = 1:N

              for i = 2:M+1

                  vhat(i-1,j) = (v(i,j)+v(i,j+1))/2;

              end

          end
    end
    
    dtu = CFL*h/(max(max(abs(2*u)))+max(max(abs(v))));
    dtv = CFL*h/(max(max(abs(2*v)))+max(max(abs(u))));
    dtY = CFL*h/(max(max(abs(uhat)))+max(max(abs(vhat))));
    dtpara = CFL*Re*h^2/4;
    
    dt = min(min(dtu,dtv),min(dtY,dtpara));
    
        %changing dt for plotting
    if (t < tPlot(count)) && (t+dt > tPlot(count))
        dt = tPlot(count)-t;
        count = count+1;
        doPlot = 1;
    end
    
    if t+dt == tPlot(count)
        doPlot = 1;
        count = count+1;
    end
    
    %solving for u^n+1
    for j = 2:N+1
        
        for i = 2:M
            
            u2(i,j) = u(i,j)+dt*(1/Re*(((u(i+1,j)-2*u(i,j)+u(i-1,j))/h^2)+((u(i,j+1)-2*u(i,j)+u(i,j-1))/h^2))- ... 
                ((((u(i+1,j)+u(i,j))/2)^2-((u(i,j)+u(i-1,j))/2)^2)/h)- ... 
                ((((v(i,j)+v(i+1,j))/2)*((u(i,j)+u(i,j+1))/2)-((v(i,j-1)+v(i+1,j-1))/2)*((u(i,j)+u(i,j-1))/2))/h));
        end
        
        
        
    end           
    
    
    %solving for v^n+1
    for j = 2:N
        
        for i = 2:M+1

            v2(i,j) = v(i,j)+dt*(1/Re*(((v(i+1,j)-2*v(i,j)+v(i-1,j))/h^2)... 
                +((v(i,j+1)-2*v(i,j)+v(i,j-1))/h^2))-((((u(i,j)+u(i,j+1))/2)*((v(i,j)+v(i+1,j))/2)... 
                -((u(i-1,j)+u(i-1,j+1))/2)*((v(i,j)+v(i-1,j))/2))/h)-((((v(i,j)+v(i,j+1))/2)^2 ... 
                -((v(i,j)+v(i,j-1))/2)^2)/h));
            
        end
        
    end
    
      %solving for Y^n+1
      for j = 2:N+1
        
        for i = 2:M+1
            
            Y2(i,j) = Y(i,j)+dt*(1/(Re*Sc)*(((Y(i+1,j)-2*Y(i,j)+Y(i-1,j))/h^2) ... 
                +((Y(i,j+1)-2*Y(i,j)+Y(i,j-1))/h^2)))+(-((u(i,j)+u(i-1,j))/2*dt/h)/2*(Y(i+1,j) ... 
                -Y(i-1,j))+((u(i,j)+u(i-1,j))/2*dt/h)^2/2*(Y(i+1,j)-2*Y(i,j)+Y(i-1,j)))+ ... 
                (-((v(i,j)+v(i,j-1))/2*dt/h)/2*(Y(i,j+1)-Y(i,j-1))+((v(i,j)+v(i,j-1))/2*dt/h)^2/2*(Y(i,j+1) ... 
                -2*Y(i,j)+Y(i,j-1)));
            
        end
        
      end
      
    %apply boundary conditions
    u = u2;
    v = v2;
    Y = Y2;
    [u,v,Y] = BC(u,v,Y,t+dt);
    
    R = Y.*(1-Y);
    
    %correct with conservation of mass
    
    v = ucorrect(v,u,t);
   
    %correct with pressure term
    rhs = zeros(M+2,N+2);
      
      for j = 2:N+1
          
          for i = 2:M+1
              
              rhs(i,j) = 1/dt*((u(i,j)-u(i-1,j))/h+(v(i,j)-v(i,j-1))/h);
              
          end
          
      end      
      
      %solve for pressure
      p = cenPoisson(p,rhs,h,100,1e-10);
   
    %correct u velocity
    for j = 2:N+1
        
        for i = 1:M+1
            
            u(i,j) = u(i,j)-dt*(p(i+1,j)-p(i,j))/h;
            
        end
        
    end    
    
    %correct v velocity
    for j = 1:N+1
        
        for i = 2:M+1
            
            v(i,j) = v(i,j)-dt*(p(i,j+1)-p(i,j))/h;
            
        end
        
    end
    
    [u,v] = uBC(u,v,t);
      
      %calculate integrals
      S(n) = 1/(Lx*Ly)*sum(sum(R(2:M+1,2:N+1)*h^2));
      
      uhat = zeros(M,N);
      vhat = zeros(M,N);
      
      for j = 2:N+1
          
          for i = 1:M
              
              uhat(i,j-1) = (u(i,j)+u(i+1,j))/2;
              
          end
          
      end
      
      for j = 1:N
          
          for i = 2:M+1
              
              vhat(i-1,j) = (v(i,j)+v(i,j+1))/2;
              
          end
          
      end
      
      k(n) = 1/(Lx*Ly)*sum(sum((1/2*(uhat.^2+vhat.^2))*h^2));
      
      t4Plot(n) = t;
      
      t = t+dt;
      
    if doPlot == 1
        figure()
        pcolor(xu,yu,u');
        shading interp
        colormap jet
        caxis([-1.5 1.5])
        colorbar 
        xlabel('x','fontsize',18)
        ylabel('y','fontsize',18)
        title(['FTCS values for u t = ',num2str(t)],'fontsize',18)
        drawnow;
   
        figure()
        pcolor(xv,yv,v');
        shading interp
        colormap jet
        caxis([-1.5 1.5])
        colorbar 
        xlabel('x','fontsize',18)
        ylabel('y','fontsize',18)
        title(['FTCS values for v t = ',num2str(t)],'fontsize',18)
        drawnow;
        
        figure()
        pcolor(xY,yY,Y');
        shading interp
        caxis([0 1])
        colormap hot
        colorbar 
        xlabel('x','fontsize',18)
        ylabel('y','fontsize',18)
        title(['FTCS values for Y t = ',num2str(t)],'fontsize',18)
        drawnow;
        
        figure()
        pcolor(xY,yY,R');
        shading interp
        caxis([0 .25])
        colormap jet
        colorbar 
        xlabel('x','fontsize',18)
        ylabel('y','fontsize',18)
        title(['FTCS values for R t = ',num2str(t)],'fontsize',18)
        drawnow;
        doPlot = 0;
        
    end
    
    n = n+1;
    
end

%plot k and S
figure()
plot(t4Plot,k,'LineWidth',3)
xlabel('t','fontsize',18)
ylabel('k(t)','fontsize',18)
title('k as a function of t','fontsize',18)

figure()
plot(t4Plot,S,'LineWidth',3)
xlabel('t','fontsize',18)
ylabel('S(t)','fontsize',18)
title('S as a function of t','fontsize',18)

%% Part 3

M = 32; 
h = Lx/M;
N = M;
n = 1;
k = 0;
S = 0;
t = 0;
count = 1;
check = 0;
n = 1;

%spacial values for each vector
xu = linspace(0,Lx,M+1);
yu = linspace(0-h(count)/2,Ly+h(count)/2,N+2);

xv = linspace(0-h(count)/2,Lx+h(count)/2,M+2);
yv = linspace(0,Ly,N+1);

xY = xv;
yY = yu;

% initial conditions
u = zeros(M+1,N+2);
v = zeros(M+2,N+1);
Y = zeros(M+2,N+2);
p = zeros(M+2,N+2);

u2 = zeros(M+1,N+2);
v2 = zeros(M+2,N+1);
Y2 = zeros(M+2,N+2);

%set boundary conditions
[u,v,Y] = BC(u,v,Y,t);

while t < tEnd
    
    if check == 0
          for j = 2:N+1

              for i = 1:M

                  uhat(i,j-1) = (u(i,j)+u(i+1,j))/2;

              end

          end

          for j = 1:N

              for i = 2:M+1

                  vhat(i-1,j) = (v(i,j)+v(i,j+1))/2;

              end

          end
    end
    
    dtu = CFL*h/(max(max(abs(2*u)))+max(max(abs(v))));
    dtv = CFL*h/(max(max(abs(2*v)))+max(max(abs(u))));
    dtY = CFL*h/(max(max(abs(uhat)))+max(max(abs(vhat))));
    dtpara = CFL*Re*h^2/4;
    
    dt = min(min(dtu,dtv),min(dtY,dtpara));
    
        %changing dt for plotting
    if (t < tPlot(count)) && (t+dt > tPlot(count))
        dt = tPlot(count)-t;
        count = count+1;
        doPlot = 1;
    end
    
    if t+dt == tPlot(count)
        doPlot = 1;
        count = count+1;
    end
    
    %solving for u^n+1
    for j = 2:N+1
        
        for i = 2:M
            
            u2(i,j) = u(i,j)+dt*(1/Re*(((u(i+1,j)-2*u(i,j)+u(i-1,j))/h^2)+((u(i,j+1)-2*u(i,j)+u(i,j-1))/h^2))- ... 
                ((((u(i+1,j)+u(i,j))/2)^2-((u(i,j)+u(i-1,j))/2)^2)/h)- ... 
                ((((v(i,j)+v(i+1,j))/2)*((u(i,j)+u(i,j+1))/2)-((v(i,j-1)+v(i+1,j-1))/2)*((u(i,j)+u(i,j-1))/2))/h));
        end
        
        
        
    end           
    
    
    %solving for v^n+1
    for j = 2:N
        
        for i = 2:M+1

            v2(i,j) = v(i,j)+dt*(1/Re*(((v(i+1,j)-2*v(i,j)+v(i-1,j))/h^2)... 
                +((v(i,j+1)-2*v(i,j)+v(i,j-1))/h^2))-((((u(i,j)+u(i,j+1))/2)*((v(i,j)+v(i+1,j))/2)... 
                -((u(i-1,j)+u(i-1,j+1))/2)*((v(i,j)+v(i-1,j))/2))/h)-((((v(i,j)+v(i,j+1))/2)^2 ... 
                -((v(i,j)+v(i,j-1))/2)^2)/h));
            
        end
        
    end
    
      %solving for Y^n+1
      for j = 2:N+1
        
        for i = 2:M+1
            
            Y2(i,j) = Y(i,j)+dt*(1/(Re*Sc)*(((Y(i+1,j)-2*Y(i,j)+Y(i-1,j))/h^2) ... 
                +((Y(i,j+1)-2*Y(i,j)+Y(i,j-1))/h^2)))+(-((u(i,j)+u(i-1,j))/2*dt/h)/2*(Y(i+1,j) ... 
                -Y(i-1,j))+((u(i,j)+u(i-1,j))/2*dt/h)^2/2*(Y(i+1,j)-2*Y(i,j)+Y(i-1,j)))+ ... 
                (-((v(i,j)+v(i,j-1))/2*dt/h)/2*(Y(i,j+1)-Y(i,j-1))+((v(i,j)+v(i,j-1))/2*dt/h)^2/2*(Y(i,j+1) ... 
                -2*Y(i,j)+Y(i,j-1)));
            
        end
        
      end
      
    %apply boundary conditions
    u = u2;
    v = v2;
    Y = Y2;
    [u,v,Y] = BC(u,v,Y,t+dt);
    
    R = Y.*(1-Y);
    
    %correct with conservation of mass
    
    v = ucorrect(v,u,t);
   
    %correct with pressure term
    rhs = zeros(M+2,N+2);
      
      for j = 2:N+1
          
          for i = 2:M+1
              
              rhs(i,j) = 1/dt*((u(i,j)-u(i-1,j))/h+(v(i,j)-v(i,j-1))/h);
              
          end
          
      end      
      
      %solve for pressure
      p = cenPoisson(p,rhs,h,100,1e-10);
   
    %correct u velocity
    for j = 2:N+1
        
        for i = 1:M+1
            
            u(i,j) = u(i,j)-dt*(p(i+1,j)-p(i,j))/h;
            
        end
        
    end    
    
    %correct v velocity
    for j = 1:N+1
        
        for i = 2:M+1
            
            v(i,j) = v(i,j)-dt*(p(i,j+1)-p(i,j))/h;
            
        end
        
    end
    
    [u,v] = uBC(u,v,t);
      
      if t >= 10
          check = 1;
          %calculate integrals
          S(n) = 1/(Lx*Ly)*sum(sum(R(2:M+1,2:N+1)*h^2));

          uhat = zeros(M,N);
          vhat = zeros(M,N);

          for j = 2:N+1

              for i = 1:M

                  uhat(i,j-1) = (u(i,j)+u(i+1,j))/2;

              end

          end

          for j = 1:N

              for i = 2:M+1

                  vhat(i-1,j) = (v(i,j)+v(i,j+1))/2;

              end

          end

          k(n) = 1/(Lx*Ly)*sum(sum((1/2*(uhat.^2+vhat.^2))*h^2));

          tInt(n) = t;
          n = n+1;
      end
      
      t = t+dt;
      
end

%T integral
T = 1/15*trapezoidal(S,tInt)
    
%KE integral
KE = 1/15*trapezoidal(k,tInt)