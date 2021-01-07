clc; clear all;

%% 
%% paramaters
%%
nu = 0.01;
Lx =1; Ly =1;
rho = 1;
U = 1;
dx = 0.1 ; dy = dx;
omega = 1.6;
gamma =0;
Re =100;
dt = 0.1;
Tolfac = 10^-4;
KE =0;
%% stencil for x
xu = 0:dx:Lx; 
yu = -dy/2:dy:Ly+dy/2;
nux = length (xu); 
nuy = length (yu); 
%%stencil for y 
xv = -dx/2:dx:Lx+dx/2;
yv = 0:dy:Ly;
nvy = length (yv);
nvx = length (xv);
%stencil for pressure
xp = -dx/2:dx:Lx+dx/2;
yp = -dy/2:dy:Ly+dy/2;
npx = length (xp);
npy = length (yp);
%The inner domain
I = length(xu);
J = length (yv);
%Meshgriding
[X1,Y1] = meshgrid (xu,yu);
[X2,Y2] = meshgrid (xv,yv);
X1 = X1'; Y1 = Y1';
X2 = X2'; Y2 = Y2';
 
%%
%% initialise
%% 

u = zeros(nux,nuy);
v = zeros(nvx,nvy);
g = zeros(I-1,J-1);

%%
%% Boundary conditions
%%
time=1;t=0;
while (time>10e-8)
    
    p=zeros(I+1,J+1);

    %Left boundary
    for j = 2:J
        u(1,j) =0;
        v(1,j) = -v(2,j);
    end

    %Right boundary
    for j = 2:J
        u(I,j)=0;
        v(I+1,j) = -v(I,j);
    end

    %Bottom boundary
    for i = 2:I
        u(i,1) = -u(i,2);
        v(i,1) = 0;
    end

    %Top boundary
    for i = 2:I
        u(i,J+1) = 2*U - u(i,J);
        v(i,J) =0;
    end

%% 
%% Marching in F and G
%%
    Udiff =u;
    Vdiff =v;
    duvdyc =u;
    duvdxc =v;
    F =u;
    G =v;
    
    for i = 2:I-1
        for j =2: J
            Udiff = ((u(i+1,j) -2*u(i,j) +u(i-1,j))/(dx*dx)+(u(i,j+1) -2*u(i,j) +u(i,j-1))/(dy*dy));
            Uadv  = (((u(i,j)+u(i+1,j))/2)^2 - ((u(i-1,j)+u(i,j))/2)^2)/dx+ (abs(u(i,j)+u(i+1,j)) * (u(i,j)-u(i+1,j))- abs(u(i-1,j)+u(i,j)) * (u(i-1,j)-u(i,j))) *gamma/(4*dx);
            Umix  = (((v(i,j)+v(i+1,j))/2) * ((u(i,j+1)+u(i,j))/2) - ((v(i,j-1)+v(i+1,j-1))/2) * ((u(i,j-1)+u(i,j))/2))/dy +(abs(v(i,j)+v(i+1,j)) * (u(i,j)-u(i,j+1)) - abs(v(i,j-1)+v(i+1,j-1)) * (u(i,j-1)-u(i,j)))*gamma/(4*dx);
            F(i,j) = u(i,j) +dt*((Udiff)/Re - Uadv -Umix);
        end
    end
    
     for i = 2:I
        for j =2: J-1
            Vdiff = ((v(i+1,j) -2*v(i,j) +v(i-1,j))/(dx*dx)+ (v(i,j+1) -2*v(i,j) +v(i,j-1))/(dy*dy));
            Vadv  = (((v(i,j)+v(i,j+1))/2)^2 - ((v(i,j-1)+v(i,j))/2)^2)/dy +(abs(v(i,j)+v(i,j+1)) * (v(i,j)-v(i,j+1)) - abs(v(i,j-1)+v(i,j)) * (v(i,j-1)-v(i,j))) *gamma/(4*dy);   
            Vmix  = (((v(i+1,j)+v(i,j))/2) * ((u(i,j+1)+u(i,j))/2) - ((v(i-1,j)+v(i,j))/2) * ((u(i-1,j+1)+u(i-1,j))/2))/dx +(abs(u(i,j+1)+u(i,j)) * (v(i,j)-v(i+1,j)) - abs(u(i-1,j+1)+u(i-1,j)) * (v(i-1,j)-v(i,j)))*gamma/(4*dx);
            G(i,j) = v(i,j) +dt*((Vdiff)/Re - Vadv -Vmix);
         end
     end
     
%%
%%PPE source term
%%

for i =2:I
    for j = 2:J
        g(i-1,j-1) = rho/dt*((F(i,j) - F(i-1,j))/dx + (G(i,j)-G(i,j-1))/dy);
    end
end

%%
%% Tolerance
%%

Z= reshape(g,[(I-1)*(J-1),1]);
Tol = norm (Z,Inf)*Tolfac;

%%
%% Initialsing epsilon 
%%

RES = zeros(I,J);
res=1; %pick any value greater than Tol so that it goes inside the loop first
 r=1;
while (res> Tol)   
             for i = 2:I
                 for j = 2:J
                     ew=1;
                     ee=1;
                     en=1;
                     if (i==2)
                         ew=0;
                     end
                     if (i == I)
                         ee=0;
                     end
                     if (j==(J))
                         en=0;
                     end
                     Q = (1/(ee+ew+en+1))*(ee*p(i+1,j) + ew*p(i-1,j) + en*p(i,j+1) + p(i,j-1) -g(i-1,j-1)*dx*dx);
                     p(i,j) = (1-omega)*p(i,j) +omega*Q;
                     RES(i,j) = (ee*(p(i+1,j)-p(i,j)) + ew*(p(i-1,j)-p(i,j)) + en*(p(i,j+1)-p(i,j)) + (p(i,j-1)-p(i,j)))/(dx*dx)-g(i-1,j-1);
                 end
             end
B = reshape (RES,[I*J,1]);
res=norm(B,Inf);
Linf(r)=res;
r=r+1;
end

%%
%% solve for u and v 
%%

for i = 2:I-1
    for j = 2:J
        u(i,j) = F(i,j) -dt/dx*(p(i+1,j)-p(i,j));
    end
end
for i = 2:I
    for j = 2:J-1
        v(i,j) = G(i,j) -dt/dy*(p(i,j+1)-p(i,j));
    end
end
%%
%% Kinetic energy 
%%

t=t+1;
KEsum(t) = sum(sum(KE));
for i = 2:I
    for j = 2:J
        KE(i,j) = (0.5*(u(i,j)^2+v(i,j)^2));
    end
end
time=abs(KEsum(t)- sum(sum(KE)));

%%
%% going to the last time step
%%

if (time<10e-8)
resplot = Linf;
end
  figure(7)
        time = t;
        uplot(1:I,1:J)  = (1/2) * (u(1:I,1:J) + u(1:I,2:J+1));
        vplot(1:I,1:J)  = (1/2) * (v(1:I,1:J) + v(2:I+1,1:J));
        Len             = sqrt(uplot.^2 + vplot.^2 + eps);
        uplot           = uplot./Len;
        vplot           = vplot./Len;
        q               = ones(I-1,J-1);
        q(2:I,2:J)  = reshape(g,I-1,J-1);       
        sx = 0:.05:2;
        sy = 0:.05:2;
        fn = stream2(xu,yv,uplot',vplot',sx,sy);
        clf, streamline(fn);
        axis equal
        axis([0 Lx 0 Ly])
        hold off
        colormap(jet)
        colorbar
        axis([0 Lx 0 Ly]) 
        title({['2D Cavity Flow with Re = ',num2str(Re)];        
        ['time(\itt) = ',num2str(time)]})
        xlabel('Spatial Coordinate (x) \rightarrow')
        ylabel('Spatial Coordinate (y) \rightarrow')
        drawnow;
end

%plots 

figure(1)
plot(KEsum,'-');
title ('Kinetic energy');
figure(2)
plot(resplot,'-');
title ('Linf RES');
% 
figure('units','normalized','position',[0 0.33 .3 .3])
surf(X1,Y1,u)
xlabel('x')
ylabel('y')
set(gca,'fontsize',26)
title('u')
% 
figure('units','normalized','position',[0 0.01 .3 .3])
surf(X2,Y2,v)
xlabel('x')
ylabel('y')
set(gca,'fontsize',26)
title('v')
%
temp = find(xu==0.5);
y_p=[1.0000 0.9766  0.9688 0.9609  0.9531 0.8516  0.7344 0.6172 0.5000 0.4531 0.2813 0.1719  0.1016 0.0703 0.0625 0.0547 0.0000];%y coordinate
 u_re100=[1.0000 0.8412 0.7887 0.7372 0.68717 0.2315 0.0033  -0.1364  -0.2058  -0.2109  -0.1566  -0.1015  -0.0643  -0.04775  -0.0419  -0.0371 0.0000];% Re=100
figure(5);
plot(y_p,u_re100, yu, u(temp,:))
title ('u velocity')
legend('Data','Computational');
%
x_p=[1.0000 0.9688 0.9609 0.9531 0.9453 0.9063 0.8594 0.8047 0.5000 0.2344 0.2266 0.1563 0.0938 0.0781 0.0703 0.0625 0.0000];
v_re100=[0.0000 -0.05906  -0.0739 -0.0886 -0.10313 -0.16914 -0.22445 -0.24533 0.05454 0.17527 0.17507 0.16077 0.12317 0.1089 0.1009 0.0923 0.0000];
temp = find(yv==0.5);
figure(6);
plot(x_p,v_re100, xv, v(:,temp))
title ('v velocity')
legend('Data','Computational');

%% contour plots

      

%% to see the animation put the contour plot inside the while loop
