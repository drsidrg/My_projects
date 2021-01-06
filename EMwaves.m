%%
%% Siddharth Raj Gupta
%%

clear all;
clc;

xdim=200; % Grid Dimension in x 
ydim=200; % Grid Dimension in y
T=3500; %Total no of time steps
x_s=100; %x cordinate of the source 
y_s=100; %y cordinate of the source
S=1/sqrt(2); % stability factor
 
epsilon0= 8.85*10^-12; % Permitivity of the free space
mu0=4*pi*10^-7; %Permiablity of free space
c=3*10^8; %speed fo light
epsilon=epsilon0*ones(xdim,ydim); %permitivity vector
mu=mu0*ones(xdim,ydim); %permiability vector

delta=1e-6; %spatial step sizes 
delta_t=S*delta/c; %temporal step sizes

Ez=zeros(xdim,ydim); % Initialization of Ez
Hy=zeros(xdim,ydim); % Initialization of Hy
Hx=zeros(xdim,ydim); % Initialization of Hx

sigma=4e-4*ones(xdim,ydim); % Initializing electric conductivity 
sigma_star=4e-4*ones(xdim,ydim);% Initializing magnetic conductivity 

source=0; %source 


A=((mu-0.5*delta_t*sigma_star)./(mu+0.5*delta_t*sigma_star));%Multiplication factor matrices for H
B=(delta_t/delta)./(mu+0.5*delta_t*sigma_star);%Multiplication factor matrices for H    
C=((epsilon-0.5*delta_t*sigma)./(epsilon+0.5*delta_t*sigma)); %Multiplication factor matrices for E
D=(delta_t/delta)./(epsilon+0.5*delta_t*sigma); %Multiplication factor matrices for E


% Update loop 
for n=1:1:T
    
    % Time dependent boundaries update
    if n<(x_s-2)
        n1=x_s-(n+1);
    else
        n1=1;
    end
    if n<xdim-(1+x_s)
        n2=x_s+n;
    else
        n2=xdim-1;
    end
    if n<(y_s-2)
        n11=y_s-(n+1);
    else
        n11=1;
    end
    if n<ydim-(1+y_s)
        n21=y_s+n;
    else
        n21=ydim-1;
    end
    
    %Hy and Hx update
    Hy(n1:n2,n11:n21)=A(n1:n2,n11:n21).*Hy(n1:n2,n11:n21)+B(n1:n2,n11:n21).*(Ez(n1+1:n2+1,n11:n21)-Ez(n1:n2,n11:n21));
    Hx(n1:n2,n11:n21)=A(n1:n2,n11:n21).*Hx(n1:n2,n11:n21)-B(n1:n2,n11:n21).*(Ez(n1:n2,n11+1:n21+1)-Ez(n1:n2,n11:n21));
    
    %Ez update
    Ez(n1+1:n2+1,n11+1:n21+1)=C(n1+1:n2+1,n11+1:n21+1).*Ez(n1+1:n2+1,n11+1:n21+1)+D(n1+1:n2+1,n11+1:n21+1).*(Hy(n1+1:n2+1,n11+1:n21+1)-Hy(n1:n2,n11+1:n21+1)-Hx(n1+1:n2+1,n11+1:n21+1)+Hx(n1+1:n2+1,n11:n21));
    
    % Perfect Electric Conductor boundary condition
    Ez(1:xdim,1)=0;
    Ez(1:xdim,ydim)=0;
    Ez(1,1:ydim)=0;
    Ez(xdim,1:ydim)=0;
    
    Hy(1:xdim,1)=0;
    Hy(1:xdim,ydim)=0;
    Hx(1,1:ydim)=0;
    Hx(xdim,1:ydim)=0;
    
    % Source conditions
        
            Ez(x_s,y_s)=1;
       
        
  
   im = imagesc(delta*(1:1:xdim)*1e+6,(1e+6*delta*(1:1:ydim))',Ez',[-1,1]);colormap(jet);
    title(['Development of Ez with PEC boundary at time (fs) = ',num2str(round(n*delta_t*1e+15))]); 
    xlabel('x (in um)');
    ylabel('y (in um)');
    set(gca,'FontSize');
    colorbar;
   frame(n) = getframe(gcf);
   drawnow;

end
 % create the video writer with 1 fps
  writerObj = VideoWriter('myVideo.avi');
  writerObj.FrameRate = 10;
  open(writerObj);
% write the frames to the video
for n=1:1:T
    % convert the image to a frame
    frame2 = frame(n) ;    
    writeVideo(writerObj, frame2);
end
% close the writer object
close(writerObj);
