function SOM(ntrials,swidth,inputangle)
%%% self-organizing map

%%% parameters

Ntrials=ntrials; % default = 10

N_grid=10; % the units form an array of N_grid-by-N_grid
s=swidth; % gaussian width of lateral interaction, default = 1.5
q=0.2; % learning rate


%%% set initial states

th=randn(N_grid)*2*pi; % initial preferred directional angles of all units
Wx0=cos(th); % preferred direction: x component
Wy0=sin(th); % preferred direction: y component

%%% learning loop & visualization

ii=meshgrid(1:N_grid)'; % row numbers mesh
jj=ii'; % column numbers mesh

Wx=Wx0;
Wy=Wy0;
figure

for t=1:Ntrials
        a=rand*inputangle*pi; % input: random directional angle, default = rand*2*pi
        Vx=cos(a);  % x velocity component
        Vy=sin(a);  % y velocity component
        MATCH=cos(th-a); % cosine of difference of angles
        [Rmax,k]=max(MATCH(:)); % maximum response unit and its linear index
        [ic,jc]=ind2sub([N_grid,N_grid],k); % matrix index of unit with max response
        H=exp(-((ii-ic).^2+(jj-jc).^2)/(2*s^2)); % gaussian falloff from the max response unit
        
        clf
        quiver(Wx,Wy)
        hold on
        plot(jc,ic,'ro')    % unit with max response
        plot([jc, jc+Vx],[ic,ic+Vy],'r-') % direction of input
        title(['Trial # ' num2str(t)])
        drawnow
        
        Wx=Wx+q*(H.*(Vx-Wx));   % update
        Wy=Wy+q*(H.*(Vy-Wy));
        W=sqrt(Wx.^2+Wy.^2);    % new vector length
        Wx=Wx./W; % normalize vector length to 1
        Wy=Wy./W;   % normalize vector length to 1
        th=angle(Wx+sqrt(-1)*Wy); %  angles of updated weights
end