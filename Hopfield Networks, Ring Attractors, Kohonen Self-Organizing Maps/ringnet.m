function ringnet(b,s,boolnoise)
% parameters
N=40;       % number of neurons
tmax=20;    % max time
dt=0.02;    % time step
k=2;        % gain function slope
B=b;     % weight baseline, default = -1.5
S=s;        % weight shift, default = 0


Nt=fix(tmax/dt);    % total number of time steps
dat=zeros(N,Nt);    % data (activity-by-time)

% define weight vector w & weight matrix W
x=linspace(0, 2*pi*(N-1)/N, N); % angular variable from 0 to 2*pi
w=exp(cos(x+S))+B; % weight vector (= 1st row of weight matrix)
W=zeros(N,N);   % weight matrix (circulant)
for i=1:N
    W(i,:)=circshift(w,[0, i-1]);
end

% add noise to weight matrix
if boolnoise == 1
    W=W+randn(size(W));
end
u=3*randn(N,1); % random initial state

% time evolution (Euler method)
for i=1:Nt
    g=1./(1+exp(-k*u));
    u=u+(-u+W*g)*dt;
    dat(:,i)=u;
end

% plotting
figure
clf

subplot(2,2,2)
plot(W(1,:))
title('1st row of W')

subplot(2,2,1)
imagesc(W)
xlabel('Column #')
ylabel('Row #')
axis equal
axis tight
title(['Weight matrix W  (shift S = ' num2str(S) ')'])
colorbar

subplot(2,1,2)
imagesc((1:Nt)*dt, 1:N, dat)
xlabel('Time')
ylabel('Neuron #')
title('Activity')
colorbar
