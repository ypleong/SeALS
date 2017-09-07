% Monte Carlo Particle Filter
clear all
%% Set parameters
dim = 2;
x = sym('x',[dim,1]); %do not change
%t = sym('t');
measure = 0;
kpen = (2*pi/3)^2;
fFPE = [x(2);-kpen*(x(1))];
n = [301 301];
    x0 = [0.7, 0.1];
    diagSigma = [0.6 0.5];
dt = 0.04;
finalt = 12;
t = 0:dt:finalt;
nparticles = 2e3;
fDYN = matlabFunction(fFPE,'Vars',x);

% init variables
particle = cell(length(t),1);
%% Generate Initial distribution
particle{1} = mvnrnd(x0,diag(diagSigma),nparticles )';
scatter(particle{1}(1,:),particle{1}(2,:),1)
axis equal
axis([-pi pi -5 5])

%% Propagate particles
time1 = tic;
for k = 2:length(t)
     particle{k} = dynUpdateParticles( particle{k-1}, fDYN,[t(k-1) t(k)]);
end
toc(time1)

%% Compute mean and covariance
xmean = zeros(dim,length(t));
covP  = zeros(dim,dim,length(t));
for k = 1:length(t)
    xmean(:,k) = mean(particle{k}');
    covP(:,:,k) = cov(particle{k}');
end
%% Plot variables
figure
plot(t,xmean)
xlabel('time(s)')
legend('x1','x2')

% kalman

particleKalman =  cell(length(t),1);
particleKalman{1} = x0';
time1 = tic;
for k = 2:length(t)
     particleKalman{k} = dynUpdateParticles( particleKalman{k-1}, fDYN,[t(k-1) t(k)]);
end
toc(time1)
xKmean = zeros(dim,length(t));
for k = 1:length(t)
    xKmean(:,k) = particleKalman{k};
end

figure
plot(t,squeeze(covP(1,1,:)))
xlabel('time(s)')
ylabel('cov')

afigure
plot(t,xmean(1,:),t,xmean(2,:),'r',t,xKmean(1,:),'--black',t,xKmean(2,:),'--red')
xlabel('time(s)')
legend('x_1 PF','x_2 PF','x_1 KF','x_2 KF')


%% Display particles animation
hf = figure;
save_to_file = 0;
saveResults = 1;
if (save_to_file && exist('saveResults','var') )
    FPS = 25;  
    pause(1)
    str_title = ['Probability Density Function Evolution'];
    writerObj = VideoWriter('pdf_gaussian_.avi');
    writerObj.FrameRate = FPS;
    myVideo.Quality = 100;
    set(hf,'Visible','on');
    open(writerObj);
    set(gcf,'Renderer','OpenGL'); %to save to file
    pause(2)
end
hold on
ht_kf = scatter(particle{1}(1,:),particle{1}(2,:),1);
xlabel('\theta')
ylabel('\dot{theta}')
axis([-pi pi -5 5])
axis equal
grid on
for k = 2:1:length(t)
    set ( ht_kf, 'XData',particle{k}(1,:),'YData',particle{k}(2,:) );
    pause(1.0/10);
    axis([-pi pi -5 5])
    drawnow
    if (save_to_file  && exist('saveResults','var') )
        M = getframe(gcf);
        writeVideo(writerObj, M);
    end        
end
if (save_to_file   && exist('saveResults','var') )
    close(writerObj);
end