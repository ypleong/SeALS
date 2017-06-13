%% Test Finite Difference PDF Evolution
clear all
% Initial distribution

dx = 0.05;
x0 = 2;
sigma02 = 0.5;
x = [-5:dx:25]';
p0 = normpdf(x,x0,sqrt(sigma02));
%p0 = p0 + normpdf(x,x0+3,sqrt(sigma02/0.8));
p0 = p0/sum(p0)/dx;

plot(x,p0)
xlabel('x')
ylabel('p(x)')
grid on


% Simulation Parameters
dt = 0.01;
t = 0:dt:15;
px = zeros(length(x),length(t));
xkalman = zeros(length(t),1);
covKalman = zeros(length(t),1);
exp = zeros(length(t),1);
exp2 = zeros(length(t),1);
cov = zeros(length(t),1);
q = 0.25; %0.25;
aspeed = 1;

xkalman(1) = x0;
covKalman(1) = sigma02;
% Equation xdot = a*x + noise
exp(1) = x0;
cov(1) = sigma02;
px(:,1) = p0;
for i = 2:length(t)
    for j=2:length(x)-1
        fx = aspeed*(px(j+1,i-1)-px(j-1,i-1))/dx/2;
        gx = q/2*(px(j+1,i-1)-2*px(j,i-1)+px(j-1,i-1))/dx^2;
        px(j,i) = px(j,i-1) - dt*fx + dt*gx;
    end
    %get expected values
    exp(i)  = px(:,i)'*x*dx;
    exp2(i) = px(:,i)'*(x.^2)*dx;
    cov(i) = exp2(i) - exp(i).^2;
    cov(i) = sum(px(:,i).*(x - exp(i)).*(x - exp(i)))*dx;
    
    %trapz(x,(x-exp(i)).*(x-exp(i)).*px(:,i) );
    
    xkalman(i) = xkalman(i-1) + aspeed*dt;
    xkalman(i);
    covKalman(i) = covKalman(i-1) + q*dt;
    px(:,i) = px(:,i)/sum(px(:,i))/dx;
end

figure
plot(x,normpdf(x,exp(end),sqrt(cov(end))),x,px(:,end),x,normpdf(x,xkalman(end),sqrt(covKalman(end))))
legend('recoverNormal','px','Kalman')


%% Error
RMS_mean = sqrt(mean((xkalman-exp).^2));
RMS_cov  = sqrt(mean((covKalman-cov).^2));

figure
plot(t,covKalman,t,cov)


%% Measurement
zmes = 11;
Rmes = 2;
pz = normpdf(x,zmes,sqrt(Rmes));
pzpos = zeros(length(x),1);
for j=1:length(x)
    pzpos(j) = px(j,end)*pz(j);
end
pzpos = pzpos/sum(pzpos)/dx;
plot(x,pzpos)
title('posterior')
xpos = pzpos'*x*dx;
covpos = sum(pzpos.*(x - xpos).*(x - xpos))*dx;
% Kalman Measurement
KalmanG = covKalman(end)/(covKalman(end)+Rmes);
xposKalman = xkalman(end) + KalmanG*(zmes - xkalman(end));
covKalmanPos = covKalman(end)*(1-KalmanG);

plot(x,px(:,end),x,pz,x,pzpos)
legend('prior','measure','posterior')

plot(sum(px)*dx)

hf = figure;
hold on
ht_pdf = plot(x,p0);
ht_gaussian = plot(x,normpdf(x,exp(1),sqrt(cov(1))));
%legend('Non Linear', 'Gaussian Approximation')
axis ([-3,15,-0.1,1])
grid on
save_to_file = 1;
if (save_to_file  )
    FPS = 25;  
    str_title = ['Probability Density Function Evolution. Simple case.'];
    writerObj = VideoWriter('pdf_gaussian.avi');
    writerObj.FrameRate = FPS;
    myVideo.Quality = 100;
    set(hf,'Visible','on');
    open(writerObj);
    set(gcf,'Renderer','OpenGL'); %to save to file
end
for i = 2:5:length(t)
     set(ht_pdf, 'YData', px(:,i) );
     set(ht_gaussian, 'YData', normpdf(x,exp(i),sqrt(cov(i))) );
     drawnow
     pause(1.0/50.0);
     if (save_to_file )
         M = getframe(gcf); %hardcopy(hf, '-dzbuffer', '-r0');
         writeVideo(writerObj, M);
     end
end
if (save_to_file  )
    close(writerObj);
end

[cov(end),covKalman(end)]


figure
plot(t,exp,t,xkalman)
xlabel('Time(s)')
ylabel('Mean Position(m)')

figure
plot(t,cov,t,covKalman)
xlabel('Time(s)')
ylabel('Mean Cov')
title('check that both lines are similar')

figure
h=pcolor(x,t,px') %,'EdgeColor','none')
set(h, 'EdgeColor', 'none');
xlabel('X(m)')
ylabel('Time(s)')
zlabel('Pdf(x,t)')
title('Pdf Evolution')
colorbar


