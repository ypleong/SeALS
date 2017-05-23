clear all
close all


% Parameters
g = 9.81;
epsilon = 0.01;


% Simulation Parameters
dt = 0.01;
t = 0:dt:20;
x = zeros(6,length(t));
%x(:,1) = [1, 0.1, 5, 0.2, 0.1, -0.01]'; % some movement
x(:,1) = [1, 0.0, 5, 0.0, 0.00, -0.00]'; % static


% Kalman Estimation
xhat = zeros(6,length(t));
xhat(:,1) = [1.1, 0.0, 5.05, 0.0, 0.00, -0.00]';
Phat = zeros(6,6,length(t));
Phat(:,:,1) = diag([3^2,0.2^2,3^2,0.2^2,0.1^2,0.01^2]);

G = [0 0 0
 1 0 0
 0 0 0
 0 1 0
 0 0 0
 0 0 1];


Q = 1*diag([0.1^2, 0.1^2, 0.01^2]);

GQGp = G*Q*G';
sigmaX2 = 0.3^2;
sigmaY2 = 0.3^2;
R = diag([sigmaX2,sigmaY2]);
zmes = zeros(2,length(t));

%radars
r1 = [0,6]';
r2 = [2,4]'; 

  
for i = 2:length(t)
    xref = 2; %+ cos(t(i)*2*pi/20);
    yref = 6; %+ sin(t(i)*2*pi/20);
    u = g - 0.2*(x(3,i-1)-yref)  - 0.7*x(4,i-1);
    uu(i) = u;
    tau = -6*x(5,i-1) - 6*x(6,i-1) + 0.2*(x(1,i-1)-xref) + 0.9*x(2,i-1);
    tautau(i) = tau;
    % Propagate 
    x(:,i) =  x(:,i-1)  + dt*fVTOL( x(:,i-1), u, tau, epsilon, g );
    % Kalman
    F = [0 1 0 0 0 0
         0 0 0 0 (-u*cos(x(5,i-1))-epsilon*tau*sin(x(5,i-1))) 0
         0 0 0 1 0 0
         0 0 0 0 (-u*sin(x(5,i-1))+epsilon*tau*cos(x(5,i-1))) 0
         0 0 0 0 0 1
         0 0 0 0 0 0];
    xhat(:,i) =  xhat(:,i-1)  + dt*fVTOL( xhat(:,i-1), u, tau, epsilon, g );
    Phat(:,:,i) = F*Phat(:,:,i-1)*F' + dt*GQGp;
    
    % assume measure with noise R
    z = [x(1,i) + randn(1)*sqrt(sigmaX2) 
         x(3,i) + randn(1)*sqrt(sigmaY2)];
     
    zradar = [(x(1,i)-r1(1))^2+(x(3,i)-r1(2))^2 + randn(1)*sqrt(sigmaX2)
              (x(1,i)-r2(1))^2+(x(3,i)-r2(2))^2 + randn(1)*sqrt(sigmaY2)];
    hradar = [(xhat(1,i)-r1(1))^2+(xhat(3,i)-r1(2))^2
              (xhat(1,i)-r2(1))^2+(xhat(3,i)-r2(2))^2];
    HK = [2*(x(1,i)-r1(1)) 0 2*(x(3,i)-r1(2)) 0 0 0
          2*(x(1,i)-r2(1)) 0 2*(x(3,i)-r2(2)) 0 0 0];
    zmes(:,i) = zradar;
    SG = HK*Phat(:,:,i)*HK' + R;
    KG = Phat(:,:,i)*HK'/(SG);
    xhat(:,i) = xhat(:,i) + KG*(zradar - hradar);
    Phat(:,:,i) = (eye(6) - KG*HK)*Phat(:,:,i);
end


%% Plot

figure
subplot(2,4,1)
plot(t,x(1,:))
grid on
xlabel('Time(s)')
ylabel('X')
subplot(2,4,2)
plot(t,x(2,:))
grid on
xlabel('Time(s)')
ylabel('VX')
subplot(2,4,3)
plot(t,x(3,:))
grid on
xlabel('Time(s)')
ylabel('Y')
subplot(2,4,4)
plot(t,x(4,:))
grid on
xlabel('Time(s)')
ylabel('VY')
subplot(2,4,5)
plot(t,x(5,:))
grid on
xlabel('Time(s)')
ylabel('\theta')
subplot(2,4,6)
plot(t,x(6,:))
grid on
xlabel('Time(s)')
ylabel('Theta Dot')
subplot(2,4,7)
plot(t,uu)
grid on
xlabel('Time(s)')
ylabel('u input')
subplot(2,4,8)
plot(t,tautau)
grid on
xlabel('Time(s)')
ylabel('\tau input')

figure
hold on
plot( x(1,:), x(3,:) ,xhat(1,:), xhat(3,:) )
scatter (r1(1),r1(2))
scatter (r2(1),r2(2))
%scatter ( zmes(1,:), zmes(2,:))
legend('truth','kalman','radar1','radar2')

figure
plot( t, sqrt( (x(1,:)-xhat(1,:)).^2 + (x(3,:)- xhat(3,:)).^2) )
grid on
xlabel('Time(s)')
ylabel('Estimation Error(m)')

detP = zeros(length(t),1);
for i=1:length(t)
   detP(i) = det(Phat(:,:,i));
end

figure
plot( t, detP )
grid on
xlabel('Time(s)')
ylabel('Covariance Determinant')


% figure
% hold on
% axis ([0,4,3,8])
% scatter(xref,yref)
% ht_pdf = quiver(x(1,1),x(3,1),-sin(x(5,1)),cos(x(5,1)));
% ht_scatter = plot(x(1,1),x(3,1));
% grid on
% for i = 2:5:length(t)
%      set ( ht_pdf, 'XData', x(1,i), ...
%                    'YData', x(3,i), ...
%                    'UData', -sin(x(5,i)), ... 
%                    'VData', cos(x(5,i)) );
%      set ( ht_scatter, 'XData', x(1,1:i), ...
%                    'YData', x(3,1:i)  );
%      drawnow
%      pause(1.0/100.0);
% end

% close all
% axis tight
% set(gca,'nextplot','replacechildren')
% set(gcf,'Renderer', 'zbuffer')
% f = getframe;
% [im,map] = rgb2ind(f.cdata,256,'nodither');
% im(1,1,1,length(t)) = 0;
% for i = 2:5:length(t)
%     quiver(x(1,i),x(3,i),-sin(x(5,i)),cos(x(5,i)));
%     axis ([0,4,3,8])
%     frame = getframe;
%     im(:,:,1,i) = rgb2ind(f.cdata,map,'nodither');
% end
%imwrite(im,map,'DancingPeaks.gif','DelayTime',0,'LoopCount',inf) %g443800