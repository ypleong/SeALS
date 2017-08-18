clear all
close all


% Parameters
g = 9.81;
epsilon = 0.01;


% Simulation Parameters
dt = 0.01;
t = 0:dt:15;
x = zeros(6,length(t)); % State: x,vx,y,vy,theta,vtheta
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

Q = 0.01*diag([0.1^2, 0.1^2, 0.01^2]);

GQGp = G*Q*G';
sigmaX2 = 0.01^2;
sigmaY2 = 0.01^2;
R = diag([sigmaX2,sigmaY2]);
zmes = zeros(2,length(t));
HK = [1 0 0 0 0 0
      0 0 1 0 0 0];

xref = 2; %+ cos(t(i)*2*pi/20);
yref = 6; %+ sin(t(i)*2*pi/20);
  
%% Time Evolution 
for i = 2:length(t)
   
    uControl(i) =  g - 0.3*(xhat(3,i-1)-yref)  - 0.9*xhat(4,i-1);
    tauControl(i) = -6*xhat(5,i-1) - 6*xhat(6,i-1) + 0.2*(xhat(1,i-1)-xref) + 1.0*xhat(2,i-1);
    % Propagate 
    x(:,i) =  x(:,i-1)  + dt*fVTOL( x(:,i-1), uControl(i), tauControl(i), epsilon, g );
    % Kalman
    F = [0 1 0 0 0 0
         0 0 0 0 (-tauControl(i)*cos(xhat(5,i-1))-epsilon*tauControl(i)*sin(xhat(5,i-1))) 0
         0 0 0 1 0 0
         0 0 0 0 (-tauControl(i)*sin(xhat(5,i-1))+epsilon*tauControl(i)*cos(xhat(5,i-1))) 0
         0 0 0 0 0 1
         0 0 0 0 0 0];
    xhat(:,i) =  xhat(:,i-1)  + dt*fVTOL( xhat(:,i-1), uControl(i), tauControl(i), epsilon, g );
    Phat(:,:,i) = F*Phat(:,:,i-1)*F' + dt*GQGp;
    
    % assume measure with noise R
    z = [x(1,i) + randn(1)*sqrt(sigmaX2) 
         x(3,i) + randn(1)*sqrt(sigmaY2)];
    zmes(:,i) = z;
    SG = HK*Phat(:,:,i)*HK' + R;
    KG = Phat(:,:,i)*HK'/(SG);
    xhat(:,i) = xhat(:,i) + KG*(z - HK*xhat(:,i));
    Phat(:,:,i) = (eye(6) - KG*HK)*Phat(:,:,i);
    
    % Propagate 
    
end


%% Plot
%afigure;

afigure
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
plot(t,uControl)
grid on
xlabel('Time(s)')
ylabel('u input')
subplot(2,4,8)
plot(t,tauControl)
grid on
xlabel('Time(s)')
ylabel('\tau input')

afigure
hold on
plot( x(1,:), x(3,:) ,xhat(1,:), xhat(3,:) )
scatter ( zmes(1,:), zmes(2,:))
legend('truth','kalman')

afigure
plot( t, sqrt( (x(1,:)-xhat(1,:)).^2 + (x(3,:)- xhat(3,:)).^2) )
grid on
xlabel('Time(s)')
ylabel('Estimation Error(m)')

hf = figure
hold on
axis ([0,4,3,8])
scatter(xref,yref)
xlabel('X')
ylabel('Y')
title('VTOL Simulation')
ht_pdf = quiver(x(1,1),x(3,1),-sin(x(5,1)),cos(x(5,1)));
set(ht_pdf,'linewidth',2);
ht_scatter = plot(x(1,1),x(3,1));
set(ht_scatter,'linewidth',2);
grid on


save_to_file = 0;
if (save_to_file  )
    FPS = 25;  
    writerObj = VideoWriter('pdf_gaussian_.avi');
    writerObj.FrameRate = FPS;
    set(hf,'Visible','on');
    open(writerObj);
    set(gcf,'Renderer','OpenGL'); %to save to file
    pause(5)
end
for i = 2:5:length(t)
     set ( ht_pdf, 'XData', x(1,i), ...
                   'YData', x(3,i), ...
                   'UData', -sin(x(5,i)), ... 
                   'VData', cos(x(5,i)) );
     set ( ht_scatter, 'XData', x(1,1:i), ...
                   'YData', x(3,1:i)  );
     drawnow
     pause(1.0/100.0);
    if (save_to_file )
        M = getframe(gcf); 
        writeVideo(writerObj, M);
    end  
end
if (save_to_file  )
    close(writerObj);
end
