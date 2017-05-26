%% Test Finite Difference PDF Evolution 2 dimensions

clear all
% Initial distribution

dim = 2;
dx = 0.1;
x0 = [7,13];%7:6+dim;
sigma02 = 2.5;
xvector = [-5:dx:20]';
%RotA = orth(randn(dim,dim));
sigma02Matrix = [2.4 0.4; 0.4 sigma02];%diag(sigma02+x0/dim/2);
n = length(xvector);
tic
[xgrid,ygrid] = meshgrid(xvector,xvector);
xyvector = [reshape(xgrid,n*n,1),reshape(ygrid,n*n,1)];
%p0 = reshape(mvnpdf(xyvector,x0,sigma02Matrix),n,n);

p0 = reshape(mvnpdf(xyvector,x0,sigma02Matrix),n,n)' + ...
     reshape(mvnpdf(xyvector,x0+[3;-4]',sigma02Matrix/0.2),n,n)' + ...
     reshape(mvnpdf(xyvector,x0+[-1;-6]',sigma02Matrix/1.4),n,n)';
 p0 = p0/trapz(xvector,trapz(xvector,p0));
toc
% tic
% for i=1:n
%     for j=1:n
%         p0(i,j) = mvnpdf([xvector(i) xvector(j)],x0,sigma02Matrix) ;
%     end
% end
% toc
%p0 = normpdf(x,x0,sqrt(sigma02));

if dim == 2
    hh = pcolor(xvector,xvector,p0);
    colorbar
    set(hh, 'EdgeColor', 'none');
    xlabel('x')
    zlabel('p(x)')
    grid on
end


% Simulation Parameters
dt = 0.01;
t = 0:dt:2;
px = zeros(n,n,length(t));
xkalman = zeros(dim,length(t));
covkalman = zeros(dim,dim,length(t));
expec = zeros(dim,length(t));
expec2 = zeros(dim,length(t));
cov = zeros(dim,dim,length(t));
q = diag([0.4,0.6]); %0.25; [0.5 0.1; 0.1 0.8]
aspeed = [0;0]; %[2;-2]; %[0,3]';

xkalman(:,1) = x0;
covKalman(:,:,1) = sigma02Matrix;
% Equation xdot = a*x + noise
expec(:,1) = x0;
cov(:,:,1) = sigma02Matrix;
px(:,:,1) = p0;
tic
for k = 2:length(t)
    
    for i=2:n-1
        for j=2:n-1
            fx = aspeed(1)*(px(i,j+1,k-1)-px(i,j-1,k-1))/dx/2;
            fy = aspeed(2)*(px(i+1,j,k-1)-px(i-1,j,k-1))/dx/2;
            gx(1,1) = q(2,2)/2*( px(i+1,j,k-1)-2*px(i,j,k-1)+px(i-1,j,k-1))/dx^2;
            gx(1,2) = q(1,2)/2*( px(i+1,j+1,k-1)-px(i-1,j,k-1)-px(i,j-1,k-1)+px(i-1,j-1,k-1))/4/dx^2;
            gx(2,1) = q(2,1)/2*( px(i+1,j+1,k-1)-px(i-1,j,k-1)-px(i,j-1,k-1)+px(i-1,j-1,k-1))/4/dx^2;
            gx(2,2) = q(1,1)/2*( px(i,j+1,k-1)-2*px(i,j,k-1)+ px(i,j-1,k-1))/dx^2;
             
            px(i,j,k) = px(i,j,k-1) - dt*(fx+fy) + dt*sum(sum(gx));
        end
    end
    %get expected values

    expec(:,k)  = [ trapz( xvector, xvector'.*trapz(xvector, px(:,:,k),1 ) )
                  trapz( xvector, xvector.*trapz(xvector, px(:,:,k),2 ) )];
              
    xgridDiff = xgrid - expec(1,k);
    ygridDiff = ygrid - expec(2,k);
   
    cov(1,1,k) = trapz( xvector, trapz(xvector,xgridDiff.*xgridDiff.*px(:,:,k) ));
    cov(1,2,k) = trapz( xvector, trapz(xvector,xgridDiff.*ygridDiff.*px(:,:,k) ));
    cov(2,1,k) = cov(1,2,k);
    cov(2,2,k) = trapz( xvector, trapz(xvector,ygridDiff.*ygridDiff.*px(:,:,k) ));

    xkalman(:,k) = xkalman(:,k-1) + aspeed*dt;
    covKalman(:,:,k) = covKalman(:,:,k-1) + q*dt;
end
toc
if dim == 2
    hh = pcolor(xvector,xvector,px(:,:,end));
    colorbar
    set(hh, 'EdgeColor', 'none');
    xlabel('x')
    zlabel('p(x)')
    grid on
end

figure
plot(t,xkalman,t,expec,'.')
xlabel('time')
zlabel('position')
legend('KalmaX','KalmanY','FPE X','FPE Y' )

figure
plot2d( xvector, reshape(mvnpdf(xyvector, xkalman(:,k)',covKalman(:,:,end)),n,n)-px(:,:,end) )
title('Difference with kalman')

figure
plot2d( xvector, reshape(mvnpdf(xyvector,expec(:,end)',cov(:,:,end)),n,n) - px(:,:,end) )


figure
plot(t,squeeze(covKalman(1,1,:)),'b',t,squeeze(cov(1,1,:)),'b--',t,squeeze(covKalman(2,2,:)),'r',t,squeeze(cov(2,2,:)),'r--')
xlabel('time')
zlabel('position')
legend('KalmaX','FPE X','KalmanY','FPE Y')



% Measure
zmes = [12,17]';
Rmes = diag([1 1]); %0.5*[3,1;1,2]*[3,1;1,2];
pz = reshape(mvnpdf(xyvector,zmes',Rmes),n,n);
pzpos = pz.*px(:,:,end);
pzposNorm = trapz(xvector,trapz(xvector,pzpos));
pzpos = pzpos/pzposNorm;

% subplot(2,2,1)
% plot2d( xvector,px(:,:,end) )
% subplot(2,2,2)
% plot2d( xvector,pz )
% subplot(2,2,3)
% plot2d( xvector,pzpos )
% subplot(2,2,4)
% hold on
% hh = surf(xvector,xvector,px(:,:,end),'b');
% set(hh, 'EdgeColor', 'none');
% hh = surf(xvector,xvector,pz);
% set(hh, 'EdgeColor', 'none');
% hh = surf(xvector,xvector,pzpos);
% set(hh, 'EdgeColor', 'none');


xpos = [ trapz( xvector, xvector'.*trapz(xvector, pzpos,1 ) )
         trapz( xvector,  xvector.*trapz(xvector, pzpos,2 ) )];
xgridDiff = xgrid - xpos(1);
ygridDiff = ygrid - xpos(2);
covpos= [ trapz( xvector, trapz(xvector,xgridDiff.*xgridDiff.*pzpos ))  trapz( xvector, trapz(xvector,xgridDiff.*ygridDiff.*pzpos ))
               trapz( xvector, trapz(xvector,xgridDiff.*ygridDiff.*pzpos ))  trapz( xvector, trapz(xvector,ygridDiff.*ygridDiff.*pzpos ))];

HK = eye(dim);

% Kalman Measurement
KalmanG = covKalman(:,:,end)*HK'*inv(HK*covKalman(:,:,end)*HK'+Rmes);
xposKalman = xkalman(:,end) + KalmanG*(zmes - HK*xkalman(:,end));
covKalmanPos = (eye(dim)-KalmanG*HK)*covKalman(:,:,end);

fprintf(['Measure update: Mean values ',num2str(xposKalman,4),' Kalman, ',num2str(xpos, 4) , ' Tensor-FPE \n'])
fprintf(['Measure update: Mean cov ',num2str(covKalmanPos,4),' Kalman, ',num2str(covpos, 4) , ' Tensor-FPE \n'])

%%  animated figure
figure
ht_pdf = surf(xvector,xvector,p0);
h = colorbar;
caxis([0 max(max(p0))])
set(ht_pdf, 'EdgeColor', 'none');
xlabel('x')
ylabel('y')
grid on
axis ([0,20,0,20])
grid on
for i = 2:1:length(t)
     set ( ht_pdf, 'ZData', px(:,:,i) );
     drawnow
     pause(1.0/100.0);
end


%% animated figure difference
figure
ht_pdf = pcolor(xvector,xvector,reshape(mvnpdf(xyvector, xkalman(:,k)',covKalman(:,:,end)),n,n)- px(:,:,1) );
h = colorbar;
%set(h, 'ylim', [0 0.25])
caxis([0 0.25])
set(ht_pdf, 'EdgeColor', 'none');
xlabel('x')
ylabel('y')
grid on
axis ([0,25,0,25])
grid on
for k = 2:10:length(t)
     set ( ht_pdf, 'CData', reshape(mvnpdf(xyvector, xkalman(:,k)',covKalman(:,:,end)),n,n)- px(:,:,k)  );
     drawnow
     pause(1.0/10.0);
end