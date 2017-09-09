clear all
close all

dim = 6;
x = sym('x',[dim,1]); %do not change

% Parameters
g = 9.81;
epsilon = 0.01;

% Simulation Parameters
dt = 0.001;
t = 0:dt:15;
xGT = zeros(dim,length(t)); % State: x,vx,y,vy,theta,vtheta
xGT(:,1) = [1, 0.0, 5, 0.0, 0.00, -0.00]'; % static

xref = 2; %+ cos(t(i)*2*pi/20);
yref = 6; %+ sin(t(i)*2*pi/20);

%% Kalman Estimation Preparation
xhat = zeros(dim,length(t));
xhat(:,1) = [1.1, 0.0, 5.05, 0.0, 0.00, -0.00]';
Phat = zeros(dim,dim,length(t));
Phat(:,:,1) = diag([0.5^2,0.1^2,0.5^2,0.1^2,0.01^2,0.005^2]);

G = [0 0 0
 1 0 0
 0 0 0
 0 1 0
 0 0 0
 0 0 1];

Q = 1.0*diag([0.1^2, 0.1^2, 0.01^2]);

GQGp = G*Q*G';
sigmaX2 = 0.04^2;
sigmaY2 = 0.04^2;
R = diag([sigmaX2,sigmaY2]);
zmes = zeros(2,length(t));
HK = [1 0 0 0 0 0
      0 0 1 0 0 0];

%% Tensor Estimation Preparation

% Initial distribution
n = [81 51 81 51 81 51];
bdim = [-1 3; -0.5 0.5;3 7; -0.6 0.6;-0.08 0.08; -0.06 0.06]; 
diagSigma = diag(Phat(:,:,1));
qdiag = diag(GQGp)+ones(dim,1)*0.01^2;
q =  diag( qdiag );

% Tensor parameters

for i=1:dim
   bcon{i} = {'d',0,0};
end
bcon{5} = {'p'}; % angle

% Grid Configuration
fprintf ('The grid size is %d in dimension %d \n', [n',(1:dim)']')
for i=1:dim
    dx(i) = (bdim(i,2)-bdim(i,1))/(n(i)-1);
    gridT{i} =  (bdim(i,1):dx(i):bdim(i,2))';
end

sigma02Matrix = diag(diagSigma); % [0.8 0.2; 0.2 sigma02];
for i=1:dim
   p0vector{i} =  normpdf(gridT{i},xhat(i,1),sqrt(diagSigma(i))); %*0.5 + normpdf(gridT{i},2*x0(i),sqrt(diagSigma(i)))*0.5;
end
pk{1} = ktensor(p0vector);
plot2DslicesAroundPoint( pk{1}, xhat(:,1), gridT);
expec = zeros(dim,length(t));
covT = zeros(dim,dim,length(t));

covT(:,:,1) = sigma02Matrix;
expec(:,1) = xhat(:,1);

bsca = []; %no manual scaling
region = [];
regval = 1;
regsca = [];
sca_ver = 1;

tol_err_op = 1e-5;
tol_err = 1e-8;
als_options = [];
als_variant = []; 

for i=1:dim
    acc(:,i) = [2,2]';
end
[D,D2,fd1,fd2] = makediffop(gridT,n,dx,acc,bcon,region);

% Create weights for integration
for i=1:dim
    for j = 1:dim
        if i == j
            weMean{i,j} = gridT{j};
        else
            weMean{i,j} = ones(n(j),1);
        end         
    end
end
for i=1:dim
    weOnes{i} = ones(n(i),1);
end
weCov = cell(dim,dim,dim);
for i=1:dim
    for j = 1:dim
        weCov(i,j,:) = weOnes;
        if i == j
            weCov{i,j,i} = gridT{j}.*gridT{j};
        else
            weCov{i,j,i} = gridT{i};
            weCov{i,j,j} = gridT{j};
        end         
    end
end



use_tensor = 1;
%% Time Evolution 
time1 = tic;
for k = 2:length(t)
   
    % Get Control
    uControl(k) =  g - 0.3*(xhat(3,k-1)-yref)  - 0.9*xhat(4,k-1);
    tauControl(k) = -6*xhat(5,k-1) - 6*xhat(6,k-1) + 0.2*(xhat(1,k-1)-xref) + 1.0*xhat(2,k-1);
    % Propagate 
    xGT(:,k) =  xGT(:,k-1)  + dt*fVTOL( xGT(:,k-1), uControl(k), tauControl(k), epsilon, g );
    % Kalman
    F = [0 1 0 0 0 0
         0 0 0 0 (-tauControl(k)*cos(xhat(5,k-1))-epsilon*tauControl(k)*sin(xhat(5,k-1))) 0
         0 0 0 1 0 0
         0 0 0 0 (-tauControl(k)*sin(xhat(5,k-1))+epsilon*tauControl(k)*cos(xhat(5,k-1))) 0
         0 0 0 0 0 1
         0 0 0 0 0 0];
    xhat(:,k) =  xhat(:,k-1)  + dt*fVTOL( xhat(:,k-1), uControl(k), tauControl(k), epsilon, g );
    Phat(:,:,k) = F*Phat(:,:,k-1)*F' + dt*GQGp;
    
    
    % Tensor
    if (use_tensor)
%         for i=1:dim % rebuild the grid
%             bdim(i,:) = [xhat(i,k)-10*sqrt(Phat(i,i,k)) xhat(i,k)+10*sqrt(Phat(i,i,k))]
%             dx(i) = (bdim(i,2)-bdim(i,1))/(n(i)-1);
%             gridT{i} =  (bdim(i,1):dx(i):bdim(i,2))';
%         end
        
        op = createOP_VTOL(D, D2, gridT, tol_err_op, uControl(k), tauControl(k), g, epsilon, x, dt, q);
        pk{k} = SRMultV( op, pk{k-1});
        %length(pk{k}.lambda)
        if length(pk{k}.lambda) > 1
            [pk{k},~] = tenid(pk{k},tol_err_op,1,9,'frob',[],fnorm(pk{k}),0);

            pk{k} = fixsigns(arrange(pk{k}));
            [pk{k}, err_op, iter_op, enrich_op, t_step_op, cond_op, noreduce] = als2(pk{k},tol_err_op);

        end

        % Get Mean and Covariance values
        for i=1:dim
            expec(i,k)  = intTens(pk{k}, [], gridT, weMean(i,:));
        end
        for i=1:dim
            for j = 1:dim
                covT(i,j,k) = intTens(pk{k}, [], gridT, weCov(i,j,:)) - expec(i,k)*expec(j,k);
            end
        end
    end
    diagT(:,k) = diag(covT(:,:,k));
    if ( mod(k,10)==0 )
        
        % assume measure wkth nokse R
        z = [xGT(1,k) + randn(1)*sqrt(sigmaX2) 
             xGT(3,k) + randn(1)*sqrt(sigmaY2)];
        zmes(:,k) = z;
        SG = HK*Phat(:,:,k)*HK' + R;
        KG = Phat(:,:,k)*HK'/(SG);
        xhat(:,k) = xhat(:,k) + KG*(z - HK*xhat(:,k));
        Phat(:,:,k) = (eye(6) - KG*HK)*Phat(:,:,k);
        
        if (use_tensor)
            for i=1:dim
                for j=1:length(z)
                    if ( HK(j,i) == 1 )
                        pz{i} =  normpdf(gridT{i},z(j),sqrt(R(j,j))); %*0.5 + normpdf(gridT{i},2*x0(i),sqrt(diagSigma(i)))*0.5;
                    elseif (all(HK(:,i)==0) )
                        pz{i} = ones(n(i),1);
                    end
                end
            end

            zkcompressed = ktensor(pz);
            % Direct PDF Bayesian Measurement
            pk{k} = HadTensProd(pk{k},zkcompressed);
            pk{k} = pk{k} *(1/  intTens(pk{k}, [], gridT, weOnes));
            

            [ expec(:,k), covT(:,:,k) ] = meanCovTensor( pk{k}, gridT, weMean, weCov, weOnes );
%             for i=1:dim
%                 expec(i,k)  = intTens(pk{k}, [], gridT, we(i,:));
%             end
%             for i=1:dim
%                 for j = 1:dim
%                     covT(i,j,k) = intTens(pk{k}, [], gridT, weCov(i,j,:)) - expec(i,k)*expec(j,k);
%                 end
%             end
            
            
            [ pk{k}, gridT, dx] = fitTensorBoundaries( pk{k}, gridT, expec(:,k), covT(:,:,k), n, 6 );
            [D,D2,~,~] = makediffop(gridT,n,dx,acc,bcon,region);
            op = createOP_VTOL(D, D2, gridT, tol_err_op, uControl(k), tauControl(k), g, epsilon, x, dt, q);
            [ weMean, weCov, weOnes ] = createWeights( gridT, n );
            %figure
            %plot2DslicesAroundPoint( pk{k}, expec(:,k), gridT )
            
            
            
        end
        
    end
    
end
toc(time1)


%% Tensor plot

% handleSlices = plot2DslicesAroundPoint( pk{1}, expec(:,1), gridT);
% for k = 2:10:length(t)
%     plot2DslicesAroundPoint( pk{k}, expec(:,k), gridT, handleSlices)
%     pause(1.0/1);      
% end
%% Plot
%afigure;

afigure
subplot(2,4,1)
plot(t,xGT(1,:))
grid on
xlabel('Time(s)')
ylabel('X')
subplot(2,4,2)
plot(t,xGT(2,:))
grid on
xlabel('Time(s)')
ylabel('VX')
subplot(2,4,3)
plot(t,xGT(3,:))
grid on
xlabel('Time(s)')
ylabel('Y')
subplot(2,4,4)
plot(t,xGT(4,:))
grid on
xlabel('Time(s)')
ylabel('VY')
subplot(2,4,5)
plot(t,xGT(5,:))
grid on
xlabel('Time(s)')
ylabel('\theta')
subplot(2,4,6)
plot(t,xGT(6,:))
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
plot( xGT(1,:), xGT(3,:) ,xhat(1,:), xhat(3,:) )
scatter ( zmes(1,:), zmes(2,:),2)
legend('truth','kalman')
xlabel('X')
ylabel('Y')

afigure
hold on
plot( t, sqrt( (xGT(1,:)-xhat(1,:)).^2 + (xGT(3,:)- xhat(3,:)).^2) )
plot( t, sqrt( (xGT(1,:)-expec(1,:)).^2 + (xGT(3,:)- expec(3,:)).^2) )
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
ht_pdf = quiver(xGT(1,1),xGT(3,1),-sin(xGT(5,1)),cos(xGT(5,1)));
set(ht_pdf,'linewidth',2);
ht_scatter = plot(xGT(1,1),xGT(3,1));
set(ht_scatter,'linewidth',2);
grid on

for k=1:length(t)
   detPhat(k) = det(Phat(:,:,k)); 
end

save_to_file = 1;
if (save_to_file  )
    FPS = 25;  
    writerObj = VideoWriter('pdf_gaussian_VTOL_.avi');
    writerObj.FrameRate = FPS;
    set(hf,'Visible','on');
    open(writerObj);
    set(gcf,'Renderer','OpenGL'); %to save to file
    pause(5)
end
for i = 2:5:length(t)
     set ( ht_pdf, 'XData', xGT(1,i), ...
                   'YData', xGT(3,i), ...
                   'UData', -sin(xGT(5,i)), ... 
                   'VData', cos(xGT(5,i)) );
     set ( ht_scatter, 'XData', xGT(1,1:i), ...
                   'YData', xGT(3,1:i)  );
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



%%
k = 1;
for i=1:dim
    bdim(i,:) = [xhat(i,k)-10*sqrt(Phat(i,i,k)) xhat(i,k)+10*sqrt(Phat(i,i,k))]
    dx(i) = (bdim(i,2)-bdim(i,1))/(n(i)-1);
    gridT{i} =  (bdim(i,1):dx(i):bdim(i,2))';
end
handleSlices = plot2DslicesAroundPointGaussian(Phat(:,:,1), xhat(:,1), gridT);
for k = 2:30:length(t)
    for i=1:dim
        bdim(i,:) = [xhat(i,k)-4*sqrt(Phat(i,i,k)) xhat(i,k)+4*sqrt(Phat(i,i,k))]
        dx(i) = (bdim(i,2)-bdim(i,1))/(n(i)-1);
        gridT{i} =  (bdim(i,1):dx(i):bdim(i,2))';
    end
    plot2DslicesAroundPointGaussian(Phat(:,:,k), xhat(:,k), gridT, handleSlices);
    drawnow
    pause(1.0/100.0); 
    k/length(t)
end