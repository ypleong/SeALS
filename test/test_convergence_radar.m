
%% Load data
load('gridVTOL.mat')
%% Analyze Measurement PDFs for non-linear functions



% Old
%zradar = [1.04,4.6996];
%r1 = [0,6]';
%r2 = [2,4]'; 
r1 = [-1,5]';
r2 = [3,3.5]'; 
xGT = [1,5];
zradar = [(xGT(1)-r1(1))^2+(xGT(2)-r1(2))^2 ;%+ randn(1)*sqrt(sigmaX2)
              (xGT(1)-r2(1))^2+(xGT(2)-r2(2))^2] ;%+ randn(1)*sqrt(sigmaX2)];

gridXY = {gridT{1},gridT{3}}; % extract only x and y grids
for i=1:length(gridXY{1})
    for j=1:length(gridXY{2})
        pzx1(i,j) = exp( -(zradar(1)-hRadar([gridXY{1}(i),0,gridXY{2}(j),0,0,0],r1))^2/sigmaX2/25/2);
        pzx2(i,j) = exp( -(zradar(2)-hRadar([gridXY{1}(i),0,gridXY{2}(j),0,0,0],r2))^2/sigmaX2/25/2);
        pzx(i,j) = pzx1(i,j)*pzx2(i,j) ;
        pzxPlus(i,j) = pzx1(i,j)+pzx2(i,j);
    end
end

%% Show basis
figure
set(gca,'fontsize',20)
hold on
hPcolor1 = pcolor(gridXY{1},gridXY{2},pzx1');
xlabel('x')
ylabel('y')
scatter(xGT(1),xGT(2),'r','f')
scatter(r1(1),r1(2),'black','f')
colorbar
set(hPcolor1, 'EdgeColor','none')
legend('PDF','Target','Radar')

figure
set(gca,'fontsize',20)
hold on
hPcolor2 = pcolor(gridXY{1},gridXY{2},pzx2');
xlabel('x')
ylabel('y')
scatter(xGT(1),xGT(2),'r','f')
scatter(r2(1),r2(2),'black','f')
set(hPcolor2, 'EdgeColor','none')
legend('PDF','Target','Radar')

figure
hold on
set(gca,'fontsize',20)
hPcolorPlus = pcolor(gridXY{1},gridXY{2},pzxPlus');
xlabel('x')
ylabel('y')
scatter(xGT(1),xGT(2),'r','f')
scatter(r1(1),r1(2),'black','f')
scatter(r2(1),r2(2),'black','f')
set(hPcolorPlus, 'EdgeColor','none')
legend('PDF','Target','Radar')

figure
set(gca,'fontsize',20)
hold on
hPcolorboth = pcolor(gridXY{1},gridXY{2},pzx');
xlabel('x')
ylabel('y')
scatter(xGT(1),xGT(2),'r')
scatter(r1(1),r1(2),'black','f')
scatter(r2(1),r2(2),'black','f')
set(hPcolorboth, 'EdgeColor','none')
legend('PDF','Target','Radar')
colorbar


%% Show basis functions with increasing order
rankV = 1:15; 
err = zeros(length(rankV),1);
for i = 1:length(rankV)
    [pKxy,Uo,out(i)] = cp_als(tensor(pzx),rankV(i),'init','nvecs');
    err(i) = norm(double(pKxy)-pzx);
    %handles = plot2DslicesAroundPoint(pKxy, [0,0], gridXY);
    %ylim(handles{1,1},[-3 3])
    %ylim(handles{2,2},[-3 3])
    %set(handles{1,1},'Color','black')
    %colormap(handles{1,1},'gray')
    %saveas(gcf,['convergence_rank_',num2str(rankV(i)),'.png'])
    comp(i) = (length(gridXY{1})+length(gridXY{1}))*rankV(i)/length(gridXY{1})/length(gridXY{2});
    %clf
end
afigure
plot(rankV,err)
xlabel('Tensor Decomposition Rank')
ylabel('Error Norm')

afigure
plot(rankV,comp)

%% Show the influence of random init
rankV = 20*ones(200,1);
err = zeros(length(rankV),1);
figure
grid on
for j=1:rankV(1)
    hf(j,1) = subplot(2,rankV(1),j);
    xlabel(['X Basis Function ',num2str(j)])
    ylabel('Value')
    hold( hf(j,1),'on')
    grid( hf(j,1),'on')
    hf(j,2) = subplot(2,rankV(1),j+rankV(1));
    hold( hf(j,2),'on')
    ylabel('Value')
    xlabel(['Y Basis Function ',num2str(j)])
    grid( hf(j,2),'on')
end
for i = 1:length(rankV)
    pKxy = cp_als(tensor(pzx),rankV(i));%,'init','nvecs');
    for j=1:rankV(1)
        plot(hf(j,1),gridXY{1},pKxy.lambda(j)*pKxy.U{1}(:,j),'black');
        plot(hf(j,2),gridXY{2},pKxy.lambda(j)*pKxy.U{2}(:,j),'black');
    end
    %handles =plot2DslicesAroundPoint(pKxy, [0,0], gridXY);
    %ylim(handles{1,1},[-8 8])
    %ylim(handles{2,2},[-8 8]) 
    %pause(1)
end
[mean(err),mean(errNVECS)]
figure
hist(err)
xlabel('Tensor Decomposition Rank')
ylabel('Error Norm')


%% Show the extended mesurement pdf to 6D for the VTOL case
pKxy = cp_als(tensor(pzx),5);
dim = 6;
dim_z = [1,3];
Ucell = cell(dim,1);
%surf(gridXY{1},gridXY{2},pzx,'EdgeColor','none')
for i=1:dim
    if all(i ~= dim_z)
        Ucell{i} = ones(length(gridT{i}),5);
    else
        
        Ucell{i} = pKxy.U{find(dim_z==i)};
    end
end
pKxyFull = ktensor(pKxy.lambda,Ucell);
for i=1:dim
    expec(i)  = intTens(pKxyFull, [], gridT, we(i,:));
end

plot2DslicesAroundPoint(pKxyFull, expec, gridT);

%% build a sparse tensor from the points
k = 1;
spValue = 0;
spSubs = [1,0];
for i=1:10:length(gridXY{1})
    for j=1:10:length(gridXY{2})
        pzx1(i,j) = exp( -(zradar(1)-hRadar([gridXY{1}(i),0,gridXY{2}(j),0,0,0],r1))^2/sigmaX2/25/2);
        pzx2(i,j) = exp( -(zradar(2)-hRadar([gridXY{1}(i),0,gridXY{2}(j),0,0,0],r2))^2/sigmaX2/25/2);
        pzx(i,j) = pzx1(i,j)*pzx2(i,j) ;
        spValue(k) = pzx(i,j);
        spSubs(k,:) = [i,j];
        k = k+1;
    end
end
dim = 2;
% Build ktensor sparse
for i=1:dim
   UcellSP{i}=zeros(length(gridXY{i}),length(spValue));
end
for i=length(spValue)
    for j=1:dim
        UcellSP{j}(spSubs(i,j),i)=1;
    end
end
sparseKtensorPzk = ktensor(spValue',UcellSP);

plot2DslicesAroundPoint(sparseKtensorPzk, [0,0], gridXY);

[spPzkComp, err_op, iter_op, enrich_op, t_step_op, cond_op, noreduce] = als2(sparseKtensorPzk,1e-3);

spPzk = sptensor(spSubs,spValue',[length(gridXY{1}) length(gridXY{2})]);

kspPxk = cp_als(full(spPzk),5,'init','nvecs');

plot2DslicesAroundPoint(kspPxk, [0,0], gridXY);

figure
set(gca,'fontsize',20)
hold on
hPcolorboth = pcolor(gridXY{1},gridXY{2},double(kspPxk)');
xlabel('x')
ylabel('y')
scatter(xGT(1),xGT(2),'r')
scatter(r1(1),r1(2),'black','f')
scatter(r2(1),r2(2),'black','f')
set(hPcolorboth, 'EdgeColor','none')
legend('PDF','Target','Radar')
colorbar
