

%%% Analyze Measurement PDFs for non-linear functions

%% Load data
load('gridVTOL.mat')

z = [1.04,4.6996];
r1 = [0,6]';
r2 = [2,4]'; 

gridXY = {gridT{1},gridT{3}}; % extract only x and y grids
for i=1:length(gridXY{1})
    for j=1:length(gridXY{2})
        pzx1(i,j) = exp( -(z(1)-hRadar([gridXY{1}(i),0,gridXY{2}(j),0,0,0],r1))^2/sigmaX2/25/2);
        pzx2(i,j) = exp( -(z(2)-hRadar([gridXY{1}(i),0,gridXY{2}(j),0,0,0],r2))^2/sigmaX2/25/2);
        pzx(i,j) = pzx1(i,j)*pzx2(i,j) ;
    end
end

%% Show basis functions with increasing order
rankV = 1:15; 
err = zeros(length(rankV),1);
for i = 1:length(rankV)
    pKxy = cp_als(tensor(pzx),rankV(i),'init','nvecs');
    err(i) = norm(double(pKxy)-pzx);
    handles =plot2DslicesAroundPoint(pKxy, [0,0], gridXY);
    ylim(handles{1,1},[-3 3])
    ylim(handles{2,2},[-3 3])
    %set(handles{1,1},'Color','black')
    %colormap(handles{1,1},'gray')
    saveas(gcf,['convergence_rank_',num2str(rankV(i)),'.png'])
    clf
end
afigure
plot(rankV,err)
xlabel('Tensor Decomposition Rank')
ylabel('Error Norm')

%% Show the influence of random init
rankV = 5*ones(200,1);
err = zeros(length(rankV),1);
afigure
for i = 1:length(rankV)
    pKxy = cp_als(tensor(pzx),rankV(i));%,'init','nvecs');
    pKxy = cp_als(tensor(pzx),rankV(i));%,'init','nvecs');
    err(i) = norm(double(pKxy)-pzx);
    %handles =plot2DslicesAroundPoint(pKxy, [0,0], gridXY);
    %ylim(handles{1,1},[-8 8])
    %ylim(handles{2,2},[-8 8])    
end
afigure
plot(rankV,err)
xlabel('Tensor Decomposition Rank')
ylabel('Error Norm')


%% Show the extended mesurement pdf to 6D for the VTOL case
pKxy = cp_als(tensor(pzx),5);
%surf(gridXY{1},gridXY{2},pzx,'EdgeColor','none')
