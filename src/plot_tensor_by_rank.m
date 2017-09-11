function plot_tensor_by_rank(F)

d = ndims(F);
r = ncomponents(F);

xlabel_list = {'Trial','x','y','z','Features'};

figure;
for i = 1:d
    F_U = F.U{i};
    n = size(F_U,1);
    
    ylimrange = [min(min(F_U)) max(max(F_U))];
    
    for j = 1:r
        subplot(r,d,i+(j-1)*d)
        scatter(1:size(F_U,1), F_U(:,j), '.')
        hold on
        plot([1, size(F_U,1)],[0 0],'k','linewidth',0.5)
        hold off
        xlim([1 n]);
        ylim(ylimrange);
        
        if j == r
            xlabel(xlabel_list(i))
        else
            set(gca,'xticklabel',[])
        end
        
        if i == 1
            text(-40, 0, num2str(F.lambda(j),'%.0f'));
        end
    end   
end

lambda = F.lambda;

figure;

for i = 1:d
    F_U = F.U{i};
    n = size(F_U,1);
    
%     ylimrange = [min(min(F_U)) max(max(F_U))];
    
    subplot(d,1,i)
    plot(bsxfun(@mtimes,F_U,lambda'))
    xlim([1 n]);
    xlabel(xlabel_list(i))
%     ylim(ylimrange);
      
end


% title(F.lambda)

end