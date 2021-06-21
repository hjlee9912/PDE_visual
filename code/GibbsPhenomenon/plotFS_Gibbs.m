function plotFS_Gibbs(x,f,SinSeries,n,xrange)
% plot the Gibbs phenomenon
plot(x,f,'w','LineWidth',2);          hold on; 
plot(x,SinSeries,'y-','LineWidth',2);
yline( (max(f)-min(f))*0.09 + max(f), 'w','LineWidth',1); 
xlim(xrange); 
legend('\color{white} Function f','\color{white} Fourier Series','\color{white} 9% overshot line','Color',[0 0 0],'location','best','EdgeColor',[1 1 1])
title(sprintf('Gibbs Phenomenon (N=%i)',n) );
end
