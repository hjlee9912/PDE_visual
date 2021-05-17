
function plotFS_Gibbs(x,f,SinSeries,n,xrange)
% plot the Gibbs phenomenon
plot(x,f,'k','LineWidth',2);          hold on; 
plot(x,SinSeries,'g-','LineWidth',1);
yline( (max(f)-min(f))*0.09 + max(f)); 
xlim(xrange); 
legend('Function f','Fourier Series','9% overshot line')
title(sprintf('Gibbs Phenomenon (n=%i)',n) );
end