function plotFS_TBTD(x,f,Series,n,xrange,yrange,titlee,leg)
% plot the Gibbs phenomenon
plot(x,f,'k','LineWidth',2);          hold on; 
plot(x,Series,'m-','LineWidth',2);
xlim(xrange); 
ylim(yrange);
legend(leg,'Fourier Series','location','southeast')
title(sprintf(append(titlee,' (n=%i)'),n) );
end