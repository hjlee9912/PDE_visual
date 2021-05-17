function animation(data_seq,Filename,plotON)
% plot density of observation: ypdf.density: Ngrid x tN
%% plot and make an animation 
 
xpdf = data_seq.xpdf; 
ypdf = obsInfo.ypdf; 
tN   = data_seq.steps;

if ~exist(avi_Filename,'file') && plotON ==1    
    v           = VideoWriter(avi_Filename);   % open file for AVI
    v.FrameRate = 4;
    open(v);    
    for t=1:tN
        figure(11);
        subplot(121); plot(ypdf.grid(1:end-1), ypdf.density(:,t),'linewidth',1);
        xlabel('y');  ylabel('PDF of Y_t');   title(['Time t = ', num2str(t)]);
        subplot(122); plot(xpdf.grid(1:end-1), xpdf.density(:,t),'linewidth',1);
        xlabel('x');  ylabel('PDF of X_t');   title(['Time t = ', num2str(t)]);
        %  pause(0.05);
        frame = getframe(gcf);  writeVideo(v,frame);
    end
    close(v);
end
  
end
