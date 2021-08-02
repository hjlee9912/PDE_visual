function im_to_gif_dt(filename,im, dt)
nframes = numel(im);
for idx=1:nframes
    [A,map] = rgb2ind(im{idx},256);
    if idx == 1
          imwrite(A,map,filename,'gif', 'Loopcount',inf,'DelayTime',dt);
    else
          imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',dt);
    end
end
end