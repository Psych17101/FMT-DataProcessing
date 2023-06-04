%Calculating Size of Matrices:
[i_im,j_im] = size(im1);
im=floor((i_im-win_size)/(win_size-1)); %Number of I.W.s in x direction
jm=floor((j_im-win_size)/(win_size-1))-1; %Number of I.W.s in y direction

x(1:im,1:jm)=0.;y(1:im,1:jm)=0.;u(1:im,1:jm)=0.;v(1:im,1:jm)=0.;vel(1:im,1:jm)=0.;
for j=1:jm
    for i=1:im
        x(i,j)=(i*(win_size-1)/2)*scale;
        y(i,j)=(j*(win_size-1)/2)*scale;
        u(i,j)=vecx(i*win_size,j*win_size)*scale/tscale;
        v(i,j)=vecy(i*win_size,j*win_size)*scale/tscale;
        vel_new(i,j)=vec(i*win_size,j*win_size)*scale/tscale;
    end
end

  u_n(iterx,itery) = u;
        v_n(iterx,itery) = v;
        u_new(iterx,floor((itery+1)/2)) = u;
        v_new(iterx,floor((itery+1)/2)) = v;


        %       
%         x(i:i+win_size-1, j:j+win_size-1) = double(((i:i+win_size-1)*(win_size-1)/2)*scale);
%         y(i:i+win_size-1, j:j+win_size-1) = ((j:j+win_size-1)*(win_size-1)/2)*scale;
        
