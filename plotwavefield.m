clear;
close all;
error = 0;
i = 1;
filename1=(['Slice',num2str(i,'%5.5d'),'.dat']);
% filename2=(['Hx_EyT',num2str(i,'%5.5d'),'.txt']);
while fopen(filename1) ~= 0
    
    fid1 = fopen(filename1);
%     fid2 = fopen(filename2);
    if fid1 == -1
        disp('No more file')
        break;
    end
    Ey1=fscanf(fid1,'%e',[101,126]);
%     Ey2=fscanf(fid2,'%e',[72,71]);
%     subplot(1,2,1)
    imagesc(-Ey1')
    axis image
    colorbar;
%     subplot(1,2,2)
%     imagesc(-Ey2')
%     axis image
%     colorbar;
    %subplot(1,3,3)
    %imagesc(-Ey1'+Ey2'[:,)
    %caxis([-0.00001 0.00001])
%     axis image
%     colorbar;
    i=i+1
    
    filename1=(['Slice',num2str(i,'%5.5d'),'.dat']);
%     filename2=(['Hx_EyFie',num2str(i,'%5.5d'),'.txt']);
    pause(0.5)
end

