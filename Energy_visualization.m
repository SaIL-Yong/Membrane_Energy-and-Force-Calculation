clc
clear all
close all


fid=fopen('/home/BU/dredwan1/Energy and Force Calculation/BigSphere.off');

while feof(fid)==0
    temp=fgetl(fid);
    temp=fgetl(fid);
    [nV,nF,nE]=strread(temp, '%d %d %d');
    
    FV.vertices=zeros(nV,3);
    FV.faces=zeros(nF,3);
    
    for i=1:nV
        temp=fgetl(fid);
        vertex=sscanf(temp, '%g %g %g');
        x(i)=vertex(1);
        y(i)=vertex(2);
        z(i)=vertex(3);    
    end
     x=x';
     y=y';
     z=z';
    
    for i=1:nF
        temp=fgetl(fid);
        tri(i,:)=sscanf(temp, '%*d %d %d %d %*d %*d');
    end
    tri=tri+1;

end
fid2=fopen('/home/BU/dredwan1/Energy and Force Calculation/Force_txt_finer.txt');

while feof(fid2)==0
    for i=1:length(x)
       temp=fgetl(fid2);
       Energy_bending(i,:)=sscanf(temp, '%f');
    end
    
end


tri3d=triangulation(tri,x,y,z);



figure
trisurf(tri3d,[0.8 0.8 1.0]);
axis equal
hold on
pp= patch('Faces',tri,'Vertices',[x,y,z],'FaceVertexCData',Energy_bending,'FaceColor','interp','EdgeColor','none');
axis equal

colorbar
title('Bending_Energy');
caxis([0.01 0.05]);
title('Bending Energy');

figure2
trisurf(tri3d,[0.8 0.8 1.0]);
axis equal
hold on


% 
% figure
% trisurf(tri3d,[0.8 0.8 1.0]);
% axis equal
% hold on
% pp= patch('Faces',tri3d,'Vertices',[x,y,z],'FaceVertexCData',Eb,'FaceColor','interp','EdgeColor','none');
% axis equal
% colorbar

