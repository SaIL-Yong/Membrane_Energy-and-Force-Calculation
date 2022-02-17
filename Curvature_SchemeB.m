tic
clc
clear all
close all
fid=fopen('/home/BU/dredwan1/Energy and Force Calculation/Coarse_Mesh_Sphere.off');

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
        
        
%         FV.vertices(vertex(1),:)=[vertex(2) vertex(3) vertex(4)];
    end
    x=x';
    y=y';
    z=z';
    
    for i=1:nF
        temp=fgetl(fid);
        tri(i,:)=sscanf(temp,'%*d %d %d %d %*d %*d');
    end
    tri=tri+1; %%% if only an entry is zero/ to avoid non-zero input

end
% fd=@(p) (sum(p.^2,2)+.8^2-.2^2).^2-4*.8^2*(p(:,1).^2+p(:,2).^2);
%     [p,t]=distmeshsurface(fd,@huniform,0.1,[-1,-1,-.25;1,1,.25]);
% tri=t;
% x=p(:,1);
% y=p(:,2);
% z=p(:,3);
%x=zeros(1,896)
%for r=1:896
    
 
tri3d=triangulation(tri,x,y,z);
[bndry_edge,dummy]=freeBoundary(tri3d);
f_normal = faceNormal(tri3d);
f_center = incenter(tri3d);

%%%% angles and edges of each triangle
% for i=1:length(tri(:,1))
%      
%     p1=tri(i,1);
%     p2=tri(i,2);
%     p3=tri(i,3);
%     
%     v1(i,:)=[x(p2)-x(p1),y(p2)-y(p1),z(p2)-z(p1)];
%     v2(i,:)=[x(p3)-x(p2),y(p3)-y(p2),z(p3)-z(p2)];
%     v3(i,:)=[x(p1)-x(p3),y(p1)-y(p3),z(p1)-z(p3)];
%     
%     l_edg(i,1)=norm(v1(i,:));
%     l_edg(i,2)=norm(v2(i,:));
%     l_edg(i,3)=norm(v3(i,:));
%     
%     ang_tri(i,1)=acos(dot(v1(i,:)/l_edg(i,1),-v3(i,:)/l_edg(i,3)));
%     ang_tri(i,2)=acos(dot(-v1(i,:)/l_edg(i,1),v2(i,:)/l_edg(i,2)));
%     ang_tri(i,3)=pi-(ang_tri(i,1)+ang_tri(i,2));
%     
% end
%%% area of all triangles
for i=1:length(tri(:,1))
    v1=[x(tri(i,2))-x(tri(i,1)),y(tri(i,2))-y(tri(i,1)),z(tri(i,2))-z(tri(i,1))];
    v2=[x(tri(i,3))-x(tri(i,1)),y(tri(i,3))-y(tri(i,1)),z(tri(i,3))-z(tri(i,1))];
    area_tri(i,1)=0.5*norm(cross(v1,v2));
    
end

V = vertexAttachments(tri3d);
%%% Finding the total area of the Neighbor
area_total=zeros(length(x),1);
for i=1:length(x)
    neighb=length(V{i});
    
    for ii=1:neighb
        area_total(i)=area_total(i)+area_tri(V{i}(ii));
    end
    
    area_total(i)=area_total(i)/3;%%%% divided by 3 to get Ai
end


sum_of_lt=zeros(length(x),1); %%% summation of length * theta
% 
E_e=edges(tri3d);
for i=1:length(x)
   
    N=length(V{i}); %% attached triangle 
    for j=1:N
        K=tri3d.ConnectivityList(V{i}(j),:)%%%% attached vertices
%         E=edgeAttachments(tri3d,i,K(1));
            
            for ii=1:3
                if K(ii)==i  %%%% avoiding zero values
                    continue
                else
                    l=sqrt(((x(i)-x(K(ii)))^2)+((y(i)-y(K(ii)))^2)+((z(i)-z(K(ii)))^2));
                    
                    E=edgeAttachments(tri3d,i,K(ii));
                    
                    theta=acos(dot(f_normal(E{1}(1),:),f_normal(E{1}(2),:)));
                    sum_of_lt(i)=sum_of_lt(i)+l*theta;
                end

                
            end

    end
       
   
end
E_e=edges(tri3d);

% for i=length(x)
%     
%     N=length(V{i});
%     for j=1:N
%        length_edges=sqrt(((x(i)-x(E(i,:,j))).^2)+((y(i)-y(E(:,j))).^2)+((z(i)-z(E(:,j))).^2))
%     end    
% end 


Mean_Curvature=zeros(length(x),1);

for iii=1:length(x)
    Mean_Curvature(iii)=sum_of_lt(iii)/(8*area_total(iii));
end

% fd=@(p) (sum(p.^2,2)+.8^2-.2^2).^2-4*.8^2*(p(:,1).^2+p(:,2).^2);
%     [p,t]=distmeshsurface(fd,@huniform,0.1,[-1,-1,-.25;1,1,.25]);
% tri=t;
% x=p(:,1);
% y=p(:,2);
% z=p(:,3);
%x=zeros(1,896)
%for r=1:896
    
 


a=max(Mean_Curvature)
b=min(Mean_Curvature)
% for i=1:length(x)
%     if (x(i)^2+y(i)^2>10^2)
%         if z(i)>=0 
%             theta(i)=asin(z(i)./2.5);
%         else
%             theta(i)=asin(z(i)./2.5)+2*pi;
%         end
%     else
%         if z(i)>=0 
%             theta(i)=pi-asin(z(i)./2.5);
%         else
%             theta(i)=pi-asin(z(i)./2.5);
%         end
%     end
% %     
% %     theta(i)=acos(z(i)/sqrt(x(i)^2+y(i)^2+z(i)^2));    
% end

% 
% theta_min=min(theta)
% theta_max=max(theta)
% figure
% plot(theta,Mean_Curvature,'bo','MarkerSize',4,'MarkerFaceColor',[0.5,0.5,0.5])
% 
% R=10; % outer radius of torus
% r=2.5; % inner tube radius
% th=linspace(0,2*pi,100); % partitions along perimeter of the tube 
% % phi=linspace(0,2*pi,100); %  partitions along azimuth of torus
% 
% phi=pi/4;


% 
% Mean_Curv=(-(R+2*r*cos(th))./(2*r*(R+r*cos(th))));
% % 
% c=max(Mean_Curv)
% d=min(Mean_Curv)


% % plot(th,Mean_Curv_P);
% hold on
% plot(th,Mean_Curv,'Linewidth',2.5);
% xlabel('Theta');ylabel('Mean_ Curvature')
% title('Torus_ Curvature')


figure
trisurf(tri3d,[0.8 0.8 1.0]);
axis equal
hold on
pp= patch('Faces',tri,'Vertices',[x,y,z],'FaceVertexCData',Mean_Curvature,'FaceColor','interp','EdgeColor','none');
axis equal
colorbar


% %%%% To check the average of all Mean Curvature 
% avg=mean(Mean_Curvature)
% std=std(Mean_Curvature)
% 
% a=max(Mean_Curvature)
% b=min(Mean_Curvature)
% 
% 
% for i=1:length(x)
%     if (x(i)^2+y(i)^2>10^2)
%         if z(i)>=0 
%             theta(i)=asin(z(i)./2.5);
%         else
%             theta(i)=asin(z(i)./2.5)+2*pi;
%         end
%     else
%         if z(i)>=0 
%             theta(i)=pi-asin(z(i)./2.5);
%         else
%             theta(i)=pi-asin(z(i)./2.5);
%         end
%     end
% %     
% %     theta(i)=acos(z(i)/sqrt(x(i)^2+y(i)^2+z(i)^2));    
% end
% min(theta)
% max(theta)
% figure
% plot(theta,Mean_Curvature,'bo','MarkerSize',4,'MarkerFaceColor',[0.5,0.5,0.5])
% 
% R=10; % outer radius of torus
% r=2.5; % inner tube radius
% th=linspace(0,2*pi,100); % partitions along perimeter of the tube 
% % phi=linspace(0,2*pi,100); %  partitions along azimuth of torus
% 
% phi=pi/4;
% 
% % x=(R+r*cos(th))*cos(phi)
% % z=r*sin(th);
% % y=(R+r*cos(th))*sin(phi);
% % surf(x,y,z); % plot surface
% 
% 
% Mean_Curv=-(R+2*r*cos(th))./(2*r*(R+r*cos(th)));
% 
% max=max(Mean_Curv)
% min=min(Mean_Curv)
% 
% 
% % plot(th,Mean_Curv_P);
% hold on
% plot(th,Mean_Curv,'Linewidth',2.5);
% xlabel('Theta');ylabel('Mean_ Curvature')
% title('Torus_ Curvature')
% 
% 
% figure
% trisurf(tri3d,[0.8 0.8 1.0]);
% axis equal
% hold on
% pp= patch('Faces',tri,'Vertices',[x,y,z],'FaceVertexCData',Mean_Curvature,'FaceColor','interp','EdgeColor','none');
% axis equal
% colorbar
% 
% 
trisurf(tri3d,[0.8 0.8 1.0]);
axis equal
hold on
pp= patch('Faces',tri,'Vertices',[x,y,z],'FaceVertexCData',Mean_Curvature,'FaceColor','interp','EdgeColor','none');
axis equal
colorbar
% 
% toc
