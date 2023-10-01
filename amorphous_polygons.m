
clear; clc; clf;

aniGif_filename = ['AmorphousPolygonWaveEquation_sub']

pi
MaxPoints=100
Angles=2*pi*[1:MaxPoints]/MaxPoints
radius=1

PolygonPoints = [cos(Angles)',sin(Angles)'];


PolygonPoints'

scatter(PolygonPoints(:,1),PolygonPoints(:,2),'k')
hold;
plot(PolygonPoints(:,1),PolygonPoints(:,2),'-k')

axis equal;
axis([-1.2,1.2,-1.2,1.2]);
axis off;
print('-depsc','Circle_noAxis')


% Calculate Area of the polygon
x_co = 0;
y_co = 0;
for i=1:MaxPoints-1,

	x_co =x_co +( PolygonPoints(i,1) * PolygonPoints(i+1,2));
	y_co =y_co +( PolygonPoints(i+1,1) * PolygonPoints(i,2));
	
end
PolyArea = (x_co - y_co)/2



dTheta = 2*pi/MaxPoints;
ds = sqrt((PolygonPoints(1,1) - PolygonPoints(2,1))^2 + (PolygonPoints(1,2) - PolygonPoints(2,2))^2)
o= 1.;
xm = (MaxPoints*ds)/2;
A=0.5;
alpha=1;
dt= ds;
v=1;
tmax=100;

for i=1:MaxPoints;
	x=i*ds;
gauss_wave(i,1) = A*exp(-((x-xm)^2) / (o^2) );
end
w(:,1) = 0*gauss_wave(:,1);


clf
plot(1:MaxPoints, gauss_wave(:,1))


%      r = u + alpha*(w-w2)
%  LeapFrog(u,w,w2,
% frog(n+1,j) = LeapFrog(frog(n,j),w(n,j+1),w(n,j)

%gauss_wave(:,3) = gauss_wave(:,1);

%{
for i=1:MaxPoints,
	gauss_wave(i,2)= gauss_wave(i,3) - (gauss_wave(i,3)-gaauss_wave(i-1,3));
end
gauss_wave(1,2) = gauss_wave(MaxPoints,2);
%}


for time = 1:tmax,
	

	
	
	for x = 1:MaxPoints-1,
	gauss_wave(x,time+1) = gauss_wave(x,time) + (w(x+1,time)-w(x,time));
	end
	
	gauss_wave(MaxPoints,time+1)=gauss_wave(1,time+1);
	
	
	for x = 1:MaxPoints-1,
	w(x+1,time+1) = w(x+1,time) + (gauss_wave(x+1,time+1)-gauss_wave(x,time+1));
	end
	w(MaxPoints,time+1)=w(1,time+1);

	

	
	
	for i=1:MaxPoints,
	GaussPoly(i,1)=(1-gauss_wave(i,time+1))*cos(Angles(i));
	GaussPoly(i,2)=(1-gauss_wave(i,time+1))*sin(Angles(i));
	end
	
	
	
	
	%pause(0.05)
	%{%

	scatter(GaussPoly(:,1),GaussPoly(:,2),'k')
	hold on;
	plot(GaussPoly(:,1),GaussPoly(:,2),'-k')
	axis equal;
	axis([-1.7,1.7,-1.7,1.7]);
	%axis off;
	hold off;
	%}
	
	%plot(1:MaxPoints, gauss_wave(:,time+1),'g')

	
	
	% Options for exporting to Animated Gif

%{	%
figure(1)
drawnow
frame = getframe(1);
im = frame2im(frame);
[imind,cm] = rgb2ind(im,256);

% Starts the animation with first frame
if  time == 1;
	imwrite(imind,cm,aniGif_filename,'gif','DelayTime',0.0005, 'Loopcount',inf);
% Addes the other frames to the animation
else
	imwrite(imind,cm,aniGif_filename,'gif','DelayTime',0.0005,'WriteMode','append');
	
	% Holds the plot for longer at the end of the simulation
	if time == tmax,
	imwrite(imind,cm,aniGif_filename,'gif','DelayTime',2,'WriteMode','append');
	end

end
%}	
	
	
end


%pause(2);
%plot(1:MaxPoints, gauss_wave(:,tmax),'r')





% Leapfrog





%{
for j=1:12000;
rand_band=0.001;
for i=1:MaxPoints,
	PolygonPoints(i,1)=PolygonPoints(i,1) + (-rand()*rand_band + rand()*rand_band);
	PolygonPoints(i,2)=PolygonPoints(i,2) + (-rand()*rand_band + rand()*rand_band);
end
end

clf;
scatter(PolygonPoints(:,1),PolygonPoints(:,2),'k')
axis equal;
axis([-1.2,1.2,-1.2,1.2]); 
hold all;
plot(PolygonPoints(:,1),PolygonPoints(:,2),'-k')
axis equal;
axis([-1.2,1.2,-1.2,1.2]);
axis off;
print('-depsc', 'JiggledCircle_noAxis_12000')




% Calculate Area of the polygon
x_co = 0;
y_co = 0;
for i=1:MaxPoints-1,

	x_co =x_co +( PolygonPoints(i,1) * PolygonPoints(i+1,2));
	y_co =y_co +( PolygonPoints(i+1,1) * PolygonPoints(i,2));
	
end
PolyArea = (x_co - y_co)/2

%}


	%{
	for i=1:MaxPoints,
		%{%
		if((PolygonPoints(i,1)<=0) && (PolygonPoints(i,2)<=0))
			opx=-1;
			opy=-1;
		elseif((PolygonPoints(i,1)<=0) && (PolygonPoints(i,2)>=0))
			opx=-1;
			opy=1;
		elseif((PolygonPoints(i,1)>=0) && (PolygonPoints(i,2)<=0))
			opx=1;
			opy=-1;
		else
			opx=1;
			opy=1;
		end
		%}%
	GaussPoly(i,1)=PolygonPoints(i,1)+opx*(gauss_wave(i,time+1));
	GaussPoly(i,2)=PolygonPoints(i,2)+opy*(gauss_wave(i,time+1));
	end
	%}
	
	%{
	for i=1:MaxPoints,
		
		if(i==1),
		p_plus = [PolygonPoints(i+1,1),PolygonPoints(i+1,2)];
		p_minus = [PolygonPoints(MaxPoints,1),PolygonPoints(MaxPoints,2)];
		elseif(i==MaxPoints),
		p_plus = [PolygonPoints(1,1),PolygonPoints(1,2)];
		p_minus = [PolygonPoints(i-1,1),PolygonPoints(i-1,2)];
		else
		p_plus = [PolygonPoints(i+1,1),PolygonPoints(i+1,2)];
		p_minus = [PolygonPoints(i-1,1),PolygonPoints(i-1,2)];
		end
		
		%{
		if((PolygonPoints(i,1)<=0) && (PolygonPoints(i,2)<=0))
			opx=-1;
			opy=-1;
		elseif((PolygonPoints(i,1)<=0) && (PolygonPoints(i,2)>=0))
			opx=-1;
			opy=1;
		elseif((PolygonPoints(i,1)>=0) && (PolygonPoints(i,2)<=0))
			opx=1;
			opy=-1;
		else
			opx=1;
			opy=1;
		end
		
		%}
		
		
		mid = [(p_plus(1)+p_minus(1))/2 , (p_plus(2)+p_minus(2))/2];
		slope = (p_minus(2)-p_plus(2))/(p_minus(1)-p_plus(1));
		recip_slope = 1/slope;
		b= recip_slope*mid(1)/mid(2);
		
		x_g = gauss_wave(i,time+1)+ PolygonPoints(i,1);
		y_g = (recip_slope*x_g)-b;
		
	GaussPoly(i,1)=x_g;
	GaussPoly(i,2)=y_g;
		
	end
	
	%}
