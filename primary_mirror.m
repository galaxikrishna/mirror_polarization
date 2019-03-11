function [pmat,pang,pdis]=primary_mirror(n,eta,f,pol)%n is R.I , eta is angle of incidence w.r.t mirror's principle axis
eta=eta*pi/180;                     %f is focal length in meters, pol is input polarisation of light in terms of [I,Q,U,V]
msum=0;
cnt=1;
for h1=0:0.01:1%the circle diameter in meters
    for tita1=0:360/(2*pi*h1*1000):360
    tita=tita1*pi/180;
i=acos(2*(4*f*cos(eta)-h1*cos(tita)*sin(eta))/sqrt(h1^2+16*f^2));
r=asin(sin(i)/n);
rp=tan(i-r)/tan(i+r);
rs=-sin(i-r)/sin(i+r);
th=angle(rp/rs);
a=rs^2;b=rp^2;
mat=0.5*[a+b,a-b,0,0;a-b,a+b,0,0;0,0,2*a*b*cos(th),-2*a*b*sin(th);0,0,2*a*b*sin(th),2*a*b*cos(th)];
%msum=(rot(-tita1)*mat*rot(tita1))+msum;
%msum=mat+msum;
pmat(:,cnt)=(rot(-tita1)*mat*rot(tita1))*pol';%to map the I,Q,U,V of 
%pmat(:,cnt)=mat*[1;0;0;0];                            %unpolarized light source
pang(cnt)=tita1;%the angle at a given point
pdis(cnt)=h1;%the radius at a given point
%using pmat,pang & pdis u can create a I,Q,U,V map of the primary
cnt=cnt+1;
    end
end
%msum=msum/(cnt-1);
end

function a=rot(th)
rot_angle=th*pi/180;
a=[1,0,0,0;0,cos(2*rot_angle),sin(2*rot_angle),0;0,-sin(2*rot_angle),cos(2*rot_angle),0;0,0,0,1];
end