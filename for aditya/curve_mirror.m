function [polo,dco]=curve_mirror(n,f,eta,pdis,oad,ca,pol)%This code is used to calculate the 
%polarisation after the 2nd parabolic mirror in evelc. The paraboloid is
%placed in the puple plane so i assume that the input is similar to primary
%miror. Run the 1st mirror code 1st and then run this code. U can run this
%code for different tilt angles of input light (eta). But be careful in
%giving pol value. It should match with the eta that u give. I have
%modified the code to match with aditya's evelc design. According to the
%design the parent parabola of the primary is 1mtr while that of the M4 or
%2nd parabolic mirror is 0.541667. This is done to illuminate the M4
%equivalent to M1. More details are written in points book
eta=eta*pi/180;                  
msum=0;
cnt=1;
ct=1;
rln=real(n);
[~,ic,~]=unique(pdis);%i am sampling my curved mirror based on the samples made on primary mirror
ansiz=diff(ic);%this gives the number of theta pointsi should have at a given radius. This
%is again determined by the sampling done on primary mirror. ansiz means
%angle size or number of points
ansiz(end+1)=1;
imn=imag(n);
for h1=(oad-ca):2*ca/(length(ic)-1):(oad+ca)%the circle diameter in meters
    tval=acosd((h1^2+oad^2-ca^2)/(2*oad*h1));
    for tita1=-tval:2*tval/(ansiz(ct)-1):tval
    tita=tita1*pi/180;
i=acos((4*f*cos(eta)-h1*cos(tita)*sin(eta))/sqrt(h1^2+16*f^2));
ii(cnt)=i;
%angmes(cnt)=i*180/pi;
aa=0.5*(rln^2-imn^2-sin(i)^2+sqrt((rln^2-imn^2-sin(i)^2)^2+(2*rln*imn)^2));
bb=0.5*(-rln^2+imn^2+sin(i)^2+sqrt((rln^2-imn^2-sin(i)^2)^2+(2*rln*imn)^2));
a=sqrt(aa);
b=sqrt(bb);
to=atan((2*b*sin(i)*tan(i))/(sin(i)^2*tan(i)^2-(aa+bb)));
toto(cnt)=to;
xx=(aa+bb-2*a*sin(i)*tan(i)+(sin(i)*tan(i))^2)/(aa+bb+2*a*sin(i)*tan(i)+(sin(i)*tan(i))^2);
x=sqrt(xx);


mat=0.5*[1+xx,1-xx,0,0;1-xx,1+xx,0,0;0,0,2*x*cos(to),2*x*sin(to);0,0,-2*x*sin(to),2*x*cos(to) ];

%--------------reflected light direction cosines---------------------------
dx=-(2*cos(tita)*cos(i))/(sqrt(1+16*f^2))-sin(eta);
dy=-(2*sin(tita)*cos(i))/(sqrt(1+16*f^2));
dz=8*f*cos(i)/(sqrt(1+16*f^2))-cos(eta);
dco(cnt,:)=[dx,dy,dz];


%msum=(rot(-tita1)*mat*rot(tita1))+msum;
%msum=mat+msum;
mat=(rot(-tita1)*mat*rot(tita1));
mulmat(:,:,cnt)=mat;
polo(cnt,:)=mat*pol(cnt,:)';%to map the I,Q,U,V of 
%pmat(:,cnt)=mat*[1;0;0;0];                            %unpolarized light source
pang(cnt)=tita1;%the angle at a given point
pdis(cnt)=h1;%the radius at a given point
%using pmat,pang & pdis u can create a I,Q,U,V map of the primary
cnt=cnt+1;
    end
    ct=ct+1;
end
%msum=msum/(cnt-1);
end

function a=rot(th)
rot_angle=th*pi/180;
a=[1,0,0,0;0,cos(2*rot_angle),sin(2*rot_angle),0;0,-sin(2*rot_angle),cos(2*rot_angle),0;0,0,0,1];
end