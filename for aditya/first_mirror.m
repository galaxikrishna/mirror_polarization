function [polo,pang,pdis,dco,mulmat]=first_mirror(n,eta,f,oad,ca,pol)%n is R.I , eta is angle of incidence w.r.t mirror's principle axis
eta=eta*pi/180;                     %f is focal length in meters, pol is input polarisation of light in terms of [I,Q,U,V], oad is 
%offaxis distance of the mirror(here off axis distance is taken from principle axis to optic centre of off axis i;e oad is also the 
%zonal radius of the off axis parabola. ca is the clear apperture size of the offaxis mirror. Note oad and ca must be in meters. All 
%these calculations are made assuming
%off axis parabola is a section of the parent taken from a circle lying on x-axis. polo is the output polarisation and dco is the 
%direction cosine of the output beam. 
msum=0;
cnt=1;
rln=real(n);

imn=imag(n);

for h1=(oad-ca):2*ca/10:(oad+ca)%the circle diameter in meters
    tval=acosd((h1^2+oad^2-ca^2)/(2*oad*h1));
    %length(-tval:1/(2*tval+1):tval)
    for tita1=-tval:10/(tval+10):tval
        
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
polo(cnt,:)=mat*pol';%to map the I,Q,U,V of 
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