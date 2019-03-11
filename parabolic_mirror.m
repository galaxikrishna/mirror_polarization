function rotsum=parabolic_mirror(n)
msum=0;
rotsum=0;
cnt=0;
for incident_angle=-60:60%from principle axis to outer edge of mirror
[rpl,rper,~,~]=fresnel_coeff(n,incident_angle);
th=angle(rpl/rper);
a=rper(2)^2;b=rpl(2)^2;
mat=0.5*[a+b,a-b,0,0;a-b,a+b,0,0;0,0,2*a*b*cos(th),-2*a*b*sin(th);0,0,2*a*b*sin(th),2*a*b*cos(th)];
msum=mat+msum;
cnt=cnt+1;
end
cnt=cnt-1;
for a=1:360 %taking the whole disk of the parabolic mirror and evaluating its effects at the focus
    rotsum=(rot(-a)*msum*rot(a))+rotsum;
    cnt=cnt+1;
end
rotsum=rotsum/(cnt-1);
end

function a=rot(th)
rot_angle=th*pi/180;
a=[1,0,0,0;0,cos(2*rot_angle),sin(2*rot_angle),0;0,-sin(2*rot_angle),cos(2*rot_angle),0;0,0,0,1];
end