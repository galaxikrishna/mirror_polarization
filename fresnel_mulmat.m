%this function can be used only at a single boundary i;e medium 1 to 2 hence
%only 2 values have to be given for n
function [mat,final_mat]=fresnel_mulmat(n,incident_angle)
[rpl,rper,~,~]=fresnel_coeff(n,incident_angle);
th=angle(rpl/rper);
a=rper(2)^2;b=rpl(2)^2;
mat=0.5*[a+b,a-b,0,0;a-b,a+b,0,0;0,0,2*a*b*cos(th),-2*a*b*sin(th);0,0,2*a*b*sin(th),2*a*b*cos(th)];
final_mat=(rot(-45)*mat*rot(45)*rot(45)*mat*rot(-45))*[1;0;0;1];
mat=(rot(-45)*mat*rot(45))*[1;0;0;1];
end

function a=rot(th)
rot_angle=th*pi/180;
a=[1,0,0,0;0,cos(2*rot_angle),sin(2*rot_angle),0;0,-sin(2*rot_angle),cos(2*rot_angle),0;0,0,0,1];
end