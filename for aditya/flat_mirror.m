function [polo,dco]=flat_mirror(n,th,dci,pol)%Here n is complex index of refraction. 'th' is the
%tilt in the flat mirror. U have to specify tilt as direction cosine of the normal i;e
%the angle which the normal to the mirror makes with x,y,z axis. So 'th' is
%a vector of 3 elements. dci is the direction cosine of the incident beam.
%if dci is n x 3 matrix then pol is a n x 4 matrix containing the
%polarisation information of each ray. polo is the output polarisation for
%a given input polarisation pol and dco is the direction cosine of the
%reflected ray.
rln=real(n);
imn=imag(n);

for i=1:size(dci,1)
    inc=acos(dot(th,dci(i,:)));%as i am cosidering direction cosines their dot
    %product itself is the cos of angle. This is because their magnitude = 1
    aa=0.5*(rln^2-imn^2-sin(inc)^2+sqrt((rln^2-imn^2-sin(inc)^2)^2+(2*rln*imn)^2));
    bb=0.5*(-rln^2+imn^2+sin(inc)^2+sqrt((rln^2-imn^2-sin(inc)^2)^2+(2*rln*imn)^2));
    a=sqrt(aa);
    b=sqrt(bb);
    to=atan((2*b*sin(inc)*tan(inc))/(sin(inc)^2*tan(inc)^2-(aa+bb)));
    xx=(aa+bb-2*a*sin(inc)*tan(inc)+(sin(inc)*tan(inc))^2)/(aa+bb+2*a*sin(inc)*tan(inc)+(sin(inc)*tan(inc))^2);
    x=sqrt(xx);
    


    mat=0.5*[1+xx,1-xx,0,0;1-xx,1+xx,0,0;0,0,2*x*cos(to),2*x*sin(to);0,0,-2*x*sin(to),2*x*cos(to) ];
    dco(i,:)=dci(i,:)-(2.*dot(dci(i,:),th).*th);%reflected ray vector or its direction cosine is r=i-2(i dot n)*n
    polo(i,:)=mat*pol(i,:)';
end

    