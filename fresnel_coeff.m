%this is to obtain transmission and reflection coefficients for electric field. the values have to
%be mag squared i;e (num)*(conjugate of num)to use it for intensity(for dielectric)
function [rpl,rper,tpl,tper]=fresnel_coeff(n,thi)%n=R.I of 1st to nth medium, thi=angle of incidence in deg
                                            %from 1st onto 2nd medium
 thi=thi*pi/180;
 rpl=1;rper=1;tpl=1;tper=1;
 
for i=2:length(n)
    ci=cos(thi);%cos of incident angle
    s=n(i-1)/n(i)*sin(thi);%to get sin of transmitted angle
    c=sqrt(1-s^2);
    n1=n(i-1);n2=n(i);
    rpl(i)=-(n2*ci-n1*c)/(n2*ci+n1*c)*rpl(i-1);%generates a 1-D array of reflection coeff at 
    rper(i)=(n1*ci-n2*c)/(n1*ci+n2*c)*rper(i-1);%different boundaries starting from the boundary
    tpl(i)=(2*ci*n1)/(n1*c+n2*ci);              %of 1st to 2nd medium
    tper(i)=(2*ci*n1)/(n1*ci+n2*c);
    thi=asin(s);
end

    