function mulmat=mast_mulmat(alt,azi,n,k)
m12=[1,0.0006430484081,0,0;0.0006430482742,0.999998317,0,0;0,0,0.999989923,0.00412808886;0,0,0.004128089532,-0.99998824];
m3=mirror_mulmat(n,k,45);
m4=mirror_mulmat(n,k,16.85);
m5=mirror_mulmat(n,k,61.85);
m6=mirror_mulmat(n,k,60);
m7=mirror_mulmat(n,k,30);
m8=mirror_mulmat(n,k,60);
m9=mirror_mulmat(n,k,45);

mulmat=m9*rotat(-(alt+azi)/2)*m8*m7*m6*rotat((azi+((alt)+(azi))/2+5))*m5*m4*rotat(-(90-alt))*m3*rotat(90)*m12;
end

function mat=mirror_mulmat(n,k,i)
a=sqrt(0.5*(n*n-k*k-sind(i)*sind(i)+sqrt((n*n-k*k-sind(i))^2+(4*n*n*k*k))));
b=sqrt(0.5*(-n*n+k*k+sind(i)*sind(i)+sqrt((n*n-k*k-sind(i))^2+(4*n*n*k*k))));
xx=(a^2+b^2-2*a*sind(i)*tand(i)+sind(i)^2*tand(i)^2)/(a^2+b^2+2*a*sind(i)*tand(i)+sind(i)^2*tand(i)^2);
tu=atand((2*b*sind(i)*tand(i))/(sind(i)^2*tand(i)^2-(a^2+b^2)));
x=sqrt(xx);
mat=0.5*[1+xx,1-xx,0,0;1-xx,1+xx,0,0;0,0,2*x*cosd(tu),2*x*sind(tu);0,0,-2*x*sind(tu),2*x*cosd(tu)];
end

function rot=rotat(th)
rot=[1,0,0,0;0,cosd(2*th),sind(2*th),0;0,-sind(2*th),cosd(2*th),0;0,0,0,1];
end