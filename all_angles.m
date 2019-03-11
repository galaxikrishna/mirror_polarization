%this function plots fresnell coeff for intensity and returns their values
%for electric field
function [rpl,rper,tpl,tper]=all_angles(n,medium)%give the R.I of the multiple layers and specify which medium 
                             %the variation of co-eff has to be seen in
  cnt=1;                           %(medium)
  for i=0:90                           
 [a,b,c,d]=fresnel_coeff(n,i);
 
 rpl(cnt)=a(medium);
 rper(cnt)=b(medium);
 tpl(cnt)=c(medium);
 tper(cnt)=d(medium);
 cnt=cnt+1;
  end
  
  x=1:cnt-1;
subplot(1,2,1)
  plot(x,rpl.*conj(rpl),x,rper.*conj(rper))
subplot(1,2,2)  
  plot(x,tpl.*conj(tpl),x,tper.*conj(tper))   
      