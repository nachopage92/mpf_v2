function [ufield] = preplotField(u_field,out,x_samp,y_samp,xsample)

  ufield=zeros(length(x_samp),length(y_samp));



  for i=1:(length(x_samp))
    for j=1:(length(y_samp))
      x_trial=[x_samp(i),y_samp(j)];
      dist_test=sqrt((x_trial(1)-xsample(:,1)).^2 + ... 
  		     (x_trial(2)-xsample(:,2)).^2 );
      if (min(dist_test)<1.e-5)
        [val,ii]=min(dist_test);
          
        ufield(i,j)=u_field(ii);

      else
        ufield(i,j)=NaN;      
      end
    end
  end



