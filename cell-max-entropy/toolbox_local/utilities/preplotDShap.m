function [dMatdx,dMatdy] = preplotDShap(dp,near_s,out,x_samp,y_samp,xsample)

  dMatdx=zeros(length(out),length(x_samp),length(y_samp));
  dMatdy=zeros(length(out),length(x_samp),length(y_samp));



  for i=1:(length(x_samp))
    for j=1:(length(y_samp))
      x_trial=[x_samp(i),y_samp(j)];
      dist_test=sqrt((x_trial(1)-xsample(:,1)).^2 + ... 
  		     (x_trial(2)-xsample(:,2)).^2 );
      if (min(dist_test)<1.e-5)
        [val,ii]=min(dist_test);
	
        for ij=1:length(out)
          aux=abs(near_s{ii}-out(ij));
             [v7,n7]=min(aux);
          if v7==0
            dMatdx(ij,i,j)=dp{ii}(n7,1);
            dMatdy(ij,i,j)=dp{ii}(n7,2);
          else
            dMatdx(ij,i,j)=NaN;
            dMatdy(ij,i,j)=NaN;
          end
          %if Mat(ij,i,j)<1.e-4
          %    Mat(ij,i,j)=NaN;
          %end
        end
      else
        dMatdx(:,i,j)=NaN;
        dMatdy(:,i,j)=NaN;
      end
    end
  end



