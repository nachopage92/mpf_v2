function x = QR_solver(A,b)
% function x = QR_solver(A,b)
%
% solve a system A*x=b, where A is m-by-n, and b is the r.h.s.

[Q,R,E]= qr(A);
  qb     = Q'*b;
  N      = length(R);
  ua     = zeros(N,1);

  Nr     = rank(R);
  
  % FIRST: the indetermined part is solvedI
  k    =0;
  indet=0;

  for i=N:-1:Nr
    for j=i:N
      if abs(R(i,j))>1.e-08
        if N-j==k %OK
          ua(j) = (qb(i) - R(i,j+1:N)*ua(j+1:N))/R(i,j);
          k = k+1;
        else
          kk = k;
          for n=j+1:N-kk
            ua(n) = 0;
            k     = k+1;
            indet = indet+1;
          end
          ua(j) = (qb(i) - R(i,j+1:N)*ua(j+1:N))/R(i,j);
          k = k+1;
        end
        break
      end
    end
  end

%  fprintf(1,'\n\tN=%d   i_cut=%d   diff=%d    indet=%d\n',N,Nr,N-Nr,indet);

  % SECOND: final computation for the well possed part
  for i=Nr-1:-1:1
    ua(i) = (qb(i) - R(i,i+1:N)*ua(i+1:N))/R(i,i);
  end

  x=E*ua;