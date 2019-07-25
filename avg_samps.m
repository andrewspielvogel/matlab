function out = avg_samps(kvh,n)

tic
for i=1:size(kvh.t)/n
    
   out.t(i) = mean(kvh.t((i-1)*n+1:i*n));
   out.mag(i,:) = mean(kvh.mag((i-1)*n+1:i*n,:));
   out.acc(i,:) = mean(kvh.acc((i-1)*n+1:i*n,:));
   out.ang(i,:) = mean(kvh.ang((i-1)*n+1:i*n,:));

end
out.t=out.t';
toc