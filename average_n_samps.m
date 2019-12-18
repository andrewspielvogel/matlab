function out = average_n_samps(samp, n)

num = size(samp.t,1);

new_num = floor(num/n);

for i=1:new_num
    
   out.mag(i,:) = mean(samp.mag(1+(i-1)*n:i*n,:),1); 
   %out.dmag(i,:) = mean(samp.dmag(1+(i-1)*n:i*n,:),1); 
   out.ang(i,:) = mean(samp.ang(1+(i-1)*n:i*n,:),1); 
   out.acc(i,:) = mean(samp.acc(1+(i-1)*n:i*n,:),1); 
   out.t(i) = mean(samp.t(1+(i-1)*n:i*n));
    
end
out.t = out.t';