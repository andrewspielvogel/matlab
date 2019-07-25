function out = calc_mag_rms(out,phins)

start = 20*60;
t = out.t-out.t(1);
indices = find(t>start);
index = indices(1);
out.mag = real(out.mag);
if isreal(out.mag)
    out.heading = calc_heading(out.t,out,phins);
    error = out.heading.att_error(:,index:end);
    error(isnan(error)) = [];
    rms_error = rms(error);
else
    rms_error=20;
end

out.rms_error = rms_error;