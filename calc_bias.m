function [bias, data] = calc_bias(file,phinsrph,lat)



lat = lat*pi/180;
Ren = [-sin(lat),0,-cos(lat);0,1,0;cos(lat),0,-sin(lat)];
r = 6371*1000;
a_e = [cos(lat);0;sin(lat)] - (7.292150/100000)^2*cos(lat)*[r;0;0]/9.81;
a_n = Ren'*a_e;
w_E_e = [0;0;1]*7.292150/100000;
w_E_n = Ren'*w_E_e;


fid = fopen(file);
data = textscan(fid,'%s %d/%d/%d %d:%d:%f %f %f %f,%f,%f,%f,%f,%f,%f,%f,%f, %f, %f, %f,%f, %d, %d, %d, %d, %d, %d');
fclose(fid);

bias.w_b = [mean(cell2mat(data(10)));mean(cell2mat(data(11)));mean(cell2mat(data(12)))] - rph2R([-pi/2;0;pi/2])'*rph2R(phinsrph)'*w_E_n;
bias.a_b = [mean(cell2mat(data(13)));mean(cell2mat(data(14)));mean(cell2mat(data(15)))] - rph2R([-pi/2;0;pi/2])'*rph2R(phinsrph)'*a_n;
bias.a_b = bias.a_b*9.81;