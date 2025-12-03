% process stellar spectrum for epacris run

clear all
close all

y=importdata('atlas.dat'); % Atlas stellar model 
y=y.data;
w4=y(:,1)*1e+4; % um
f4=y(:,2)*1e-7*1e+4*1e-4; % W/m2/um at the stellar surface
w4r=0.1; % bin to R=1000
while w4r(end)<200
    w4r(end+1)=w4r(end)*(1+1/1000);
end
w4r=w4r';
for i=1:length(w4r)-1
    bidx=(w4>=w4r(i)).*(w4<w4r(i+1));
    f4r(i)=sum(f4.*bidx)/sum(bidx);
end
w4rr = 0.5*w4r(1:end-1) + 0.5*w4r(2:end);

% scale for effective temperature
l = trapz(w4rr,f4r);
lref = 5.67e-8*5214^4;
f4rr = f4r/l*lref;

% flux at 1 AU
f4rr = f4rr*(0.98*695700/1.495978707e+8)^2;

% print out
fp = fopen('55cnc.dat','w');
for i=1:length(w4rr)
    fprintf(fp,'%.6f\t%.6e\n',w4rr(i)*1000,f4rr(i)/1000);
end
fclose(fp);