% compute E2

clear all

z = 10.^(-3:0.01:5);

for i=1:length(z)
    zz=z(i);
    xmax = max(100/zz,100);
    x=linspace(1,xmax,1000000);
    f=x.^(-2).*exp(-zz*x);
    E(i)=trapz(x,f);
end

fp=fopen('E2.dat','w');
for i=1:length(z)
    fprintf(fp,'%2.6e\t%2.6e\n',z(i),E(i));
end
fclose(fp);