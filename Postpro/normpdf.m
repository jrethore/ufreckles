function ygam=normpdf(x,ugam,tgam)

ygam=exp(-((x-ugam)/tgam).^2/2)/(sqrt(2*pi)*tgam);