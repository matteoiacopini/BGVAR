function y=smth(x,s,p)
% s: smoot/do not smooht
%    if s=1 then smooth
%    otherwise do not smooth
% x: the series
% p: the lags before and after (central smoothing)
y=x;
if s==1
    for i=1:p
        y(p+1:end-p)=y(p+1:end-p)+(x(p+1-i:end-p-i)+...
            x(p+1+i:end-p+i));
    end
    y(p+1:end-p)=y(p+1:end-p)/(2*p+1);
else
    y=x;
end