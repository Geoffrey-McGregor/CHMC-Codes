function [y]= AGSamChi(n,d,p)
PropArea=0;
maxI=0;
IArea=0;
PointsPerInt=0;
y=0;

GG=@(x)x.^(d-1).*exp(-x.^p/p);
LGG=@(x)(d-1)*log(abs(x))-abs(x).^p/p;

zmid=(d-1)^(1/p);

if p==2
    width=3;
end
if p==4
    width=.5;
end
if p==6
    width=0.3;
end

xmin=zmid-width;
xmax=zmid+width;


SubInts=30;
s=linspace(zmid-width,zmid+width,SubInts+1);

for i=1:length(s)-1
    IArea(i)=gamcdf(s(i+1)^p/p,d/p,1)-gamcdf(s(i)^p/p,d/p,1);
    if (s(i)<zmid && s(i+1)<zmid) || (s(i)>zmid && s(i+1)>zmid)
        maxI(i)=max(LGG(s(i+1)),LGG(s(i)));
    else
        maxI(i)=LGG(zmid);
    end
end

TotalArea=sum(IArea);
PropArea=IArea(:)/TotalArea;
PointsPerInt=round(n*PropArea,0);

ss=linspace(xmin,xmax,10000);

for j=1:length(ss)
    for i=1:length(s)-1
        if s(i)<=ss(j) && ss(j)<=s(i+1)
            MaxFunc(j)=maxI(i);
        end
    end
end


count=0;
for i=1:SubInts
    if PointsPerInt(i)>0
        Points=s(i)+rand(PointsPerInt(i),1)*(s(i+1)-s(i));
        Prob1=GG(Points)/maxI(i);
        Prob2=log(rand(PointsPerInt(i),1));

        for j=1:length(Prob1)
            if Prob2(j)<=Prob1(j)
                count=count+1;
                y(count)=Points(j);
            end
        end
    end

end

end

