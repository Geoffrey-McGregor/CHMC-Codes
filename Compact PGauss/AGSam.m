function [y]= AGSam(n,p)
PropArea=0;
MaxI=0;
IArea=0;
PointsPerInt=0;
y=0;

GG=@(x)exp(-x.^p/p);
xmin=-3;
xmax=3;
SubInts=30;
s=linspace(xmin,xmax,SubInts+1);

for i=1:length(s)-1
    IArea(i)=PGenCDF(s(i+1),p)-PGenCDF(s(i),p);

    if (s(i)<0 && s(i+1)<0) || (s(i)>0 && s(i+1)>0)
        maxI(i)=max(GG(s(i+1)),GG(s(i)));
    else
        maxI(i)=GG(0);
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
        Prob2=rand(PointsPerInt(i),1);

        for j=1:length(Prob1)
            if Prob2(j)<=Prob1(j)
                count=count+1;
                y(count)=Points(j);
            end
        end
    end

end

end

