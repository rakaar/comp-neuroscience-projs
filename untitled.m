fr=50;
dt=0.001;
nBins=18000;
spikeMat = zeros(1, nBins);
rng('shuffle');
for i=1:nBins
    if rand<=1-exp(-fr*dt);
        spikeMat(1,i)=1;
    else
        spikeMat(1,i)=0;
    end
end
spkt=find(spikeMat==1);
ds=diff(spkt);