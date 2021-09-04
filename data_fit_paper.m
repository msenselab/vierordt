% fit model to data

load('rep_randwalk.mat')
nsub=max(allArray.NSub);

secondhalf=false;
if secondhalf
    idx=allArray.NT>600 | (allArray.NT>200 & allArray.NT<400);
    %idx=~idx; % take first instead of second half
    idx=allArray.RandLevel==2 | idx; % take half, but only for random walk
    allArray=allArray(idx,:);
end

% round physical presented duration (100 ms)
allArray.duration = round(allArray.dur *10)/10;
durs=unique(allArray.duration(allArray.NSub==1));
allArray.sim=zeros(size(allArray.duration));

% take only valid data
% allArray = allArray(allArray.valid==1,:);
% alternative, leaves temporal structure
allArray.repDur(allArray.valid==0)=NaN;

cond={'random walk','randomized'};
parfit=zeros(nsub,2,1);
linfit=zeros(nsub,2,2);
simfit=zeros(nsub,2,2);


tic
for j=2:-1:1
    figure
    for i=1:nsub
        idx=(allArray.NSub==i) & (allArray.RandLevel==j);
        
        % original model, just 1 parameter
        [px,ci,resnorm,simres]=fitdata1pv([allArray.dur(idx) allArray.repDur(idx)]);
        parfit(i,j,1)=px;
        err=simres-allArray.repDur(idx);
        idn=isfinite(err);
        
        stid=allArray.dur(idx);
        repd=allArray.repDur(idx);
        linfit(i,j,:)=polyfit(stid(idn),repd(idn),1);
        
        if j==1
            % predict response from parameter fitted for other condition
            [~, ~, ~, simres, ~]=kmodel1pv(parfit(i,3-j,1),[allArray.dur(idx) allArray.repDur(idx)]);
        else
            % simulate response with fitted parameters for this condition
            [~, ~, ~, simres, ~]=kmodel1pv(parfit(i,j,1),[allArray.dur(idx) allArray.repDur(idx)]);
        end
        
        allArray.sim(idx)=simres;
        
        simfit(i,j,:)=polyfit(stid(idn),simres(idn),1);
        
        subplot(4,4,i)
        prcterr=(allArray.repDur(idx)-allArray.dur(idx))./allArray.dur(idx)*100;
        prcterrsim=(simres-allArray.dur(idx))./allArray.dur(idx)*100;
        
        plot(allArray.dur(idx),allArray.repDur(idx),'.',allArray.dur(idx),simres,'.')
        ylim([0 3])
        
        %plot(allArray.dur(idx),prcterr,'.',allArray.dur(idx),prcterrsim,'.')
        %ylim([-50 50])
        
        meanprct=grpstats(prcterr,allArray.duration(idx));
        [meanrepdur,stdrepdur]=grpstats(allArray.repDur(idx),allArray.duration(idx),{'mean','std'});
        durs=unique(allArray.duration(idx));
        meanprctsim=grpstats(prcterrsim,allArray.duration(idx));
        [meanrepsim,stdrepsim]=grpstats(simres,allArray.duration(idx),{'mean','std'});
    end
end
toc


%%
i=nsub;
% i=5;
% i=8;
figure('name',['data from subject ' int2str(i)])
%clr=colororder;
% randomized
j=2;
idx=(allArray.NSub==i) & (allArray.RandLevel==j);
subplot(1,2,1)
hold on
plot([0 2],linfit(i,j,1)*[0 2]+linfit(i,j,2),'linewidth',2)
plot([0 2],simfit(i,j,1)*[0 2]+simfit(i,j,2),'--','linewidth',2)
plot(allArray.dur(idx),allArray.repDur(idx),'.')
plot(allArray.dur(idx),allArray.sim(idx),'.')
plot([0 3],[0 3],'--k')
hold off
xlim([0 2])
ylim([0 2.5])
xlabel('stimulus duration (s)')
ylabel('reproduced duration (s)')
legend('experiment','model simulation')
title('randomized condition')
set(gca,'Fontsize',16)
text('Parent',gca,'FontSize',36,'String','A','Position',[-0.4 2.5 0]);

% calculate error (response-actual)
errA=allArray.repDur(idx)-allArray.dur(idx);
presdurA=allArray.dur(idx);
simerrA=allArray.sim(idx)-allArray.dur(idx);

% now for random walk
j=1;
idx=(allArray.NSub==i) & (allArray.RandLevel==j);
subplot(1,2,2)
hold on
plot([0 2],linfit(i,j,1)*[0 2]+linfit(i,j,2),'linewidth',2)
plot([0 2],simfit(i,j,1)*[0 2]+simfit(i,j,2),'--','linewidth',2)
plot(allArray.dur(idx),allArray.repDur(idx),'.')
plot(allArray.dur(idx),allArray.sim(idx),'.')
plot([0 3],[0 3],'--k')
hold off
xlim([0 2])
ylim([0 2.5])
xlabel('stimulus duration (s)')
ylabel('reproduced duration (s)')
legend('experiment','model prediction')
title('random walk condition')
set(gca,'Fontsize',16)
text('Parent',gca,'FontSize',36,'String','B','Position',[-0.4 2.5 0]);
set(gcf,'Position',[560   556   836   392])

% calculate error (response-actual)
errB=allArray.repDur(idx)-allArray.dur(idx);
presdurB=allArray.dur(idx);
simerrB=allArray.sim(idx)-allArray.dur(idx);

%%
figure('name',['data from subject ' int2str(i) ' over trial'])
subplot(2,1,1)
j=2;
idx=(allArray.NSub==i) & (allArray.RandLevel==j);
plot([allArray.dur(idx),allArray.repDur(idx),allArray.sim(idx)],'.-','Markersize',10)
xlabel('trial')
ylabel('duration (s)')
legend('stimulus','reproduction','simulation')
set(gca,'Fontsize',16)
title('randomized')
subplot(2,1,2)
j=1;
idx=(allArray.NSub==i) & (allArray.RandLevel==j);
plot([allArray.dur(idx),allArray.repDur(idx),allArray.sim(idx)],'.-','Markersize',10)
xlabel('trial')
ylabel('duration (s)')
legend('stimulus','reproduction','prediction')
set(gca,'Fontsize',16)
title('random walk')

%%
allArray.repError = allArray.repDur - allArray.dur;
allArray.simError = allArray.sim - allArray.dur;

avgdata = grpstats(allArray, {'NSub','RandLevel','duration'},{'mean','sem'}, 'DataVars','repError');
avgdatadur = grpstats(allArray, {'NSub','RandLevel','duration'},{'mean','sem'}, 'DataVars','repDur');

avgsim = grpstats(allArray, {'NSub','RandLevel','duration'},{'mean','sem'}, 'DataVars','simError');
avgsimdur = grpstats(allArray, {'NSub','RandLevel','duration'},{'mean','sem'}, 'DataVars','sim');

avgdata.prcterr=avgdata.mean_repError./avgdata.duration*100;
avgsim.prcterr=avgsim.mean_simError./avgsim.duration*100;

gavgdata = grpstats(avgdata, {'RandLevel','duration'},{'mean','sem'}, 'DataVars','mean_repError');
gavgdatadur = grpstats(avgdatadur, {'RandLevel','duration'},{'mean','sem'}, 'DataVars','mean_repDur');
gavgsim = grpstats(avgsim, {'RandLevel','duration'},{'mean','sem'}, 'DataVars','mean_simError');
gavgsimdur = grpstats(avgsimdur, {'RandLevel','duration'},{'mean','sem'}, 'DataVars','mean_sim');

gavgdata.prcterr=gavgdata.mean_mean_repError./gavgdata.duration*100;
gavgsim.prcterr=gavgsim.mean_mean_simError./gavgsim.duration*100;

gavgprct=grpstats(avgdata, {'RandLevel','duration'},{'mean','sem'}, 'DataVars','prcterr');
gavgsimprct=grpstats(avgsim, {'RandLevel','duration'},{'mean','sem'}, 'DataVars','prcterr');


%% plot for Vierordt paper
figure

idx=gavgdata.GroupCount>2;
subplot(1,2,1)
hold on
errorbar(gavgprct.duration(gavgprct.RandLevel == 1 & idx)-0.01, gavgprct.mean_prcterr(gavgprct.RandLevel == 1 & idx),gavgprct.sem_prcterr(gavgprct.RandLevel == 1 & idx),'k-o','MarkerFaceColor','k','MarkerSize',10);
errorbar(gavgprct.duration(gavgprct.RandLevel == 2 & idx)+0.01, gavgprct.mean_prcterr(gavgprct.RandLevel == 2 & idx),gavgprct.sem_prcterr(gavgprct.RandLevel == 2 & idx),'k-o','MarkerSize',10);
plot([0 10],[0 0],'--k')
xlim([0 2])
ylim([-30 60])
title('experiment')
legend('random walk','randomized')
xlabel('mean target time (s)')
ylabel('percentage error of reproduction')
set(gca,'Fontsize',16)
text('Parent',gca,'FontSize',36,'String','A','Position',[-0.4 63 0]);
hold off

subplot(1,2,2)
hold on
errorbar(gavgsimprct.duration(gavgsimprct.RandLevel == 1 & idx)-0.01, gavgsimprct.mean_prcterr(gavgsimprct.RandLevel == 1 & idx),gavgsimprct.sem_prcterr(gavgsimprct.RandLevel == 1 & idx),'k-o','MarkerFaceColor','k','MarkerSize',10);
errorbar(gavgsimprct.duration(gavgsimprct.RandLevel == 2 & idx)+0.01, gavgsimprct.mean_prcterr(gavgsimprct.RandLevel == 2 & idx),gavgsimprct.sem_prcterr(gavgsimprct.RandLevel == 2 & idx),'k-o','MarkerSize',10);
plot([0 10],[0 0],'--k')
xlim([0 2])
ylim([-30 60])
title('simulation')
xlabel('mean target time (s)')
ylabel('percentage error of reproduction')
set(gca,'Fontsize',16)
text('Parent',gca,'FontSize',36,'String','B','Position',[-0.4 63 0]);
hold off
set(gcf,'Position',[560   556   836   392])

