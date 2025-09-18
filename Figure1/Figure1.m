clear 
clc
Table=readtable('Source Data/SynCom 67 Data.xlsx','Sheet',1);

SynComSize=table2array(Table(:,2));
DesignMatrix=table2array(Table(:,4:17));

[~,SortSize]=sort(SynComSize);

Table(SortSize,:);

[UniqueDesign,ia,ic]=unique(DesignMatrix(SortSize,:),'rows','stable');

UniqueDesignOrdered=UniqueDesign(:,[3,8,6,12,10,4,9,5,2,13,14,11,7,1]);



for i=1:7
    NumberSynComm(i)=sum(sum(UniqueDesignOrdered')==i*2);
end


NumberSynCommSum=cumsum(NumberSynComm);

%% Figure 1B

BoxColor(1,:)=[1,1,1];
BoxColor(2,:)=[0.3,0.3,0.3];

figure('Position',[100,100,1200,800],'Color',[1,1,1])
hold on
imagesc(UniqueDesignOrdered');   
colormap(BoxColor)
for i=1:13
    plot([0 68],[i+0.5,i+0.5],'-','Color',[0.6 0.6 0.6])
end

for i=1:67
    plot([i+0.5,i+0.5],[0 15],'-','Color',[0.6 0.6 0.6])
end
for i=1:6
    plot([NumberSynCommSum(i)+0.5 NumberSynCommSum(i)+0.5],[0,15],'-','LineWidth',1,'Color',[0.1 0.1 0.1])
end
box on
axis([0.5 67.5 0.5 14.5])
yticks([])
yticklabels('')
xticks([])
xlabel('Community design')
set(gca,'FontSize',20,'LineWidth',1.5)
set(gca,'layer','top')
set(gcf,'PaperType','A3')
print('Figure1B.pdf','-dpdf','-r300')

%% Figure 1D

BoxColor2=[0.9:(-0.75/13):0.15;1:(-0.45/13):0.55;0.8:(-0.75/13):0.05]';

TableForLasso=[];
for i=1:67
    TableSele=table2array(Table(sum(table2array(Table(:,4:17))==UniqueDesign(i,:),2)==14,:));
    TableForLasso(i,1:17)=TableSele(1,1:17);
    TableForLasso(i,18)=mean(TableSele(:,18));
    TableForLasso(i,19)=std(TableSele(:,18));
end

figure('Position',[100,100,1200,800],'Color',[1,1,1])
hold on
for i=1:5
    MeasuredCFUPlot=TableForLasso(TableForLasso(:,2)==2*i & ~isnan(TableForLasso(:,18)),18);
    MeasuredCFUSTDPlot=TableForLasso(TableForLasso(:,2)==2*i & ~isnan(TableForLasso(:,18)),19);

    boxchart(ones(size(MeasuredCFUPlot))*2*i,MeasuredCFUPlot,'BoxFaceColor',BoxColor2(i*2,:),'MarkerColor',BoxColor2(i*2,:),'MarkerStyle','none','LineWidth',1.5,'BoxEdgeColor',[0,0,0],'WhiskerLineColor',[0,0,0],'LineWidth',1.5,'BoxWidth',0.6,'BoxFaceAlpha',1);
    errorbar(ones(size(MeasuredCFUPlot)).*2.*i+(rand(size(MeasuredCFUPlot))-0.5).*0.5,MeasuredCFUPlot,MeasuredCFUSTDPlot,'s','Color',BoxColor2(i*2,:).*0.6,'LineWidth',1.5,'MarkerSize',8);
end
i=7;
MeasuredCFUPlot=TableForLasso(TableForLasso(:,2)==2*i & ~isnan(TableForLasso(:,19)),18);
MeasuredCFUSTDPlot=TableForLasso(TableForLasso(:,2)==2*i & ~isnan(TableForLasso(:,18)),19);

errorbar(ones(size(MeasuredCFUPlot)).*2.*i+(rand(size(MeasuredCFUPlot))-0.5).*0.5,MeasuredCFUPlot,MeasuredCFUSTDPlot,'s','Color',BoxColor2(i*2,:).*0.6,'LineWidth',1.5,'MarkerSize',8);

axis([0 15 -5 1])
set(gca,'FontSize',24,'LineWidth',1.5)
set(gca,'layer','top')
xticks(2:2:14)
xlabel('Number of strains in SynComs')
ylabel('Change in Pst DC3000 Load (log 10)')
box on
set(gcf,'PaperType','A3')
print('Figure1D.pdf','-dpdf','-r300')

save('Figure1.mat')

