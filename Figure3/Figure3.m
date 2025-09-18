clear
clc
load('Figure2.mat')

% Single Strain Data

StrainListTable=readtable('Single Strain CFU.xlsx','Sheet',5);
StrainList14=table2cell(StrainListTable(:,2));
StrainGenus14=table2cell(StrainListTable(:,3));
StrainListFull=['PstDC3000';StrainList14];
StrainGenusFull=['PstDC3000';StrainGenus14];

[GroupName,~,~]=unique(table2cell(StrainListTable(:,4)));
GroupName=GroupName([1,3,4,2]);
GroupNum=[3,2,3,6];


Test1Table=readtable('Single Strain CFU.xlsx','Sheet',1);
Test2Table=readtable('Single Strain CFU.xlsx','Sheet',2);
Test3Table=readtable('Single Strain CFU.xlsx','Sheet',3);
Test4Table=readtable('Single Strain CFU.xlsx','Sheet',4);

DCCFU1=table2array(Test1Table(CellCompare(table2cell(Test1Table(:,2)),'PstDC3000')>0,4));
DCCFU2=table2array(Test2Table(CellCompare(table2cell(Test2Table(:,2)),'PstDC3000')>0,4));
DCCFU3=table2array(Test3Table(CellCompare(table2cell(Test3Table(:,2)),'PstDC3000')>0,4));
DCCFU4=table2array(Test4Table(CellCompare(table2cell(Test4Table(:,2)),'PstDC3000')>0,4));

DCFW1=table2array(Test1Table(CellCompare(table2cell(Test1Table(:,2)),'PstDC3000')>0,5));
DCFW2=table2array(Test2Table(CellCompare(table2cell(Test2Table(:,2)),'PstDC3000')>0,5));
DCFW3=table2array(Test3Table(CellCompare(table2cell(Test3Table(:,2)),'PstDC3000')>0,5));


for i=1:15
    CFU{i,1}=table2array(Test1Table(CellCompare(table2cell(Test1Table(:,2)),StrainListFull{i})>0,4))-mean(DCCFU1);
    CFU{i,2}=table2array(Test2Table(CellCompare(table2cell(Test2Table(:,2)),StrainListFull{i})>0,4))-mean(DCCFU2);
    CFU{i,3}=table2array(Test3Table(CellCompare(table2cell(Test3Table(:,2)),StrainListFull{i})>0,4))-mean(DCCFU3);
    CFU{i,4}=table2array(Test4Table(CellCompare(table2cell(Test4Table(:,2)),StrainListFull{i})>0,4))-mean(DCCFU4);


    AbCFU{i,1}=table2array(Test1Table(CellCompare(table2cell(Test1Table(:,2)),StrainListFull{i})>0,4));
    AbCFU{i,2}=table2array(Test2Table(CellCompare(table2cell(Test2Table(:,2)),StrainListFull{i})>0,4));
    AbCFU{i,3}=table2array(Test3Table(CellCompare(table2cell(Test3Table(:,2)),StrainListFull{i})>0,4));


    FW{i,1}=table2array(Test1Table(CellCompare(table2cell(Test1Table(:,2)),StrainListFull{i})>0,5));
    FW{i,2}=table2array(Test2Table(CellCompare(table2cell(Test2Table(:,2)),StrainListFull{i})>0,5));
    FW{i,3}=table2array(Test3Table(CellCompare(table2cell(Test3Table(:,2)),StrainListFull{i})>0,5));

end

MeanCFU=cellfun(@mean,CFU);
MeanCFUPre=mean(MeanCFU(1:15,:),2,'omitnan');




% LASSO CV-10
IndexOrder=nchoosek(1:14,2);
Order2CorrIndex=zeros(14);
Order2CorrValue=zeros(14);

for i=1:100

    [B1,FitInfo1]=lasso(LassoOrder1X, LassoOrder1y,'cv',10);
    [B2,FitInfo2]=lasso(LassoOrder2X, LassoOrder1y,'cv',10);
    
    LA1Result(i,:)=sum(LassoOrder1X.*B1(:,FitInfo1.Index1SE)',2)+FitInfo1.Intercept(FitInfo1.Index1SE);
    LA2Result(i,:)=sum(LassoOrder2X.*B2(:,FitInfo2.Index1SE)',2)+FitInfo2.Intercept(FitInfo2.Index1SE);
    
    LA1ParaNum(i)=sum(B1(:,FitInfo1.Index1SE)~=0);
    LA2ParaNum(i)=sum(B2(:,FitInfo2.Index1SE)~=0);
    
    R2_Test_LA1(i)=fitlm(LA1Result(i,:),LassoOrder1y).Rsquared.Ordinary;
    R2_Test_LA2(i)=fitlm(LA2Result(i,:),LassoOrder1y).Rsquared.Ordinary;

    RegressionParametersL1(i,:)=B1(:,FitInfo1.Index1SE)';
    RegressionParametersL2(i,:)=B2(:,FitInfo2.Index1SE)';

    RegressionInterceptL1(i,:)=FitInfo1.Intercept(FitInfo1.Index1SE);
    RegressionInterceptL2(i,:)=FitInfo2.Intercept(FitInfo2.Index1SE);

    NonZeroOrder2Index=IndexOrder(B2(15:end,FitInfo2.Index1SE)'~=0,:);
    NonZeroOrder2Index2=find(B2(15:end,FitInfo2.Index1SE)'~=0);

    for j=1:size(NonZeroOrder2Index,1)
        Order2CorrIndex(NonZeroOrder2Index(j,1),NonZeroOrder2Index(j,2))=Order2CorrIndex(NonZeroOrder2Index(j,1),NonZeroOrder2Index(j,2))+1;
        Order2CorrValue(NonZeroOrder2Index(j,1),NonZeroOrder2Index(j,2),i)=B2(14+NonZeroOrder2Index2(j),FitInfo2.Index1SE);
    end
end

[~,RegSortID]=sort(mean(RegressionParametersL1));

% Figure 3A

figure('Position',[100,100,1000,1000],'Color',[1,1,1])
hold on
x = 1:length(RegSortID);

L1xgroupdata = repelem(x-0.15, 100)';
L2xgroupdata = repelem(x+0.15, 100)';
L1ydata = RegressionParametersL1(:,RegSortID);
L2ydata = RegressionParametersL2(:,RegSortID);

boxchart(L1xgroupdata,L1ydata(:),'WhiskerLineColor',[0,0,0],'LineWidth',1.5,'BoxWidth',0.6,'BoxFaceAlpha',0.6,'boxwidth',0.2)
boxchart(L2xgroupdata,L2ydata(:),'WhiskerLineColor',[0,0,0],'LineWidth',1.5,'BoxWidth',0.6,'BoxFaceAlpha',0.6,'boxwidth',0.2)

% plot(x-0.15,LOORegressionParametersL1(RegSortID),'^','MarkerSize',10,'LineWidth',1.5,'Color',[0.0660 0.4430 0.7450])
% plot(x+0.15,LOORegressionParametersL2(RegSortID),'^','MarkerSize',10,'LineWidth',1.5,'Color',[0.8660 0.3290 0])

axis([0 15 -0.9 0.1])
set(gca,'FontSize',24,'LineWidth',1.5)
set(gca,'layer','top')
xticks(1:14)
xticklabels(StrainGenusFull(RegSortID+1))
box on
legend('LASSO regression model (order=1)','LASSO regression model (order=2)', 'Location','southeast')

set(gcf,'PaperType','A3')
print('Figure3A.pdf','-dpdf','-r300')


[MeanCFUSort,MeanCFUID]=sort(MeanCFUPre);


BoxColorGroup=GroupColorGenerator(GroupNum);
BoxColor=ColorGenerator(15);
BoxColor(2:15,:)=BoxColorGroup(table2array(StrainListTable(:,5))-1,:);
BoxColor=BoxColor-0.02;
BoxColor(BoxColor<0)=0;
BoxColor(1,:)=[0.5 0.5 0.5];


% Figure 3B

figure('Position',[100,100,1000,1000],'Color',[1,1,1])
hold on
plot(-1,-1,'o','MarkerEdgeColor',BoxColorGroup(2,:),'MarkerFaceColor',BoxColorGroup(2,:),'LineWidth',1.5)
plot(-1,-1,'o','MarkerEdgeColor',BoxColorGroup(7,:),'MarkerFaceColor',BoxColorGroup(7,:),'LineWidth',1.5)
plot(-1,-1,'o','MarkerEdgeColor',BoxColorGroup(5,:),'MarkerFaceColor',BoxColorGroup(5,:),'LineWidth',1.5)
plot(-1,-1,'o','MarkerEdgeColor',BoxColorGroup(11,:),'MarkerFaceColor',BoxColorGroup(11,:),'LineWidth',1.5)

for i=1:7
    fill([i*2-0.5,i*2-0.5,i*2+0.5,i*2+0.5],[-5 2 2 -5],[0.85 0.85 0.85],'LineStyle','none')
end
plot([0.5,15.5],[0,0],'--','Color',[0.3 0.3 0.3],'LineWidth',1)
% for i=1:15
%     for j=1:4
%     boxchart(ones(size(CFU{MeanCFUID(i),j}))*(i-0.2*(j-2)),CFU{MeanCFUID(i),j},'BoxFaceColor',BoxColor(MeanCFUID(i),:),'WhiskerLineColor',BoxColor(MeanCFUID(i),:),'LineWidth',1.5,'BoxWidth',0.16,'BoxFaceAlpha',0.5);
%     end
% end
for i=1:15
    boxchart(ones(size(cell2mat(CFU(MeanCFUID(i),:)'))).*i,cell2mat(CFU(MeanCFUID(i),:)'),'BoxFaceColor',BoxColor(MeanCFUID(i),:),'WhiskerLineColor',BoxColor(MeanCFUID(i),:),'LineWidth',1.5,'BoxWidth',0.5,'BoxFaceAlpha',0.5,'MarkerStyle','none');
    
    scatter(ones(size(CFU{MeanCFUID(i),1})).*i+(rand(size(CFU{MeanCFUID(i),1}))-0.5).*0.5,CFU{MeanCFUID(i),1},30,'o','MarkerEdgeColor',BoxColor(MeanCFUID(i),:),'MarkerFaceColor',BoxColor(MeanCFUID(i),:),'LineWidth',1.5,'MarkerFaceAlpha',0.5);
    scatter(ones(size(CFU{MeanCFUID(i),2})).*i+(rand(size(CFU{MeanCFUID(i),2}))-0.5).*0.5,CFU{MeanCFUID(i),2},30,'s','MarkerEdgeColor',BoxColor(MeanCFUID(i),:),'MarkerFaceColor',BoxColor(MeanCFUID(i),:),'LineWidth',1.5,'MarkerFaceAlpha',0.5);
    scatter(ones(size(CFU{MeanCFUID(i),3})).*i+(rand(size(CFU{MeanCFUID(i),3}))-0.5).*0.5,CFU{MeanCFUID(i),3},30,'d','MarkerEdgeColor',BoxColor(MeanCFUID(i),:),'MarkerFaceColor',BoxColor(MeanCFUID(i),:),'LineWidth',1.5,'MarkerFaceAlpha',0.5);
    scatter(ones(size(CFU{MeanCFUID(i),4})).*i+(rand(size(CFU{MeanCFUID(i),4}))-0.5).*0.5,CFU{MeanCFUID(i),4},30,'<','MarkerEdgeColor',BoxColor(MeanCFUID(i),:),'MarkerFaceColor',BoxColor(MeanCFUID(i),:),'LineWidth',1.5,'MarkerFaceAlpha',0.5);
end

axis([0.5 15.5 -4.5 1.5])
set(gca,'FontSize',24,'LineWidth',1.5)
set(gca,'layer','top')
xticks(1:15)
xticklabels(StrainGenusFull(MeanCFUID))
box on
ylabel('Change in Pst DC3000 load (log 10)')
legend(GroupName,'Location','SouthEast','FontSize',18)
set(gcf,'PaperType','A3')
print('Figure3B.pdf','-dpdf','-r300')



LassoOrderTest1=cell(0);
LassoOrderTest2=cell(0);

for i=1:2
    StrainNum=i;
    StrainCombination=nchoosek(1:14,StrainNum);

    for j=1:size(StrainCombination,1)
        LassoOrderTest1{i}(j,:)=ones(1,14).*-1;
        LassoOrderTest1{i}(j,StrainCombination(j,1:StrainNum))=1;
        OrderTest2Matrix=nchoosek(LassoOrderTest1{i}(j,:),2);
        LassoOrderTest2{i}(j,:)=[LassoOrderTest1{i}(j,:),[OrderTest2Matrix(:,1).*OrderTest2Matrix(:,2)]'];
    end
end

LassoOrderTest1Blank=ones(1,14).*-1;
LassoOrderTest2Blank=[ones(1,14).*-1,ones(1,91).*1];

LOORegressionParametersL1=LOOB1(:,LOOFitInfo1.Index1SE)';
LOORegressionParametersL2=LOOB2(:,LOOFitInfo2.Index1SE)';


MeanCFUPreLa1=sum(LassoOrderTest1{1}.*LOORegressionParametersL1,2)+LOOFitInfo1.Intercept(LOOFitInfo1.Index1SE);
MeanCFUPreLa2=sum(LassoOrderTest2{1}.*LOORegressionParametersL2,2)+LOOFitInfo2.Intercept(LOOFitInfo2.Index1SE);


R2MeanCFULa1=fitlm(mean(MeanCFUPreLa1,2),mean(MeanCFU(2:end,:),2,'omitnan')).Rsquared.Ordinary;
R2MeanCFULa2=fitlm(mean(MeanCFUPreLa2,2),mean(MeanCFU(2:end,:),2,'omitnan')).Rsquared.Ordinary;


%% Figure 3C

figure('Position',[100,100,1200,1000],'Color',[1,1,1])
hold on
for i=1:14
    errorbar(mean(MeanCFUPreLa1(i,:),2),mean(MeanCFU(i+1,:),2,'omitnan'),std(MeanCFU(i+1,:),[],2,'omitnan'),std(MeanCFU(i+1,:),[],2,'omitnan'),std(MeanCFUPreLa1(i,:),[],2),std(MeanCFUPreLa1(i,:),[],2),'o','Color',BoxColor(i+1,:),'LineWidth',1.5,'MarkerSize',12)
end

plot([-7.5 0.5],[-7.5 0.5],'--','Color',[0.5 0.5 0.5])

box on
axis equal
axis([-4.5 0.5 -4.5 0.5])
set(gca,'FontSize',28,'LineWidth',1.5)
set(gca,'layer','top')
xlabel('Predicted by Regression (log 10)')
ylabel('Measured CFU Decrease (log 10)')
legend(StrainListFull(2:15),'Location','SouthEast','FontSize',18)

text(-4,0.25,['R^2 = ',num2str(R2MeanCFULa1,3)],'FontSize',24)
set(gcf,'PaperType','A3')
print('Figure3C_Order1.pdf','-dpdf','-r300')


figure('Position',[100,100,1200,1000],'Color',[1,1,1])

hold on

for i=1:14
    errorbar(mean(MeanCFUPreLa2(i,:),2),mean(MeanCFU(i+1,:),2,'omitnan'),std(MeanCFU(i+1,:),[],2,'omitnan'),std(MeanCFU(i+1,:),[],2,'omitnan'),std(MeanCFUPreLa2(i,:),[],2),std(MeanCFUPreLa2(i,:),[],2),'o','Color',BoxColor(i+1,:),'LineWidth',1.5,'MarkerSize',12)
end

plot([-7.5 0.5],[-7.5 0.5],'--','Color',[0.5 0.5 0.5])

box on
axis equal
axis([-4.5 0.5 -4.5 0.5])
set(gca,'FontSize',28,'LineWidth',1.5)
set(gca,'layer','top')
xlabel('Predicted by Regression (log 10)')
ylabel('Measured CFU Decrease (log 10)')
legend(StrainListFull(2:15),'Location','SouthEast','FontSize',18)

text(-4,0.25,['R^2 = ',num2str(R2MeanCFULa2,3)],'FontSize',24)

set(gcf,'PaperType','A3')
print('Figure3C_Order2.pdf','-dpdf','-r300')

LassoOrderTest1=cell(0);
LassoOrderTest2=cell(0);

for i=1:2
    StrainNum=i;
    StrainCombination=nchoosek(1:14,StrainNum);

    for j=1:size(StrainCombination,1)
        LassoOrderTest1{i}(j,:)=ones(1,14).*-1;
        LassoOrderTest1{i}(j,StrainCombination(j,1:StrainNum))=1;
        OrderTest2Matrix=nchoosek(LassoOrderTest1{i}(j,:),2);
        LassoOrderTest2{i}(j,:)=[LassoOrderTest1{i}(j,:),[OrderTest2Matrix(:,1).*OrderTest2Matrix(:,2)]'];
    end
end

LassoOrderTest1Blank=ones(1,14).*-1;
LassoOrderTest2Blank=[ones(1,14).*-1,ones(1,91).*1];

LOORegressionParametersL1=LOOB1(:,LOOFitInfo1.Index1SE)';
LOORegressionParametersL2=LOOB2(:,LOOFitInfo2.Index1SE)';


MeanCFUPreLa1=sum(LassoOrderTest1{1}.*LOORegressionParametersL1,2)+LOOFitInfo1.Intercept(LOOFitInfo1.Index1SE);
MeanCFUPreLa2=sum(LassoOrderTest2{1}.*LOORegressionParametersL2,2)+LOOFitInfo2.Intercept(LOOFitInfo2.Index1SE);


R2MeanCFULa1=fitlm(mean(MeanCFUPreLa1,2),mean(MeanCFU(2:end,:),2,'omitnan')).Rsquared.Ordinary;
R2MeanCFULa2=fitlm(mean(MeanCFUPreLa2,2),mean(MeanCFU(2:end,:),2,'omitnan')).Rsquared.Ordinary;


PairInterLa2=sum(LassoOrderTest2{2}.*LOORegressionParametersL2,2)+LOOFitInfo2.Intercept(LOOFitInfo2.Index1SE);

PairZero=sum(LassoOrderTest2Blank.*LOORegressionParametersL2,2)+LOOFitInfo2.Intercept(LOOFitInfo2.Index1SE);


k=0;
Delta=zeros(14);

for i=1:13
    for j=(i+1):14
        k=k+1;
        Delta(i,j)=PairInterLa2(k,:)-(MeanCFUPreLa2(i,:)+MeanCFUPreLa2(j,:)-PairZero);
        kIndex(i,j)=k;
    end
end

cmap=[];
cmap(1:128,1)=0.85:0.15/127:1;
cmap(1:128,2)=0.23:0.77/127:1;
cmap(1:128,3)=0.10:0.90/127:1;
cmap(129:256,1)=1:-1.00/127:0.00;
cmap(129:256,2)=1:-0.65/127:0.35;
cmap(129:256,3)=1:-0.15/127:0.85;

%% Figure 3D

figure('Position',[100,100,1050,1000],'Color',[1,1,1])
hold on
imagesc(Delta,[-1 1])
% h=heatmap(StrainGenus14,StrainGenus14,mean(Order2CorrValue,3),'CellLabelFormat','%0.00d');
colorbar off
axis equal
axis([0.5 14.5 0.5 14.5])
colormap(cmap)
xticks(1:14)
yticks(1:14)
xticklabels(StrainGenus14)
yticklabels(StrainGenus14)
box on
set(gca,'YDir','reverse')
set(gca,'yaxislocation','right')
set(gca,'xaxislocation','top')

for i=1:13
plot([i+0.5,i+0.5],[0.5,i+1.5],'-k','LineWidth',1)
plot([i-0.5,14.5],[i+0.5,i+0.5],'-k','LineWidth',1)
end
plot([0.5,0.5],[0.5,1.5],'-k','LineWidth',1)
plot([13.5,14.5],[14.5,14.5],'-k','LineWidth',1)

set(gca,'FontSize',16,'LineWidth',1.5)
box off
set(gcf,'PaperType','A3')
set(gcf,'Renderer','painters')
set(gca,'layer','top')
colorbar
print('Figure3D.pdf','-dpdf')


%% Fiugre 3E is move to SI
% Prediction
k=0;
Delta=zeros(14);

for i=1:13
    for j=(i+1):14
        k=k+1;
        Delta(i,j)=PairInterLa2(k,:)-(MeanCFUPreLa2(i,:)+MeanCFUPreLa2(j,:)-PairZero);
        kIndex(i,j)=k;
    end
end
PairMean=PairInterLa2;
SingleMean=MeanCFUPreLa2;


i=11;
j=12;

figure('Position',[100,100,800,600],'Color',[1,1,1])
hold on
plot(1.25,SingleMean(i),'o','MarkerSize',10,'LineWidth',1.5,'Color',BoxColor(i+1,:))
plot(1.5,SingleMean(j),'o','MarkerSize',10,'LineWidth',1.5,'Color',BoxColor(j+1,:))
plot(2.5,SingleMean(i)+SingleMean(j)-mean(PairZero),'o','MarkerSize',10,'LineWidth',1.5,'Color',[0.85,0.325,0.098])
plot([0.25 1.25],[mean(PairZero) SingleMean(i)],':','LineWidth',1,'color',[0.5,0.5,0.5])
plot([0.25 1.5],[mean(PairZero) SingleMean(j)],':','LineWidth',1,'color',[0.5,0.5,0.5])
plot([1.25 2.5],[SingleMean(i),SingleMean(i)+SingleMean(j)-mean(PairZero)],':','LineWidth',1,'color',[0.5,0.5,0.5])
plot([1.5 2.5],[SingleMean(j),SingleMean(i)+SingleMean(j)-mean(PairZero)],':','LineWidth',1,'color',[0.5,0.5,0.5])

plot(2.5,PairMean(kIndex(i,j)),'s','MarkerSize',10,'LineWidth',1.5,'Color',[0 0.45 0.74])

plot([2.5 2.65],[SingleMean(i)+SingleMean(j)-mean(PairZero),SingleMean(i)+SingleMean(j)-mean(PairZero)],':','LineWidth',1,'color',[0.5,0.5,0.5])
plot([2.5 2.65],[PairMean(kIndex(i,j)),PairMean(kIndex(i,j))],':','LineWidth',1,'color',[0.5,0.5,0.5])

text(2.7,mean([PairMean(kIndex(i,j)),(SingleMean(i)+SingleMean(j))]),['\delta=',num2str(Delta(i,j))],'FontSize',18)

axis([0 3 -6 0.5])
box on
set(gca,'FontSize',24,'LineWidth',1.5)

xticks([])
xlabel([StrainGenusFull{i+1}, ' v.s. ',StrainGenusFull{j+1}])
print('Figure3E_Pred_1.pdf','-dpdf','-r300')


i=7;
j=13;

figure('Position',[100,100,800,600],'Color',[1,1,1])
hold on
plot(1.25,SingleMean(i),'o','MarkerSize',10,'LineWidth',1.5,'Color',BoxColor(i+1,:))
plot(1.5,SingleMean(j),'o','MarkerSize',10,'LineWidth',1.5,'Color',BoxColor(j+1,:))
plot(2.5,SingleMean(i)+SingleMean(j)-mean(PairZero),'o','MarkerSize',10,'LineWidth',1.5,'Color',[0.85,0.325,0.098])
plot([0.25 1.25],[mean(PairZero) SingleMean(i)],':','LineWidth',1,'color',[0.5,0.5,0.5])
plot([0.25 1.5],[mean(PairZero) SingleMean(j)],':','LineWidth',1,'color',[0.5,0.5,0.5])
plot([1.25 2.5],[SingleMean(i),SingleMean(i)+SingleMean(j)-mean(PairZero)],':','LineWidth',1,'color',[0.5,0.5,0.5])
plot([1.5 2.5],[SingleMean(j),SingleMean(i)+SingleMean(j)-mean(PairZero)],':','LineWidth',1,'color',[0.5,0.5,0.5])

plot(2.5,PairMean(kIndex(i,j)),'s','MarkerSize',10,'LineWidth',1.5,'Color',[0 0.45 0.74])

plot([2.5 2.65],[SingleMean(i)+SingleMean(j)-mean(PairZero),SingleMean(i)+SingleMean(j)-mean(PairZero)],':','LineWidth',1,'color',[0.5,0.5,0.5])
plot([2.5 2.65],[PairMean(kIndex(i,j)),PairMean(kIndex(i,j))],':','LineWidth',1,'color',[0.5,0.5,0.5])

text(2.7,mean([PairMean(kIndex(i,j)),(SingleMean(i)+SingleMean(j))]),['\delta=',num2str(Delta(i,j))],'FontSize',18)

axis([0 3 -6 0.5])
box on
set(gca,'FontSize',24,'LineWidth',1.5)

xticks([])
xlabel([StrainGenusFull{i+1}, ' v.s. ',StrainGenusFull{j+1}])
print('Figure3E_Pred_2.pdf','-dpdf','-r300')


i=4;
j=14;

figure('Position',[100,100,800,600],'Color',[1,1,1])
hold on
plot(1.25,SingleMean(i),'o','MarkerSize',10,'LineWidth',1.5,'Color',BoxColor(i+1,:))
plot(1.5,SingleMean(j),'o','MarkerSize',10,'LineWidth',1.5,'Color',BoxColor(j+1,:))
plot(2.5,SingleMean(i)+SingleMean(j)-mean(PairZero),'o','MarkerSize',10,'LineWidth',1.5,'Color',[0.85,0.325,0.098])
plot([0.25 1.25],[mean(PairZero) SingleMean(i)],':','LineWidth',1,'color',[0.5,0.5,0.5])
plot([0.25 1.5],[mean(PairZero) SingleMean(j)],':','LineWidth',1,'color',[0.5,0.5,0.5])
plot([1.25 2.5],[SingleMean(i),SingleMean(i)+SingleMean(j)-mean(PairZero)],':','LineWidth',1,'color',[0.5,0.5,0.5])
plot([1.5 2.5],[SingleMean(j),SingleMean(i)+SingleMean(j)-mean(PairZero)],':','LineWidth',1,'color',[0.5,0.5,0.5])

plot(2.5,PairMean(kIndex(i,j)),'s','MarkerSize',10,'LineWidth',1.5,'Color',[0 0.45 0.74])

plot([2.5 2.65],[SingleMean(i)+SingleMean(j)-mean(PairZero),SingleMean(i)+SingleMean(j)-mean(PairZero)],':','LineWidth',1,'color',[0.5,0.5,0.5])
plot([2.5 2.65],[PairMean(kIndex(i,j)),PairMean(kIndex(i,j))],':','LineWidth',1,'color',[0.5,0.5,0.5])

text(2.7,mean([PairMean(kIndex(i,j)),(SingleMean(i)+SingleMean(j))]),['\delta=',num2str(Delta(i,j))],'FontSize',18)

axis([0 3 -6 0.5])
box on
set(gca,'FontSize',24,'LineWidth',1.5)

xticks([])
xlabel([StrainGenusFull{i+1}, ' v.s. ',StrainGenusFull{j+1}])
print('Figure3E_Pred_3.pdf','-dpdf','-r300')



% Measurement 

PairwiseIndex=TableForLasso(TableForLasso(:,2)==2,4:17);
PairwiseResult=LassoOrder1y(TableForLasso(:,2)==2);
PairwiseResultSTD=TableForLasso(TableForLasso(:,2)==2,19);
MeanCFUPre=mean(MeanCFU(1:15,:),2,'omitnan');
STDCFUPre=std(MeanCFU(1:15,:),[],2,'omitnan');

for i=1:size(PairwiseIndex,1)
    StrainIndex(i,:)=find(PairwiseIndex(i,:)==1);
end


PairwiseInteraction=[];
for i=[2,22,19]

    figure('Position',[100,100,800,600],'Color',[1,1,1])
    hold on
    errorbar(1.25,MeanCFUPre(StrainIndex(i,1)+1),STDCFUPre(StrainIndex(i,1)+1),'o','MarkerSize',10,'LineWidth',1.5,'Color',BoxColor((StrainIndex(i,1)+1),:))
    errorbar(1.5,MeanCFUPre(StrainIndex(i,2)+1),STDCFUPre(StrainIndex(i,2)+1),'o','MarkerSize',10,'LineWidth',1.5,'Color',BoxColor((StrainIndex(i,2)+1),:))
    errorbar(2.5,sum(MeanCFUPre(StrainIndex(i,:)+1)),sum(STDCFUPre(StrainIndex(i,:)+1)),'o','MarkerSize',10,'LineWidth',1.5,'Color',[0.85,0.325,0.098])
    plot([0.25 1.25],[0 MeanCFUPre(StrainIndex(i,1)+1)],':','LineWidth',1,'color',[0.5,0.5,0.5])
    plot([0.25 1.5],[0 MeanCFUPre(StrainIndex(i,2)+1)],':','LineWidth',1,'color',[0.5,0.5,0.5])
    plot([1.25 2.5],[MeanCFUPre(StrainIndex(i,1)+1),sum(MeanCFUPre(StrainIndex(i,:)+1))],':','LineWidth',1,'color',[0.5,0.5,0.5])
    plot([1.5 2.5],[MeanCFUPre(StrainIndex(i,2)+1),sum(MeanCFUPre(StrainIndex(i,:)+1))],':','LineWidth',1,'color',[0.5,0.5,0.5])

    errorbar(2.5,PairwiseResult(i),PairwiseResultSTD(i),'s','MarkerSize',10,'LineWidth',1.5,'Color',[0 0.45 0.74])

    plot([2.5 2.65],[sum(MeanCFUPre(StrainIndex(i,:)+1)),sum(MeanCFUPre(StrainIndex(i,:)+1))],':','LineWidth',1,'color',[0.5,0.5,0.5])
    plot([2.5 2.65],[PairwiseResult(i),PairwiseResult(i)],':','LineWidth',1,'color',[0.5,0.5,0.5])
    Delta(i)=abs(PairwiseResult(i)-sum(MeanCFUPre(StrainIndex(i,:)+1)));
    Delta(i)=abs(PairwiseResult(i)-sum(MeanCFUPre(StrainIndex(i,:)+1)));

    text(2.7,mean([PairwiseResult(i),sum(MeanCFUPre(StrainIndex(i,:)+1))]),['\delta=',num2str(Delta(i))],'FontSize',18)

    axis([0 3 -6 0.5])
    box on
    set(gca,'FontSize',24,'LineWidth',1.5)

    xticks([])
    xlabel([StrainGenusFull{StrainIndex(i,1)+1}, ' v.s. ',StrainGenusFull{StrainIndex(i,2)+1}])
    
    
    print(['Figure3E_Measure_',num2str(i),'.pdf'],'-dpdf','-r300')

end

save('Figure3.mat')




