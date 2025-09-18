clear

StrainListTable=readtable('Single Strain CFU.xlsx','Sheet',5);
StrainList14=table2cell(StrainListTable(:,2));
StrainGenus14=table2cell(StrainListTable(:,3));

StrainListFull=['Root401';StrainList14];
StrainGenusFull=['Root401';StrainGenus14];

[GroupName,~,~]=unique(table2cell(StrainListTable(:,4)));
GroupName=GroupName([1,3,4,2]);
GroupNum=[3,2,3,6];


Test1Table=readtable('Root401 Single Strain CFU Data1.xlsx','Sheet',1);
Test2Table=readtable('Root401 Single Strain CFU Data2.xlsx','Sheet',1);

Root401CFU1=table2array(Test1Table(CellCompare(table2cell(Test1Table(:,7)),'Root401K')>0,6));
for i=1:15
    CFU{i,1}=table2array(Test1Table(CellCompare(table2cell(Test1Table(:,7)),StrainListFull{i})>0,6))-mean(Root401CFU1);
end

StrainListFull=['Root401';StrainList14];
Root401CFU2=table2array(Test2Table(CellCompare(table2cell(Test2Table(:,8)),'Root401')>0,7));

for i=1:15
    CFU{i,2}=table2array(Test2Table(CellCompare(table2cell(Test2Table(:,8)),StrainListFull{i})>0,7))-mean(Root401CFU2);
end

MeanCFU=cellfun(@mean,CFU);
MeanCFUPre=mean(MeanCFU(1:15,1),2,'omitnan');

MeanCFUPre(1)=0;


BoxColorGroup=GroupColorGenerator(GroupNum);
BoxColor=ColorGenerator(15);
BoxColor(2:15,:)=BoxColorGroup(table2array(StrainListTable(:,5))-1,:);
BoxColor=BoxColor-0.02;
BoxColor(BoxColor<0)=0;
BoxColor(1,:)=[0.5 0.5 0.5];

[MeanCFUSort,MeanCFUID]=sort(MeanCFUPre);

if abs(MeanCFUPre(1))<10^-6
    MeanCFUPre(1)=[];
end

SelectStrainNum=[2,4,6,8,10];

for i=1:5

    StrainCombination=table2array(readtable('SelectStrain.xlsx','Sheet',i));
    SyncomDataTable=readtable('Root401 SyncomData.xlsx','Sheet',i);

    StrainNum=SelectStrainNum(i);
    MeasuredCFUAll{i}=[];
    for j=1:size(StrainCombination,1)
        LassoOrder1{i}(j,:)=ones(1,14).*-1;
        LassoOrder1{i}(j,StrainCombination(j,1:StrainNum))=1;
        Order2Matrix=nchoosek(LassoOrder1{i}(j,:),2);
        LassoOrder2{i}(j,:)=[LassoOrder1{i}(j,:),[Order2Matrix(:,1).*Order2Matrix(:,2)]'];      
        MeasuredCFU{i}(j)=mean(table2array(SyncomDataTable(table2array(SyncomDataTable(:,1))==j,4))-mean(table2array(SyncomDataTable(table2array(SyncomDataTable(:,1))==0,4))));
        MeasuredCFUSTD{i}(j)=std(table2array(SyncomDataTable(table2array(SyncomDataTable(:,1))==j,4))-mean(table2array(SyncomDataTable(table2array(SyncomDataTable(:,1))==0,4))));
        MeasuredCFUAll{i}=[MeasuredCFUAll{i};table2array(SyncomDataTable(table2array(SyncomDataTable(:,1))==j,4))-mean(table2array(SyncomDataTable(table2array(SyncomDataTable(:,1))==0,4)))];
    end
end

Syn14CFU=mean(table2array(SyncomDataTable(table2array(SyncomDataTable(:,1))==1,4))-mean(table2array(SyncomDataTable(table2array(SyncomDataTable(:,1))==0,4))));
Syn14CFUAll=table2array(SyncomDataTable(table2array(SyncomDataTable(:,1))==1,4))-mean(table2array(SyncomDataTable(table2array(SyncomDataTable(:,1))==0,4)));

SingleAll=MeanCFU(:);

SingleAll(abs(SingleAll)<10^-6)=[];
SingleAll(isnan(SingleAll))=[];

BoxColor2=[0.9:(-0.65/13):0.25;0.8:(-0.65/13):0.15;1:(-0.35/13):0.65]';


%% Figure 4A
figure('Position',[100,100,1200,800],'Color',[1,1,1])
hold on

for i=1:5
    boxchart(ones(size(MeasuredCFU{i}))*2*i,MeasuredCFU{i},'BoxFaceColor',BoxColor2(i*2,:),'MarkerColor',BoxColor2(i*2,:),'MarkerSize',10,'LineWidth',1.5,'BoxEdgeColor',[0,0,0],'WhiskerLineColor',[0,0,0],'LineWidth',1.5,'BoxWidth',0.6,'BoxFaceAlpha',1);
    errorbar(ones(size(MeasuredCFU{i})).*2.*i+(rand(size(MeasuredCFU{i}))-0.5).*0.5,MeasuredCFU{i},MeasuredCFUSTD{i},'s','Color',BoxColor2(i*2,:).*0.6,'LineWidth',1.5,'MarkerSize',8);
end

% boxchart([14,14,14],table2array(SyncomDataTable(table2array(SyncomDataTable(:,1))==1,4))-mean(table2array(SyncomDataTable(table2array(SyncomDataTable(:,1))==0,4))),'BoxFaceColor',BoxColor2(14,:),'BoxEdgeColor',[0,0,0],'WhiskerLineColor',[0,0,0],'LineWidth',1.5,'BoxWidth',0.6,'BoxFaceAlpha',1);
errorbar(14,mean(Syn14CFUAll),std(Syn14CFUAll),'s','Color',BoxColor2(i*2,:).*0.6,'LineWidth',1.5,'MarkerSize',8);

axis([0 15 -3.5 0.5])
set(gca,'FontSize',24,'LineWidth',1.5)
set(gca,'layer','top')
xticks(2:2:14)
xlabel('SynCom Strain Number')
ylabel('Root 401 CFU Decrease (log 10)')
box on

set(gcf,'PaperType','A3')

print('Figure4A.pdf','-dpdf','-r300')

LassoOrder1X=[cell2mat(LassoOrder1');ones(1,14)];
LassoOrder2X=[cell2mat(LassoOrder2');ones(1,105)];

LassoOrder1y=[cell2mat(MeasuredCFU)';mean(Syn14CFU)];


[LOOB1,LOOFitInfo1]=lasso(LassoOrder1X, LassoOrder1y,'cv',size(LassoOrder1y,1));
[LOOB2,LOOFitInfo2]=lasso(LassoOrder2X, LassoOrder1y,'cv',size(LassoOrder1y,1));


LOOLA2Result=sum(LassoOrder2X.*LOOB2(:,LOOFitInfo2.Index1SE)',2)+LOOFitInfo2.Intercept(LOOFitInfo2.Index1SE);
LOOR2_Test_LA2=fitlm(LOOLA2Result,LassoOrder1y).Rsquared.Ordinary;

SynCommSize=[];
for i=1:5
    SynCommSize(end+1:end+size(LassoOrder1{i},1))=i*2;
end
SynCommSize(end+1)=14;


%% Figure 4B

figure('Position',[100,100,1200,1000],'Color',[1,1,1])
hold on
for i=1:5
    errorbar(LOOLA2Result(SynCommSize==2*i),MeasuredCFU{i},MeasuredCFUSTD{i},'o','Color',BoxColor2(i*2,:),'MarkerSize',10,'LineWidth',1.5)
end
i=7;
errorbar(LOOLA2Result(SynCommSize==2*i),mean(Syn14CFU),std(Syn14CFU),'o','Color',BoxColor2(i*2,:),'MarkerSize',10,'LineWidth',1.5)


plot([-4.5 1],[-4.5 1],'--','Color',[0.5 0.5 0.5])

box on
axis equal
axis([-3.6 0.6 -3.6 0.6])
set(gca,'FontSize',28,'LineWidth',1.5)
set(gca,'layer','top')
xlabel('Predicted by Regression (log 10)')
ylabel('Measured CFU Decrease (log 10)')
legend({'2 Strain','4 Strain','6 Strain','8 Strain','10 Strain','14 Strain'},'Location','SouthEast','FontSize',22)


set(gcf,'PaperType','A3')

print('Figure4B.pdf','-dpdf','-r300')


IndexOrder=nchoosek(1:14,2);
Order2CorrIndex=zeros(14);
Order2CorrValue=zeros(14);

for i=1:100
    [B1,FitInfo1]=lasso(LassoOrder1X, LassoOrder1y,'cv',10);
    [CorrLa1All(i),PLaAll1(i)]=corr(sum(LassoOrder1X.*B1(:,FitInfo1.Index1SE)',2)+FitInfo1.Intercept(FitInfo1.Index1SE),LassoOrder1y);

    [B2,FitInfo2]=lasso(LassoOrder2X, LassoOrder1y,'cv',10);
    [CorrLa2All(i),PLaAll2(i)]=corr(sum(LassoOrder2X.*B2(:,FitInfo2.Index1SE)',2)+FitInfo2.Intercept(FitInfo2.Index1SE),LassoOrder1y);

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


%% Figure 4C


[~,RegSortID]=sort(mean(RegressionParametersL1));

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

axis([0 15 -0.5 0.1])
set(gca,'FontSize',24,'LineWidth',1.5)
set(gca,'layer','top')
xticks(1:14)
xticklabels(StrainGenusFull(RegSortID+1))
box on
legend('LASSO regression model (order=1)','LASSO regression model (order=2)', 'Location','southeast')

set(gcf,'PaperType','A3')
print('Figure4C.pdf','-dpdf','-r300')


%% Figure 4D

LassoOrderTest1Blank=ones(1,14).*-1;
LassoOrderTest2Blank=[ones(1,14).*-1,ones(1,91).*1];

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


figure('Position',[100,100,1050,1000],'Color',[1,1,1])
hold on
imagesc(Delta,[-0.5 0.5])
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
print('Figure4D.pdf','-dpdf')

save('Figure4.mat')