%% The SI Figures for "Invasion resistance in 14-member synthetic bacterial communities is controlled by multiple interacting factors"
% The index for the figure may not exatly the same as the order in the manuscript.

%% Figrue S1
% clear
% load('Figure1.mat')

DC3000FWvsCFU=readtable('FigureSI_FW_Data.xlsx','Sheet',1);
DC3000FWvsCFU=table2array(DC3000FWvsCFU(:,2:3));
[DC3000C,DC3000p]=corr(DC3000FWvsCFU(:,2),DC3000FWvsCFU(:,1));

x=DC3000FWvsCFU(:,2);
y=DC3000FWvsCFU(:,1);

[xData, yData] = prepareCurveData( x, y );
ft = fittype( 'poly1' );
[fitresult, gof] = fit( xData, yData, ft );
[intervals] = confint(fitresult);

figure('Position',[100,100,1200,800],'Color',[1,1,1])
hold on
plot(DC3000FWvsCFU(:,2),DC3000FWvsCFU(:,1),'.','MarkerSize',10)
plot([0,7],fitresult.p1.*[0,7]+fitresult.p2,'-','LineWidth',1.5,'Color',[0.5 0.5 0.5])
box on
axis([0 7 0 50])
set(gca,'FontSize',24,'LineWidth',1.5)
set(gcf,'PaperType','A3')
set(gcf,'Renderer','painters')
set(gca,'layer','top')
yticks(0:10:50)
xlabel('CFU of DC3000')
ylabel('Fresh Weight (mg)')
text(0.2,48,['Pearson = ',num2str(DC3000C,3)],'FontSize',20)
text(0.2,45,'p-Value < 0.001','FontSize',20)

print('FigureS1.pdf','-dpdf','-r300')



%% Figure S2
% clear
% load('Figure2.mat')

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

% R^2

figure('Position',[100,100,600,1000],'Color',[1,1,1])
hold on

boxchart(ones(100,1).*1,R2_Test_LA1','LineWidth',1.5,'BoxWidth',0.6,'BoxFaceAlpha',0.6,'MarkerSize',12)
boxchart(ones(100,1).*2,R2_Test_LA2','LineWidth',1.5,'BoxWidth',0.6,'BoxFaceAlpha',0.6,'MarkerSize',12)

axis([0 3 0.5 1.03])
set(gca,'FontSize',28,'LineWidth',1.5)
set(gca,'layer','top')
xticks('')
xlabel('')
ylabel('R^2')
box on
legend({'Order = 1','Order = 2'},'Location','SouthEast','FontSize',22)

set(gcf,'PaperType','A3')
print('FigureS2A.pdf','-dpdf','-r300')

% Number of Parameters
figure('Position',[100,100,600,1000],'Color',[1,1,1])
hold on

boxchart(ones(100,1).*1,LA1ParaNum','LineWidth',1.5,'BoxWidth',0.6,'BoxFaceAlpha',0.6,'MarkerSize',12)
boxchart(ones(100,1).*2,LA2ParaNum','LineWidth',1.5,'BoxWidth',0.6,'BoxFaceAlpha',0.6,'MarkerSize',12)

axis([0 3 0 30])
set(gca,'FontSize',28,'LineWidth',1.5)
set(gca,'layer','top')
xticks('')
xlabel('')
ylabel('Number of Parameters')
box on
legend({'Order = 1','Order = 2'},'Location','SouthEast','FontSize',22)

set(gcf,'PaperType','A3')
print('FigureS2B.pdf','-dpdf','-r300')


%% Figure S3 Traning Splits
% clear
% load('Figure2.mat')


for i=1:5
    TestSetIndex=find(TableForLasso(:,2)==2.*i);
    TrainningSetIndex=1:NumAll;
    TrainningSetIndex(TestSetIndex)=[];

    FullTestSetL1=LassoOrder1X(TestSetIndex,:);
    FullTestSetL2=LassoOrder2X(TestSetIndex,:);
    FullTestSetY=LassoOrder1y(TestSetIndex,:);

    FullTrainningSetL1=LassoOrder1X(TrainningSetIndex,:);
    FullTrainningSetL2=LassoOrder2X(TrainningSetIndex,:);
    FullTrainningSetY=LassoOrder1y(TrainningSetIndex,:);

    for j=1:100
        TrainingSetIndexSele=1:size(FullTrainningSetY,1);

        [B1,FitInfo1]=lasso(FullTrainningSetL1(TrainingSetIndexSele,:), FullTrainningSetY(TrainingSetIndexSele,:),'cv',10);
        [B2,FitInfo2]=lasso(FullTrainningSetL2(TrainingSetIndexSele,:), FullTrainningSetY(TrainingSetIndexSele,:),'cv',10);

        R2TrainLa1(j,i) = fitlm(sum(FullTrainningSetL1(TrainingSetIndexSele,:).*B1(:,FitInfo1.Index1SE)',2)+FitInfo1.Intercept(FitInfo1.Index1SE),FullTrainningSetY(TrainingSetIndexSele,:)).Rsquared.Ordinary;
        R2TestLa1(j,i) = fitlm(sum(FullTestSetL1.*B1(:,FitInfo1.Index1SE)',2)+FitInfo1.Intercept(FitInfo1.Index1SE),FullTestSetY).Rsquared.Ordinary;

        R2TrainLa2(j,i) = fitlm(sum(FullTrainningSetL2(TrainingSetIndexSele,:).*B2(:,FitInfo2.Index1SE)',2)+FitInfo2.Intercept(FitInfo2.Index1SE),FullTrainningSetY(TrainingSetIndexSele,:)).Rsquared.Ordinary;
        R2TestLa2(j,i) = fitlm(sum(FullTestSetL2.*B2(:,FitInfo2.Index1SE)',2)+FitInfo2.Intercept(FitInfo2.Index1SE),FullTestSetY).Rsquared.Ordinary;
    end
end


save('Figure_S3_Traning_Splits.mat')

% clear
% load('Figure_S3_Traning_Splits.mat')

% Figure SI

BoxColor=lines;
BoxColor=BoxColor(1:4,:);
BoxColor(3:4,:)=BoxColor(1:2,:)+0.3;
BoxColor(BoxColor>1)=1;

figure('Position',[100,100,1200,800],'Color',[1,1,1])
hold on
for i=2:4

    TrainLa2=R2TrainLa2(:,i);
    TestLa2=R2TestLa2(:,i);

    boxchart(ones(1,100).*(-0.25+i*2),TrainLa2(:),'MarkerStyle','.','MarkerColor',BoxColor(3,:),'BoxFaceColor',BoxColor(3,:),'BoxEdgeColor',[0,0,0],'WhiskerLineColor',[0,0,0],'LineWidth',1.5,'BoxWidth',0.4,'BoxFaceAlpha',0.6)
    boxchart(ones(1,100).*(0.25+i*2),TestLa2(:),'MarkerStyle','.','MarkerColor',BoxColor(4,:),'BoxFaceColor',BoxColor(4,:),'BoxEdgeColor',[0,0,0],'WhiskerLineColor',[0,0,0],'LineWidth',1.5,'BoxWidth',0.4,'BoxFaceAlpha',0.6)

end
axis([3 9 -0.05 1.05])
set(gca,'FontSize',20,'LineWidth',1.5)
set(gca,'layer','top')
xticks(2:2:10)
xticklabels({'2','4','6','8','10'});
xlabel('SynCom size of test sets (number of strains)');
ylabel('R^2')
legend([{'Training'},{'Test'}],'Location','SouthWest')
box on

set(gcf,'PaperType','A3')
print('FigureS3.pdf','-dpdf','-r300')


%% Figure S4 PairBeta
% clear
% load('Figure3.mat')


figure('Position',[100,100,1050,1000],'Color',[1,1,1])
hold on
imagesc(mean(Order2CorrValue,3),[-0.1 0.1])
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
set(gca,'layer','top')
set(gcf,'Renderer','painters')
colorbar
print('FigureS4A.pdf','-dpdf')



Beta=mean(Order2CorrValue,3);

figure('Position',[100,100,800,600],'Color',[1,1,1])
hold on
plot(Beta(Delta~=0),Delta(Delta~=0),'o','MarkerSize',12,'LineWidth',1.5)
plot([-0.1,0.25],[-0.4,1],'--k','LineWidth',1)
set(gca,'FontSize',16,'LineWidth',1.5)
axis([-0.1,0.25,-0.4,1])
box on
xlabel('\beta_{ij}')
ylabel('\delta_{ij}')

set(gcf,'PaperType','A3')
set(gca,'layer','top')
set(gcf,'Renderer','painters')
print('FigureS4B.pdf','-dpdf')




%% Figrue S5 Root401FWvsCFU
% clear

Root401FWvsCFU=readtable('FigureSI_FW_Data.xlsx','Sheet',2);
Root401FWvsCFU=table2array(Root401FWvsCFU(:,2:3));

x=Root401FWvsCFU(:,2);
y=Root401FWvsCFU(:,1);
[Root401C,Root401p]=corr(Root401FWvsCFU(:,2),Root401FWvsCFU(:,1));

[xData, yData] = prepareCurveData( x, y );
ft = fittype( 'poly1' );
[fitresult, gof] = fit( xData, yData, ft );

figure('Position',[100,100,1200,800],'Color',[1,1,1])
hold on
plot(Root401FWvsCFU(:,2),Root401FWvsCFU(:,1),'.','MarkerSize',10)
plot([0,8],fitresult.p1.*[0,8]+fitresult.p2,'-','LineWidth',1.5,'Color',[0.5 0.5 0.5])
box on
axis([3.5 7.5 0 50])
set(gca,'FontSize',24,'LineWidth',1.5)
set(gcf,'PaperType','A3')
set(gcf,'Renderer','painters')
set(gca,'layer','top')
yticks(0:10:50)
xticks(4:8)
xlabel('CFU of Root401')
ylabel('Fresh Weight (mg)')

text(3.6,48,['Pearson = ',num2str(Root401C,3)],'FontSize',20)
text(3.6,45,'p-Value < 0.001','FontSize',20)

print('FigureS5.pdf','-dpdf','-r300')

%% Figure S6 Root401
% clear
% load('Figure4.mat')


figure('Position',[100,100,600,1000],'Color',[1,1,1])
hold on

boxchart(ones(100,1).*1,R2_Test_LA1','LineWidth',1.5,'BoxWidth',0.6,'BoxFaceAlpha',0.6,'MarkerSize',12)
boxchart(ones(100,1).*2,R2_Test_LA2','LineWidth',1.5,'BoxWidth',0.6,'BoxFaceAlpha',0.6,'MarkerSize',12)

axis([0 3 0.5 1.03])
set(gca,'FontSize',28,'LineWidth',1.5)
set(gca,'layer','top')
xticks('')
xlabel('')
ylabel('R^2')
box on
legend({'Order = 1','Order = 2'},'Location','SouthEast','FontSize',22)

set(gcf,'PaperType','A3')
print('FigureS6A.pdf','-dpdf','-r300')

% Number of Parameters
figure('Position',[100,100,600,1000],'Color',[1,1,1])
hold on

boxchart(ones(100,1).*1,LA1ParaNum','LineWidth',1.5,'BoxWidth',0.6,'BoxFaceAlpha',0.6,'MarkerSize',12)
boxchart(ones(100,1).*2,LA2ParaNum','LineWidth',1.5,'BoxWidth',0.6,'BoxFaceAlpha',0.6,'MarkerSize',12)

axis([0 3 0 35])
set(gca,'FontSize',28,'LineWidth',1.5)
set(gca,'layer','top')
xticks('')
xlabel('')
ylabel('Number of Parameters')
box on
legend({'Order = 1','Order = 2'},'Location','SouthEast','FontSize',22)

set(gcf,'PaperType','A3')
print('FigureS6B.pdf','-dpdf','-r300')

%% Figure S7 Training Root401

clear
load('Figure4.mat')


NumberSynComm=[30,30,15,10,5,1];
NumAll=sum(NumberSynComm);

% Training 1
[~,LAYSortIndex]=sort(LassoOrder1y);

TestSetSize=10;

TrainingSetSizeList=[30,40,50,60,70,80];


for i=1:50

    TestSetIndex=LAYSortIndex([randperm(9,1),randperm(9,1)+9,randperm(9,1)+18,randperm(9,1)+27,randperm(9,1)+36,randperm(9,1)+45,randperm(9,1)+54,randperm(9,1)+63,randperm(9,1)+72,randperm(9,1)+81]);
    TrainningSetIndex=1:NumAll;
    TrainningSetIndex(TestSetIndex)=[];

    FullTestSetL1=LassoOrder1X(TestSetIndex,:);
    FullTestSetL2=LassoOrder2X(TestSetIndex,:);
    FullTestSetY=LassoOrder1y(TestSetIndex,:);

    FullTrainningSetL1=LassoOrder1X(TrainningSetIndex,:);
    FullTrainningSetL2=LassoOrder2X(TrainningSetIndex,:);
    FullTrainningSetY=LassoOrder1y(TrainningSetIndex,:);

    for j=1:50
        for k=1:6
            TrainingSetSize=TrainingSetSizeList(k);


            TrainingSetIndexSele=randperm(80,TrainingSetSize);

            [B1,FitInfo1]=lasso(FullTrainningSetL1(TrainingSetIndexSele,:), FullTrainningSetY(TrainingSetIndexSele,:),'cv',10);
            [B2,FitInfo2]=lasso(FullTrainningSetL2(TrainingSetIndexSele,:), FullTrainningSetY(TrainingSetIndexSele,:),'cv',10);

            R2TrainLa1(i,j,k) = fitlm(sum(FullTrainningSetL1(TrainingSetIndexSele,:).*B1(:,FitInfo1.Index1SE)',2)+FitInfo1.Intercept(FitInfo1.Index1SE),FullTrainningSetY(TrainingSetIndexSele,:)).Rsquared.Ordinary;
            R2TestLa1(i,j,k) = fitlm(sum(FullTestSetL1.*B1(:,FitInfo1.Index1SE)',2)+FitInfo1.Intercept(FitInfo1.Index1SE),FullTestSetY).Rsquared.Ordinary;

            R2TrainLa2(i,j,k) = fitlm(sum(FullTrainningSetL2(TrainingSetIndexSele,:).*B2(:,FitInfo2.Index1SE)',2)+FitInfo2.Intercept(FitInfo2.Index1SE),FullTrainningSetY(TrainingSetIndexSele,:)).Rsquared.Ordinary;
            R2TestLa2(i,j,k) = fitlm(sum(FullTestSetL2.*B2(:,FitInfo2.Index1SE)',2)+FitInfo2.Intercept(FitInfo2.Index1SE),FullTestSetY).Rsquared.Ordinary;

        end
    end

end


save('FigureS7_Training_Size.mat')
load('FigureS7_Training_Size.mat')

BoxColor=lines;
BoxColor=BoxColor(1:4,:);
BoxColor(3:4,:)=BoxColor(1:2,:)+0.3;
BoxColor(BoxColor>1)=1;

figure('Position',[100,100,1500,800],'Color',[1,1,1])
hold on
for i=1:6

    TrainLa1=R2TrainLa1(:,:,i);
    TestLa1=R2TestLa1(:,:,i);
    TrainLa2=R2TrainLa2(:,:,i);
    TestLa2=R2TestLa2(:,:,i);

    boxchart(ones(1,50).*(-1+i*10+20),mean(TrainLa2,2),'MarkerStyle','.','MarkerColor',BoxColor(3,:),'BoxFaceColor',BoxColor(3,:),'BoxEdgeColor',[0,0,0],'WhiskerLineColor',[0,0,0],'LineWidth',1.5,'BoxWidth',1.2,'BoxFaceAlpha',0.6)
    boxchart(ones(1,50).*(1+i*10+20),mean(TestLa2,2),'MarkerStyle','.','MarkerColor',BoxColor(4,:),'BoxFaceColor',BoxColor(4,:),'BoxEdgeColor',[0,0,0],'WhiskerLineColor',[0,0,0],'LineWidth',1.5,'BoxWidth',1.2,'BoxFaceAlpha',0.6)

end
axis([25 85 -0.05 1.05])
set(gca,'FontSize',20,'LineWidth',1.5)
set(gca,'layer','top')
xticks(20:10:80)
xlabel('Size of Trainning Set')
ylabel('R^2')
legend([{'Training'},{'Test'}],'Location','SouthEast')

box on
set(gcf,'PaperType','A1')
print('FigureS7A.pdf','-dpdf','-r300')


NumberSynCommSum=0;
NumberSynCommSum(2:7)=cumsum(NumberSynComm);


for i=1:5
    TestSetIndex=NumberSynCommSum(i)+1:NumberSynCommSum(i+1);
    TrainningSetIndex=1:NumAll;
    TrainningSetIndex(TestSetIndex)=[];

    FullTestSetL1=LassoOrder1X(TestSetIndex,:);
    FullTestSetL2=LassoOrder2X(TestSetIndex,:);
    FullTestSetY=LassoOrder1y(TestSetIndex,:);

    FullTrainningSetL1=LassoOrder1X(TrainningSetIndex,:);
    FullTrainningSetL2=LassoOrder2X(TrainningSetIndex,:);
    FullTrainningSetY=LassoOrder1y(TrainningSetIndex,:);

    for j=1:100

        TrainingSetIndexSele=1:size(FullTrainningSetY,1);

        [B1,FitInfo1]=lasso(FullTrainningSetL1(TrainingSetIndexSele,:), FullTrainningSetY(TrainingSetIndexSele,:),'cv',10);
        [B2,FitInfo2]=lasso(FullTrainningSetL2(TrainingSetIndexSele,:), FullTrainningSetY(TrainingSetIndexSele,:),'cv',10);

        [CorrTrainLa1(j,i),PTrainLa1(j,i)]=corr(sum(FullTrainningSetL1(TrainingSetIndexSele,:).*B1(:,FitInfo1.Index1SE)',2)+FitInfo1.Intercept(FitInfo1.Index1SE),FullTrainningSetY(TrainingSetIndexSele,:));
        [CorrTestLa1(j,i),PTestLa1(j,i)]=corr(sum(FullTestSetL1.*B1(:,FitInfo1.Index1SE)',2)+FitInfo1.Intercept(FitInfo1.Index1SE),FullTestSetY);

        R2TrainLa1(j,i) = fitlm(sum(FullTrainningSetL1(TrainingSetIndexSele,:).*B1(:,FitInfo1.Index1SE)',2)+FitInfo1.Intercept(FitInfo1.Index1SE),FullTrainningSetY(TrainingSetIndexSele,:)).Rsquared.Ordinary;
        R2TestLa1(j,i) = fitlm(sum(FullTestSetL1.*B1(:,FitInfo1.Index1SE)',2)+FitInfo1.Intercept(FitInfo1.Index1SE),FullTestSetY).Rsquared.Ordinary;

        [CorrTrainLa2(j,i),PTrainLa2(j,i)]=corr(sum(FullTrainningSetL2(TrainingSetIndexSele,:).*B2(:,FitInfo2.Index1SE)',2)+FitInfo2.Intercept(FitInfo2.Index1SE),FullTrainningSetY(TrainingSetIndexSele,:));
        [CorrTestLa2(j,i),PTestLa2(j,i)]=corr(sum(FullTestSetL2.*B2(:,FitInfo2.Index1SE)',2)+FitInfo2.Intercept(FitInfo2.Index1SE),FullTestSetY);

        R2TrainLa2(j,i) = fitlm(sum(FullTrainningSetL2(TrainingSetIndexSele,:).*B2(:,FitInfo2.Index1SE)',2)+FitInfo2.Intercept(FitInfo2.Index1SE),FullTrainningSetY(TrainingSetIndexSele,:)).Rsquared.Ordinary;
        R2TestLa2(j,i) = fitlm(sum(FullTestSetL2.*B2(:,FitInfo2.Index1SE)',2)+FitInfo2.Intercept(FitInfo2.Index1SE),FullTestSetY).Rsquared.Ordinary;

    end

end


BoxColor=lines;
BoxColor=BoxColor(1:4,:);
BoxColor(3:4,:)=BoxColor(1:2,:)+0.3;
BoxColor(BoxColor>1)=1;
BoxColor=BoxColor([1,3,2,4],:);
save('FigureS7_Training_Splits.mat')


load('FigureS7_Training_Splits.mat')

figure('Position',[100,100,1200,600],'Color',[1,1,1])
hold on
for i=1:5

    TrainLa1=R2TrainLa1(:,i);
    TestLa1=R2TestLa1(:,i);
    TrainLa2=R2TrainLa2(:,i);
    TestLa2=R2TestLa2(:,i);

    boxchart(ones(1,100).*(-0.15+i*2),TrainLa2(:),'MarkerStyle','.','MarkerColor',BoxColor(3,:),'BoxFaceColor',BoxColor(3,:),'BoxEdgeColor',[0,0,0],'WhiskerLineColor',[0,0,0],'LineWidth',1.5,'BoxWidth',0.15,'BoxFaceAlpha',0.6)
    boxchart(ones(1,100).*(0.15+i*2),TestLa2(:),'MarkerStyle','.','MarkerColor',BoxColor(4,:),'BoxFaceColor',BoxColor(4,:),'BoxEdgeColor',[0,0,0],'WhiskerLineColor',[0,0,0],'LineWidth',1.5,'BoxWidth',0.15,'BoxFaceAlpha',0.6)

end
axis([1 11 -0.05 1.05])
set(gca,'FontSize',20,'LineWidth',1.5)
set(gca,'layer','top')
xticks(2:2:10)
xlabel('Size of Trainning Set')
ylabel('R^2')
legend([{'Training'},{'Test'}],'Location','SouthEast')

box on

set(gcf,'PaperType','A3')
print('FigureS7.pdf','-dpdf','-r300')



%% Figure S8 Root401 Individual
clear
load('Figure4.mat')

figure('Position',[100,100,1000,1000],'Color',[1,1,1])
hold on

hold on
plot(-1,-1,'o','MarkerEdgeColor',BoxColorGroup(2,:),'MarkerFaceColor',BoxColorGroup(2,:),'LineWidth',1.5)
plot(-1,-1,'o','MarkerEdgeColor',BoxColorGroup(7,:),'MarkerFaceColor',BoxColorGroup(7,:),'LineWidth',1.5)
plot(-1,-1,'o','MarkerEdgeColor',BoxColorGroup(5,:),'MarkerFaceColor',BoxColorGroup(5,:),'LineWidth',1.5)
plot(-1,-1,'o','MarkerEdgeColor',BoxColorGroup(11,:),'MarkerFaceColor',BoxColorGroup(11,:),'LineWidth',1.5)

for i=1:7
    fill([i*2-0.5,i*2-0.5,i*2+0.5,i*2+0.5],[-5 2 2 -5],[0.85 0.85 0.85],'LineStyle','none')
end
plot([0.5,15.5],[0,0],'--','Color',[0.3 0.3 0.3],'LineWidth',1)

for i=1:15
    boxchart(ones(size(cell2mat(CFU(MeanCFUID(i),:)'))).*i,cell2mat(CFU(MeanCFUID(i),:)'),'BoxFaceColor',BoxColor(MeanCFUID(i),:),'WhiskerLineColor',BoxColor(MeanCFUID(i),:),'LineWidth',1.5,'BoxWidth',0.5,'BoxFaceAlpha',0.5,'MarkerStyle','none');

    scatter(ones(size(CFU{MeanCFUID(i),1})).*i+(rand(size(CFU{MeanCFUID(i),1}))-0.5).*0.5,CFU{MeanCFUID(i),1},30,'o','MarkerEdgeColor',BoxColor(MeanCFUID(i),:),'MarkerFaceColor',BoxColor(MeanCFUID(i),:),'LineWidth',1.5,'MarkerFaceAlpha',0.5);
    scatter(ones(size(CFU{MeanCFUID(i),2})).*i+(rand(size(CFU{MeanCFUID(i),2}))-0.5).*0.5,CFU{MeanCFUID(i),2},30,'s','MarkerEdgeColor',BoxColor(MeanCFUID(i),:),'MarkerFaceColor',BoxColor(MeanCFUID(i),:),'LineWidth',1.5,'MarkerFaceAlpha',0.5);
end

axis([0.5 15.5 -2.5 0.5])
set(gca,'FontSize',24,'LineWidth',1.5)
set(gca,'layer','top')
xticks(1:15)
xticklabels(StrainGenusFull(MeanCFUID))
box on
ylabel('Root 401 CFU Decrease (log 10)')
legend(GroupName,'Location','SouthEast','FontSize',18)
set(gcf,'PaperType','A3')

print('FigureS8A.pdf','-dpdf','-r300')


figure('Position',[100,100,1200,1000],'Color',[1,1,1])
hold on
for i=1:14
    errorbar(mean(MeanCFUPreLa1(i,:),2),mean(MeanCFU(i+1,:),2,'omitnan'),std(MeanCFU(i+1,:),[],2,'omitnan'),std(MeanCFU(i+1,:),[],2,'omitnan'),std(MeanCFUPreLa1(i,:),[],2),std(MeanCFUPreLa1(i,:),[],2),'o','Color',BoxColor(i+1,:),'LineWidth',1.5,'MarkerSize',12)
end

plot([-7.5 0.5],[-7.5 0.5],'--','Color',[0.5 0.5 0.5])

box on
axis equal
axis([-2.4 0.4 -2.4 0.4])
set(gca,'FontSize',28,'LineWidth',1.5)
set(gca,'layer','top')
xlabel('Predicted by Regression (log 10)')
ylabel('Measured CFU Decrease (log 10)')
legend(StrainListFull(2:15),'Location','SouthEast','FontSize',18)

text(-2,0.25,['R^2 = ',num2str(R2MeanCFULa1,3)],'FontSize',24)
set(gcf,'PaperType','A3')
print('Figure8B.pdf','-dpdf','-r300')


figure('Position',[100,100,1200,1000],'Color',[1,1,1])

hold on

for i=1:14
    errorbar(mean(MeanCFUPreLa2(i,:),2),mean(MeanCFU(i+1,:),2,'omitnan'),std(MeanCFU(i+1,:),[],2,'omitnan'),std(MeanCFU(i+1,:),[],2,'omitnan'),std(MeanCFUPreLa2(i,:),[],2),std(MeanCFUPreLa2(i,:),[],2),'o','Color',BoxColor(i+1,:),'LineWidth',1.5,'MarkerSize',12)
end

plot([-7.5 0.5],[-7.5 0.5],'--','Color',[0.5 0.5 0.5])

box on
axis equal
axis([-2.4 0.4 -2.4 0.4])
set(gca,'FontSize',28,'LineWidth',1.5)
set(gca,'layer','top')
xlabel('Predicted by Regression (log 10)')
ylabel('Measured CFU Decrease (log 10)')
legend(StrainListFull(2:15),'Location','SouthEast','FontSize',18)

text(-2,0.25,['R^2 = ',num2str(R2MeanCFULa2,3)],'FontSize',24)

set(gcf,'PaperType','A3')
print('Figure8C.pdf','-dpdf','-r300')

%% Figure S9 Batch Effect

clear
load('Figure3.mat')

Test1Mean=[];
Test2Mean=[];
Test1STD=[];
Test2STD=[];


Test1Mean=mean(MeanCFU(2:end,1:3),2,'omitnan');
Test1STD=std(MeanCFU(2:end,1:3),[],2,'omitnan');
Test2Mean=MeanCFU(2:end,4);
Test2STD=cellfun(@std,CFU(2:end,4));



for i=1:67
    Combination=TableForLasso(i,4:17);
    Index=sum(table2array(Table(:,4:17))==Combination,2)==14;
    if size(unique(table2array(Table(Index,1))),1)==2

        Test1Mean(end+1)=mean(table2array(Table(table2array(Table(Index,1))==1,18)));
        Test1STD(end+1)=std(table2array(Table(table2array(Table(Index,1))==1,18)));
        Test2Mean(end+1)=mean(table2array(Table(table2array(Table(Index,1))==2,18)));
        Test2STD(end+1)=std(table2array(Table(table2array(Table(Index,1))==2,18)));

    else
    end
end

[corr_batch,p]=corr(Test1Mean,Test2Mean);


figure('Position',[100,100,1000,800],'Color',[1,1,1])
hold on

errorbar(Test1Mean,Test2Mean,Test2STD,Test2STD,Test1STD,Test1STD,'Color',BoxColor2(1,:).*0.8,'LineStyle','none','MarkerSize',10,'LineWidth',1.5)


plot([-5,1],[-5,1],'--','Color',[0.3 0.3 0.3],'LineWidth',1)
axis equal
axis([-5 1 -5 1])
set(gca,'FontSize',24,'LineWidth',1.5)
set(gca,'layer','top')
xlabel('Batch 1')
ylabel('Batch 2')
text(-4.7,0.65,['Pearson = ',num2str(corr_batch)],'FontSize',24)
text(-4.7,0.3,'p-Value < 0.001','FontSize',24)

box on
% print(['C:\Users\dell\Desktop\Batch Effect.jpg'],'-djpeg','-r300')
set(gcf,'PaperType','A3')
print('FigureS9.pdf','-dpdf','-r300')





%% Figure S10 FW

FreshWeightTable=readtable('FigureSI_FW_Data.xlsx','Sheet',3);

FreshWeight=table2array(FreshWeightTable(:,2:5));


IndexAll=cell(0);
IndexAll{1}=93;
IndexAll{2}=95:108;
IndexAll{3}=1:30;
IndexAll{4}=31:60;
IndexAll{5}=62:76;
IndexAll{6}=77:86;
IndexAll{7}=87:91;
IndexAll{8}=94;


BoxColor2=[0.9:(-0.75/8):0.15;1:(-0.45/8):0.55;0.8:(-0.75/8):0.05]';


figure('Position',[100,100,1600,600],'Color',[1,1,1])
hold on

X_0=1;
for k=1:8

    Number=size(IndexAll{k},2);

    FreshWeightPlot=FreshWeight(IndexAll{k},1);
    FreshWeightPlotSTD=FreshWeight(IndexAll{k},2);

    [FreshWeightPlot,FreshWeightPlotIndex]=sort(FreshWeightPlot);


    bar(X_0:(X_0+Number-1),FreshWeightPlot,'FaceColor',BoxColor2(k,:))
    errorbar(X_0:(X_0+Number-1),FreshWeightPlot,FreshWeightPlotSTD(FreshWeightPlotIndex),'Color',[0,0,0],'LineStyle','none')

    X_0=X_0+Number;

    plot([X_0-0.5,X_0-0.5],[0,45],':k')

end

axis([0,107,0,45])
box on
xticks([])
set(gca,'FontSize',24,'LineWidth',1.5)
set(gcf,'PaperType','A3')
set(gcf,'Renderer','painters')
set(gca,'layer','top')
ylabel('Fresh Weight (mg)')

print('FigureS10.pdf','-dpdf','-r300')



%% Figure S16

clear

% Set seed for reproducibility
rng(123);

% Parameters of the dynamics
NumberSpecies = 15; % number of samples
steps = 5000;
dt = 0.01;
MeanA = 0.3; % mean of pairwise matrix A
GrowthRate = normrnd(1, 0); % growth rate


A = -ones(2,2);

Aij=-1:0.02:1;
Aji=-1:0.02:1;


Collection = [];
Abundance = [];

FinalAbundance=[];


for i=1:101
    for j=1:101
        A(1,2)=Aij(i);
        A(2,1)=Aji(j);

        Abundance = zeros(steps+1, 2);
        Abundance(1, :) = 0.01;

        for Time = 1:steps
            Abundance(Abundance < 0) = 0; % Ensure abundance cannot be negative
            pair_term = sum(A .* Abundance(Time, :),2)'; % Pairwise interaction term
            % Ncomb = Abundance(t, :)' * Abundance(t, :); % Cross product
            Abundance(Time + 1, :) = Abundance(Time, :) + dt * Abundance(Time, :) .* (GrowthRate + pair_term);
        end

        FinalAbundance(i,j) = Abundance(end, 1); % Final abundance
    end
    i
end


cmap=[];
cmap(1:128,1)=0.85:0.15/127:1;
cmap(1:128,2)=0.23:0.77/127:1;
cmap(1:128,3)=0.10:0.90/127:1;
cmap(129:256,1)=1:-1.00/127:0.00;
cmap(129:256,2)=1:-0.65/127:0.35;
cmap(129:256,3)=1:-0.15/127:0.85;

save('FigureS16A.mat')
load('FigureS16A.mat')


figure('Position',[100,100,1050,1000],'Color',[1,1,1])
imagesc([-1,0.98],[-1,0.98],log10(FinalAbundance(2:end-1,2:end-1)),[-1,1])
colormap(cmap)
colorbar
axis equal
axis([-1.01 0.99 -1.01 0.99])
colorbar
xlabel('a_{i0}')
ylabel('a_{0i}')
set(gca,'FontSize',24,'LineWidth',1.5)
set(gcf,'PaperType','A3')
set(gca,'layer','top')
set(gcf,'Renderer','painters')
print('FigureS16A.pdf','-dpdf')


load('FigureS16A.mat')

Additive=log10(FinalAbundance(26,26))+log10(FinalAbundance(26,26));

% Set seed for reproducibility
rng(123);

% Parameters of the dynamics
NumberSpecies = 15; % number of samples
steps = 5000;
dt = 0.01;
MeanA = 0.3; % mean of pairwise matrix A
GrowthRate = normrnd(1, 0); % growth rate


A = -ones(3,3);

Ajk=-1:0.02:1;
Akj=-1:0.02:1;

A(1,2:3)=[-0.5,-0.5];
A(2:3,1)=-0.5;

Collection = [];
Abundance = [];

FinalAbundance=[];


for i=1:101
    for j=1:101

        A(2,3)=Ajk(i);
        A(3,2)=Akj(j);

        Abundance = zeros(steps+1, 3);
        Abundance(1, :) = 0.01;

        for Time = 1:steps
            Abundance(Abundance < 0) = 0; % Ensure abundance cannot be negative
            pair_term = sum(A .* Abundance(Time, :),2)'; % Pairwise interaction term
            % Ncomb = Abundance(t, :)' * Abundance(t, :); % Cross product
            Abundance(Time + 1, :) = Abundance(Time, :) + dt * Abundance(Time, :) .* (GrowthRate + pair_term);
        end

        FinalAbundance(i,j) = Abundance(end, 1); % Final abundance
    end
    
end


cmap=[];
cmap(1:128,1)=0.85:0.15/127:1;
cmap(1:128,2)=0.23:0.77/127:1;
cmap(1:128,3)=0.10:0.90/127:1;
cmap(129:256,1)=1:-1.00/127:0.00;
cmap(129:256,2)=1:-0.65/127:0.35;
cmap(129:256,3)=1:-0.15/127:0.85;

save('FigureS16B2.mat')

load('FigureS16B1.mat')

figure('Position',[100,100,1050,1000],'Color',[1,1,1])
imagesc([-1,1],[-1,1],log10(FinalAbundance)-Additive,[-0.4,0.4])
colormap(cmap)
colorbar
axis equal
axis([-0.81 0.11 -0.81 0.11])
colorbar
xlabel('a_{ij}')
ylabel('a_{ji}')
set(gca,'FontSize',24,'LineWidth',1.5)
set(gcf,'PaperType','A3')
set(gca,'layer','top')
set(gcf,'Renderer','painters')
print('FigureS16B1.pdf','-dpdf')


load('FigureS16B2.mat')

figure('Position',[100,100,1050,1000],'Color',[1,1,1])
imagesc([-1,1],[-1,1],log10(FinalAbundance)-Additive,[-0.4,0.4])
colormap(cmap)
colorbar
axis equal
axis([-0.81 0.11 -0.81 0.11])
colorbar
xlabel('a_{ij}')
ylabel('a_{ji}')
set(gca,'FontSize',24,'LineWidth',1.5)
set(gcf,'PaperType','A3')
set(gca,'layer','top')
set(gcf,'Renderer','painters')
print('FigureS16B2.pdf','-dpdf')


%% Figure S14

clear
clc

% Parameters of the dynamics
NumberSpecies = 15; % number of samples
steps = 5000;
dt = 0.01;

% Set seed for reproducibility
rng(123);
GrowthRate = normrnd(1,0.2,1,15); % growth rate
GrowthRate(1) = 1;


MeanAAll = [-0.1,-0.3,-0.5];


for MeanAIndex = 2
    MeanA = MeanAAll(MeanAIndex); % mean of pairwise matrix A
    A_STD = 0.1; % standard deviation of pairwise matrix A

    rng(123);
    A = normrnd(MeanA, A_STD, NumberSpecies, NumberSpecies); % Initial interaction matrix

    for i = 1:NumberSpecies
        A(i, i) = -1;
    end


    Collection = [];
    Abundance = [];

    FinalAbundance=[];
    for SynComSize = 1:14
        Combination = nchoosek(2:NumberSpecies, SynComSize);

        if SynComSize == 2
            % Combination = Combination;
        elseif size(Combination, 1) > 50
            Combination = Combination(randsample(1:size(Combination, 1), 50),:);
        end
        for j = 1:size(Combination, 1)
            Abundance = zeros(steps+1, NumberSpecies);
            Abundance(1, :) = 0;
            % rng(SynComSize);
            Abundance(1, [1 Combination(j,:)]) = 0.01;

            for Time = 1:steps
                Abundance(Abundance < 0) = 0; % Ensure abundance cannot be negative
                pair_term = sum(A .* Abundance(Time, :),2)'; % Pairwise interaction term
                Abundance(Time + 1, :) = Abundance(Time, :) + dt * Abundance(Time, :) .* (GrowthRate + pair_term);
            end

            FinalAbundance = [FinalAbundance; Abundance(end, 1)]; % Final abundance
            Collection_sub = zeros(1, NumberSpecies);
            Collection_sub(Combination(j,:)) = 1;
            Collection = [Collection; Collection_sub]; % Species combination
        end
    end


    MatrixATable = A;
    AbundanceGroup{MeanAIndex} = cell(0);


    for i = 1:15
        AbundanceGroup{MeanAIndex}{i}=FinalAbundance(sum(Collection,2)==i);
    end
end


LassoOrder1=Collection(15:end,2:end);
LassoOrder1(LassoOrder1==0)=-1;
LassoOrder2=[];
LassoOrder3=[];

for i = 1:size(LassoOrder1,1)
    Order2Matrix=nchoosek(LassoOrder1(i,:),2);
    LassoOrder2(i,:)=[LassoOrder1(i,:),(Order2Matrix(:,1).*Order2Matrix(:,2))'];
    Order3Matrix=nchoosek(LassoOrder1(i,:),3);
    LassoOrder3(i,:)=[LassoOrder2(i,:),(Order3Matrix(:,1).*Order3Matrix(:,2).*Order3Matrix(:,3))'];
end

FinalAbundance=cell2mat(AbundanceGroup{2}');
LassoOrderY=log10(FinalAbundance(15:end));

R2La1All=[];
R2La2All=[];
R2La1TrainAll=[];
R2La1TestAll=[];
R2La2TrainAll=[];
R2La2TestAll=[];

R2SingleLA2All=[];
DeltaR2All=[];

TrainingSetSize=[30,50,75,100,150,200,250,300];

for TrainingSetIndex=1:8
    for test=1:100
        Index=randperm(size(LassoOrderY,1),TrainingSetSize(TrainingSetIndex));

        LassoOrder1Test=LassoOrder1;
        LassoOrder1Test(Index,:)=[];
        LassoOrder2Test=LassoOrder2;
        LassoOrder2Test(Index,:)=[];
        LassoOrderYTest=LassoOrderY;
        LassoOrderYTest(Index,:)=[];

        LassoOrder1Train=LassoOrder1(Index,:);
        LassoOrder2Train=LassoOrder2(Index,:);
        LassoOrderYTrain=LassoOrderY(Index,:);


        [B1,FitInfo1]=lasso(LassoOrder1Train, LassoOrderYTrain,'cv',10);
        [B2,FitInfo2]=lasso(LassoOrder2Train, LassoOrderYTrain,'cv',10);


        LassoOrder1Single=Collection(1:14,2:end);
        LassoOrder1Single(LassoOrder1Single==0)=-1;

        LassoOrder2Single=[];

        for i = 1:size(LassoOrder1Single,1)
            Order2Matrix=nchoosek(LassoOrder1Single(i,:),2);
            LassoOrder2Single(i,:)=[LassoOrder1Single(i,:),(Order2Matrix(:,1).*Order2Matrix(:,2))'];
        end


        LA2PredPair=sum(LassoOrder2.*B2(:,FitInfo2.Index1SE)',2)+FitInfo2.Intercept(FitInfo2.Index1SE);
        LA2PredPair=LA2PredPair(1:91);
        gLVPari=LassoOrderY(1:91);
        PairCollection=Collection(sum(Collection,2)==2,:);

        gLVSingle=log10(FinalAbundance(1:14));
        LA2Single=sum(LassoOrder2Single.*B2(:,FitInfo2.Index1SE)',2)+FitInfo2.Intercept(FitInfo2.Index1SE);

        for i=1:91
            Delta_gLV(i)=gLVPari(i)-sum(gLVSingle(PairCollection(i,2:14)==1));
            Delta_LA2(i)=LA2PredPair(i)-sum(LA2Single(PairCollection(i,2:14)==1));
        end

        R2La1=fitlm(sum(LassoOrder1.*B1(:,FitInfo1.Index1SE)',2)+FitInfo1.Intercept(FitInfo1.Index1SE),LassoOrderY).Rsquared.Ordinary;
        R2La2=fitlm(sum(LassoOrder2.*B2(:,FitInfo2.Index1SE)',2)+FitInfo2.Intercept(FitInfo2.Index1SE),LassoOrderY).Rsquared.Ordinary;

        R2La1Train=fitlm(sum(LassoOrder1Train.*B1(:,FitInfo1.Index1SE)',2)+FitInfo1.Intercept(FitInfo1.Index1SE),LassoOrderYTrain).Rsquared.Ordinary;
        R2La1Test=fitlm(sum(LassoOrder1Test.*B1(:,FitInfo1.Index1SE)',2)+FitInfo1.Intercept(FitInfo1.Index1SE),LassoOrderYTest).Rsquared.Ordinary;
        R2La2Train=fitlm(sum(LassoOrder2Train.*B2(:,FitInfo2.Index1SE)',2)+FitInfo2.Intercept(FitInfo2.Index1SE),LassoOrderYTrain).Rsquared.Ordinary;
        R2La2Test=fitlm(sum(LassoOrder2Test.*B2(:,FitInfo2.Index1SE)',2)+FitInfo2.Intercept(FitInfo2.Index1SE),LassoOrderYTest).Rsquared.Ordinary;

        R2SingleLA2=fitlm(sum(LassoOrder2Single.*B2(:,FitInfo2.Index1SE)',2)+FitInfo2.Intercept(FitInfo2.Index1SE),log10(FinalAbundance(1:14))).Rsquared.Ordinary;
        DeltaR2=fitlm(Delta_gLV,Delta_LA2).Rsquared.Ordinary;

        R2La1All(TrainingSetIndex,test)=R2La1;
        R2La2All(TrainingSetIndex,test)=R2La2;
        R2La1TrainAll(TrainingSetIndex,test)=R2La1Train;
        R2La1TestAll(TrainingSetIndex,test)=R2La1Test;
        R2La2TrainAll(TrainingSetIndex,test)=R2La2Train;
        R2La2TestAll(TrainingSetIndex,test)=R2La2Test;

        R2SingleLA2All(TrainingSetIndex,test)=R2SingleLA2;
        DeltaR2All(TrainingSetIndex,test)=DeltaR2;

    end
end

Colors=lines;

BoxPlotColor(1:3,:)=Colors(1:3,:)+0.2;
BoxPlotColor(4:6,:)=Colors(1:3,:)-0.2;
BoxPlotColor(BoxPlotColor<0)=0;
BoxPlotColor(BoxPlotColor>1)=1;

save('FigureS14.mat')

figure('Position',[100,100,1600,800],'Color',[1,1,1])
hold on
for i =1:8
    boxchart(ones(1,100).*(-0.3+i),R2La2TrainAll(i,:),'MarkerStyle','.','MarkerColor',BoxPlotColor(4,:),'BoxFaceColor',BoxPlotColor(4,:),'BoxEdgeColor',[0,0,0],'WhiskerLineColor',[0,0,0],'LineWidth',1.5,'BoxWidth',0.2,'BoxFaceAlpha',0.6)
    boxchart(ones(1,100).*(0+i),R2La2TestAll(i,:),'MarkerStyle','.','MarkerColor',BoxPlotColor(5,:),'BoxFaceColor',BoxPlotColor(5,:),'BoxEdgeColor',[0,0,0],'WhiskerLineColor',[0,0,0],'LineWidth',1.5,'BoxWidth',0.2,'BoxFaceAlpha',0.6)
    boxchart(ones(1,100).*(0.3+i),R2La2All(i,:),'MarkerStyle','.','MarkerColor',BoxPlotColor(6,:),'BoxFaceColor',BoxPlotColor(6,:),'BoxEdgeColor',[0,0,0],'WhiskerLineColor',[0,0,0],'LineWidth',1.5,'BoxWidth',0.2,'BoxFaceAlpha',0.6)
end
axis([0.5 8.5 0.68 1.02])
set(gca,'FontSize',20,'LineWidth',1.5)

set(gca,'layer','top')
set(gcf,'PaperType','A3')
set(gcf,'Renderer','painters')
xlabel('Trainning Set Size')
xticklabels(TrainingSetSize)
ylabel('R^2')
box on
legend({'Trainning','Test','All'},'Location','SouthEast','FontSize',22)
print('FigureS14A.pdf','-dpdf','-r300')



figure('Position',[100,100,1600,800],'Color',[1,1,1])
hold on
for i =1:8
    boxchart(ones(1,100).*(-0+i),R2SingleLA2All(i,:),'MarkerStyle','.','MarkerColor',BoxPlotColor(1,:),'BoxFaceColor',BoxPlotColor(1,:),'BoxEdgeColor',[0,0,0],'WhiskerLineColor',[0,0,0],'LineWidth',1.5,'BoxWidth',0.2,'BoxFaceAlpha',0.6)
end
axis([0.5 8.5 0 1])
set(gca,'FontSize',24,'LineWidth',1.5)
set(gca,'layer','top')
set(gcf,'PaperType','A3')
xticklabels(TrainingSetSize)
set(gcf,'Renderer','painters')
xlabel('Trainning Set Size')
ylabel('R^2')
box on
legend({'Individual'},'Location','SouthEast','FontSize',22)

print('FigureS14B.pdf','-dpdf','-r300')
