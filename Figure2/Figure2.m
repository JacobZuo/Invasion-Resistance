clear
clc
load('Figure1.mat')

NumAll=sum(NumberSynComm);


LassoOrder1X=TableForLasso(:,4:17).*2-1;
LassoOrder2X=[];

for j=1:NumAll
    Order2Matrix=nchoosek(LassoOrder1X(j,:),2);
    LassoOrder2X(j,:)=[LassoOrder1X(j,:),(Order2Matrix(:,1).*Order2Matrix(:,2))'];
end

LassoOrder1y=TableForLasso(:,18);

% LASSO LOO
[LOOB1,LOOFitInfo1]=lasso(LassoOrder1X, LassoOrder1y,'cv',size(LassoOrder1y,1));
[LOOB2,LOOFitInfo2]=lasso(LassoOrder2X, LassoOrder1y,'cv',size(LassoOrder1y,1));

LOOLA1Result=sum(LassoOrder1X.*LOOB1(:,LOOFitInfo1.Index1SE)',2)+LOOFitInfo1.Intercept(LOOFitInfo1.Index1SE);
LOOLA2Result=sum(LassoOrder2X.*LOOB2(:,LOOFitInfo2.Index1SE)',2)+LOOFitInfo2.Intercept(LOOFitInfo2.Index1SE);

LOOLA1ParaNum=sum(LOOB1(:,LOOFitInfo1.Index1SE)~=0);
LOOLA2ParaNum=sum(LOOB2(:,LOOFitInfo2.Index1SE)~=0);

LOOR2_LA1=fitlm(LOOLA1Result,LassoOrder1y).Rsquared.Ordinary;
LOOR2_LA2=fitlm(LOOLA2Result,LassoOrder1y).Rsquared.Ordinary;

BoxColor2=[0.6:(-0.45/13):0.15;0.9:(-0.35/13):0.55;0.6:(-0.55/13):0.05]';


% Figure 2A
figure('Position',[100,100,1200,1000],'Color',[1,1,1])
hold on
for i=[1:5,7]
    errorbar(LOOLA1Result(TableForLasso(:,2)==2*i),TableForLasso(TableForLasso(:,2)==2*i,18),TableForLasso(TableForLasso(:,2)==2*i,19),'o','Color',BoxColor2(i*2,:),'MarkerSize',10,'LineWidth',1.5)
end

plot([-7.5 0.5],[-7.5 0.5],'--','Color',[0.5 0.5 0.5])

box on
axis equal
axis([-5.5 0.25 -5.5 0.25])
set(gca,'FontSize',28,'LineWidth',1.5)
set(gca,'layer','top')
xlabel('Predicted by Regression (log 10)')
ylabel('Measured CFU Decrease (log 10)')
legend({'2 Strain','4 Strain','6 Strain','8 Strain','10 Strain','14 Strain'},'Location','SouthEast','FontSize',22)
text(-5,-0.2,['R^2 = ',num2str(LOOR2_LA1,3)],'FontSize',24)
set(gcf,'PaperType','A3')
print('Figure2A.pdf','-dpdf','-r300')


% Figure 2B
figure('Position',[100,100,1200,1000],'Color',[1,1,1])
hold on
for i=[1:5,7]
    errorbar(LOOLA2Result(TableForLasso(:,2)==2*i),TableForLasso(TableForLasso(:,2)==2*i,18),TableForLasso(TableForLasso(:,2)==2*i,19),'o','Color',BoxColor2(i*2,:),'MarkerSize',10,'LineWidth',1.5)
end

plot([-7.5 0.5],[-7.5 0.5],'--','Color',[0.5 0.5 0.5])

box on
axis equal
axis([-5.5 0.25 -5.5 0.25])
set(gca,'FontSize',28,'LineWidth',1.5)
set(gca,'layer','top')
xlabel('Predicted by Regression (log 10)')
ylabel('Measured CFU Decrease (log 10)')
legend({'2 Strain','4 Strain','6 Strain','8 Strain','10 Strain','14 Strain'},'Location','SouthEast','FontSize',22)
text(-5,-0.2,['R^2 = ',num2str(LOOR2_LA2,3)],'FontSize',24)
set(gcf,'PaperType','A3')
print('Figure2B.pdf','-dpdf','-r300')

save('Figure2.mat')


%% Trainning Split 1
% clear
% load('Figure2.mat')

[~,LAYSortIndex]=sort(LassoOrder1y);

TestSetSize=7;

TrainingSetSizeList=[20,30,40,50,60];

for i=1:50

    TestSetIndex=LAYSortIndex([randperm(10,1),randperm(9,1)+10,randperm(10,1)+19,randperm(9,1)+29,randperm(10,1)+38,randperm(9,1)+48,randperm(9,1)+57]);
    % TestSetIndex=randperm(67,TestSetSize);
    TrainningSetIndex=1:NumAll;
    TrainningSetIndex(TestSetIndex)=[];

    FullTestSetL1=LassoOrder1X(TestSetIndex,:);
    FullTestSetL2=LassoOrder2X(TestSetIndex,:);
    FullTestSetY=LassoOrder1y(TestSetIndex,:);

    FullTrainningSetL1=LassoOrder1X(TrainningSetIndex,:);
    FullTrainningSetL2=LassoOrder2X(TrainningSetIndex,:);
    FullTrainningSetY=LassoOrder1y(TrainningSetIndex,:);

    for j=1:50
        for k=1:5
            TrainingSetSize=TrainingSetSizeList(k);

            TrainingSetIndexSele=randperm(60,TrainingSetSize);
            CoverRate=sum(sum(FullTrainningSetL1(TrainingSetIndexSele,:)==1)>0);
            Try=1;
            while (CoverRate<14) && (Try<=5)
                TrainingSetIndexSele=randperm(60,TrainingSetSize);
                CoverRate=sum(sum(FullTrainningSetL1(TrainingSetIndexSele,:)==1)>0);
                Try=Try+1;
            end

            [B1,FitInfo1]=lasso(FullTrainningSetL1(TrainingSetIndexSele,:), FullTrainningSetY(TrainingSetIndexSele,:),'cv',10);
            [B2,FitInfo2]=lasso(FullTrainningSetL2(TrainingSetIndexSele,:), FullTrainningSetY(TrainingSetIndexSele,:),'cv',10);

            R2TrainLa1(i,j,k) = fitlm(sum(FullTrainningSetL1(TrainingSetIndexSele,:).*B1(:,FitInfo1.Index1SE)',2)+FitInfo1.Intercept(FitInfo1.Index1SE),FullTrainningSetY(TrainingSetIndexSele,:)).Rsquared.Ordinary;
            R2TestLa1(i,j,k) = fitlm(sum(FullTestSetL1.*B1(:,FitInfo1.Index1SE)',2)+FitInfo1.Intercept(FitInfo1.Index1SE),FullTestSetY).Rsquared.Ordinary;

            R2TrainLa2(i,j,k) = fitlm(sum(FullTrainningSetL2(TrainingSetIndexSele,:).*B2(:,FitInfo2.Index1SE)',2)+FitInfo2.Intercept(FitInfo2.Index1SE),FullTrainningSetY(TrainingSetIndexSele,:)).Rsquared.Ordinary;
            R2TestLa2(i,j,k) = fitlm(sum(FullTestSetL2.*B2(:,FitInfo2.Index1SE)',2)+FitInfo2.Intercept(FitInfo2.Index1SE),FullTestSetY).Rsquared.Ordinary;

        end
    end

end

save('Figure_2C_Traning_Size.mat')
% clear
% load('Figure_2C_Traning_Size.mat')

% Figure 2C

BoxColor=lines;
BoxColor=BoxColor(1:4,:);
BoxColor(3:4,:)=BoxColor(1:2,:)+0.3;
BoxColor(BoxColor>1)=1;

figure('Position',[100,100,1200,1000],'Color',[1,1,1])
hold on
for i=1:5

    TrainLa2=R2TrainLa2(:,:,i);
    TestLa2=R2TestLa2(:,:,i);

    boxchart(ones(1,50).*(-1.2+i*10+10),mean(TrainLa2,2),'MarkerStyle','.','MarkerColor',BoxColor(3,:),'BoxFaceColor',BoxColor(3,:),'BoxEdgeColor',[0,0,0],'WhiskerLineColor',[0,0,0],'LineWidth',1.5,'BoxWidth',1.6,'BoxFaceAlpha',0.6)
    boxchart(ones(1,50).*(1.2+i*10+10),mean(TestLa2,2),'MarkerStyle','.','MarkerColor',BoxColor(4,:),'BoxFaceColor',BoxColor(4,:),'BoxEdgeColor',[0,0,0],'WhiskerLineColor',[0,0,0],'LineWidth',1.5,'BoxWidth',1.6,'BoxFaceAlpha',0.6)

end
axis([15 65 -0.05 1.05])
set(gca,'FontSize',20,'LineWidth',1.5)
set(gca,'layer','top')
xticks(20:10:60)
xticklabels([{'20'},{'30'},{'40'},{'50'},{'60'}])

xlabel('Size of Trainning Set')
ylabel('R^2')
legend([{'Training'},{'Test'}],'Location','SouthEast')

box on
set(gcf,'PaperType','A1')
print('Figure2C.pdf','-dpdf','-r300')



