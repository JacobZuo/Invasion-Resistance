close all
clc
MeanAAll = 0.1:-0.025:-0.5;
A_STDAll = 0.02:0.02:0.4;

for  MeanAIndex=14:size(MeanAAll,2) % mean of pairwise matrix A
    for A_STDIndex=1:size(A_STDAll,2) % deviation of pairwise matrix A
        for TestIndex=1:3

        % Parameters of the dynamics
        NumberSpecies = 15; % number of samples
        steps = 5000;
        dt = 0.01;

        % Set seed for reproducibility
        % rng(123);
        GrowthRate = normrnd(1,0.2,1,15); % growth rate
        GrowthRate(1) = 1;

        MeanA = MeanAAll(MeanAIndex); % mean of pairwise matrix A
        A_STD = A_STDAll(A_STDIndex); % standard deviation of pairwise matrix A

        % rng(123);
        A = normrnd(MeanA, A_STD, NumberSpecies, NumberSpecies); % Initial interaction matrix

        for i = 1:NumberSpecies
            A(i, i) = -1;
        end

        Collection = [];
        Abundance = [];

        FinalAbundance=[];
        for SynComSize = 1:14
            Combination = nchoosek(2:NumberSpecies, SynComSize);
            if size(Combination, 1) > 30
                Combination = Combination(randsample(1:size(Combination, 1), 30),:);
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
        AbundanceGroup = cell(0);

        for i = 1:15
            AbundanceGroup{i}=FinalAbundance(sum(Collection,2)==i);
        end

        LassoOrder1=Collection(15:end,2:end);
        LassoOrder1(LassoOrder1==0)=-1;
        LassoOrder2=[];
        for i = 1:size(LassoOrder1,1)
            Order2Matrix=nchoosek(LassoOrder1(i,:),2);
            LassoOrder2(i,:)=[LassoOrder1(i,:),(Order2Matrix(:,1).*Order2Matrix(:,2))'];
        end

        LassoOrderY=log10(FinalAbundance(15:end));

        [B1,FitInfo1]=lasso(LassoOrder1, LassoOrderY,'cv',10);

        [B2,FitInfo2]=lasso(LassoOrder2, LassoOrderY,'cv',10);

        Pred1 = sum(LassoOrder1 .* B1(:, FitInfo1.IndexMinMSE)', 2) + FitInfo1.Intercept(FitInfo1.IndexMinMSE);
        Pred2 = sum(LassoOrder2 .* B2(:, FitInfo2.IndexMinMSE)', 2) + FitInfo2.Intercept(FitInfo2.IndexMinMSE);
        validIdx = isfinite(LassoOrderY);

        R2La1=fitlm(Pred1(validIdx),LassoOrderY(validIdx)).Rsquared.Ordinary;
        R2La2=fitlm(Pred2(validIdx),LassoOrderY(validIdx)).Rsquared.Ordinary;
        Lasso1ParaNum=sum(B1(:,FitInfo1.IndexMinMSE)~=0);
        Lasso2ParaNum=sum(B2(:,FitInfo2.IndexMinMSE)~=0);


        R2La2All(MeanAIndex,A_STDIndex,TestIndex)=R2La2;
        Lasso2ParaNumAll(MeanAIndex,A_STDIndex,TestIndex)=Lasso2ParaNum;

        FileName=['Simulation\Results\gLV-',num2str(MeanA),'-',num2str(A_STD),'-',num2str(TestIndex),'.mat'];
        disp(FileName)
        disp(R2La2)
        disp(Lasso2ParaNum)
        save(FileName)
        end
    end
end
