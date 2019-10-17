% Code belonging to figure 3 of Ma et al. (2011), Behavior and neural basis
% of near-optimal visual search, Nature Neuroscience 14 (6), 783-90

% Generates distributions of decision variable and ROCs
% Change line 20 to switch between homogeneous and heterogeneous
% distractors

clear all; close all;
set(0,'DefaultLineLineWidth',2);
sT = 10;
sD = 0;
Ntrials = 100000;
sigma_high = 2;
sigma_low = 6;
sigma_inference = 4;
kappa_high = 10;
kappa_low = 5;
kappa_inference = 5;

experiment = 1; % 1 for homogeneous distractors, 2 for heterogeneous distractors

switch experiment
    case 1 % homogeneous distractors
        Nvec = [4 6 8 ];

        for Nind = 1:length(Nvec)
            N = Nvec(Nind);
            high = (rand(Ntrials,N)>0.5);
            xA = zeros(Ntrials,N);
            xP = zeros(Ntrials,N);


            sigma = sigma_high * high + sigma_low * (1-high);

            % Generating distractors
            xA(high) = sD + sigma_high * randn(sum(sum(high)),1);
            xA(~high) = sD + sigma_low * randn(sum(sum(~high)),1);
            xP(high) = sD + sigma_high * randn(sum(sum(high)),1);
            xP(~high) = sD + sigma_low * randn(sum(sum(~high)),1);

            % Generating targets
            xP(:,1) = sT + sigma(:,1) .* randn(Ntrials,1);

            % Local log likelihood ratio
            diA = (sT-sD) * (xA-(sT+sD)/2)./sigma.^2;
            diP = (sT-sD) * (xP-(sT+sD)/2)./sigma.^2;        

            % Global log likelihood ratio
            dA = log(mean(exp(diA),2));
            dP = log(mean(exp(diP),2));

            % Max_x model
            dA_max = max(xA,[],2);
            dP_max = max(xP,[],2);

            % Single-reliability (1r) model
            diA_sr = (sT-sD) * (xA-(sT+sD)/2)/sigma_inference^2;
            diP_sr = (sT-sD) * (xP-(sT+sD)/2)/sigma_inference^2;        
            dA_sr = log(mean(exp(diA_sr),2));
            dP_sr = log(mean(exp(diP_sr),2));

            % Estimate distributions and ROC
            [h_old c] = hist([dA dP],500);
            h(:,:,Nind) = h_old/Ntrials;
            [h_max_old c_max] = hist([dA_max dP_max],100);
            h_max(:,:,Nind) = h_max_old/Ntrials;
            [h_sr_old c_sr] = hist([dA_sr dP_sr],100);
            h_sr(:,:,Nind) = h_sr_old/Ntrials;
        end


        for Nind = 1:length(Nvec)
            figure; hold on;
            plot(c,h(:,1,Nind),'r');  
            plot(c,h(:,2,Nind),'Color',[0 0.5 0]); 
            xlim([-15 25])
            set(gca,'ytick',[]);
            xlabel('Log likelihood ratio'); ylabel('Probability')
            title(strcat('N=', int2str(Nvec(Nind))));
            
            figure; hold on;
            plot(cumsum(h_max(end:-1:1,1,Nind)),cumsum(h_max(end:-1:1,2,Nind)),'Color',[0 1 1]); axis([0 1 0 1]);
            plot(cumsum(h_sr(end:-1:1,1,Nind)),cumsum(h_sr(end:-1:1,2,Nind)),'Color',[.5 .25 .25]); axis([0 1 0 1]);
            plot(cumsum(h(end:-1:1,1,Nind)),cumsum(h(end:-1:1,2,Nind)),'b'); axis([0 1 0 1]);
            plot([0 1], [0 1],'k--')
            set(gca,'xtick',[0 0.5 1]); set(gca,'ytick',[0 0.5 1]); box on; 
            xlabel('False-alarm rate'); ylabel('Detection rate');
            legend('Max_x','1r','Optimal','Location','Best');
            title(strcat('N=', int2str(Nvec(Nind))));
        end

      figure; hold on;
      plot(cumsum(squeeze(h(end:-1:1,1,1))),squeeze(cumsum(h(end:-1:1,2,1))),'r'); 
      plot(cumsum(squeeze(h(end:-1:1,1,2))),squeeze(cumsum(h(end:-1:1,2,2))),'b'); 
      plot(cumsum(squeeze(h(end:-1:1,1,3))),squeeze(cumsum(h(end:-1:1,2,3))),'g'); 
      axis([0 1 0 1]);
      plot([0 1], [0 1],'k--')
      xlabel('False-alarm rate'); ylabel('Detection rate');
      set(gca,'xtick',[0 0.5 1]); set(gca,'ytick',[0 0.5 1]); box on; 
      legend('N=4', 'N=6', 'N=8');


    case 2 % heterogeneous distractors
        N = 4;
        high = (rand(Ntrials,N)>0.5);
        xA = zeros(Ntrials,N);
        xP = zeros(Ntrials,N);
        
        sT = sT * pi/180; % everything in radians
        kappa = kappa_high * high + kappa_low * (1-high);
        
        % Generating distractors
        xA = pi * rand(Ntrials,N);
        xP = pi * rand(Ntrials,N);

        % Generating targets
        xP(high(:,1),1) = sT + randraw('vonmises', [0 kappa_high], sum(high(:,1)),1)/2;        
        xP(~high(:,1),1) = sT + randraw('vonmises', [0 kappa_low], sum(~high(:,1)),1)/2;        
        
        % Local log likelihood ratio
        diA = -log(besseli(0,kappa)) + kappa .* cos(2*(xA-sT));
        diP = -log(besseli(0,kappa)) + kappa .* cos(2*(xP-sT));
        
        % Global log likelihood ratio
        dA = log(mean(exp(diA),2));
        dP = log(mean(exp(diP),2));
        
        % Max_x model
        dA_max = max(xA,[],2);
        dP_max = max(xP,[],2);

        % Single-reliability (1r) model
        diA_sr = -log(besseli(0,kappa_inference)) + kappa_inference * cos(2*(xA-sT));
        diP_sr = -log(besseli(0,kappa_inference)) + kappa_inference * cos(2*(xP-sT));
        dA_sr = log(mean(exp(diA_sr),2));
        dP_sr = log(mean(exp(diP_sr),2));
        
        % Conditioning on target contrast
        I = (high(:,1)==0);
        J = (high(:,1)==1);
        
        % Estimate distributions and ROC
        [h c] = hist([dA dP],100);
        h = h/Ntrials;
        [h_max c_max] = hist([dA_max dP_max],100);
        h_max = h_max/Ntrials;
        [h_sr c_sr] = hist([dA_sr dP_sr],100);
        h_sr = h_sr/Ntrials;
        figure; hold on;
        plot(c,h(:,1),'r');  
        plot(c,h(:,2),'Color',[0 0.5 0]); 
        xlim([-10 2]); set(gca,'ytick',[]); 
        xlabel('Log likelihood ratio'); ylabel('Probability')

        figure; hold on;
        plot(cumsum(h_max(end:-1:1,1)),cumsum(h_max(end:-1:1,2)),'Color', [0 1 1]); axis([0 1 0 1]);
        plot(cumsum(h_sr(end:-1:1,1)),cumsum(h_sr(end:-1:1,2)),'Color', [.5 .25 .25]); axis([0 1 0 1]);
        plot(cumsum(h(end:-1:1,1)),cumsum(h(end:-1:1,2)),'b'); axis([0 1 0 1]);
        plot([0 1], [0 1],'k--')
        set(gca,'xtick',[0 0.5 1]); set(gca,'ytick',[0 0.5 1]); box on; 
        xlabel('False-alarm rate'); ylabel('Detection rate');    
        legend('max_x','1r','Optimal','Location','Best');
        title('Comparison between models')

        [h_all c_all] = hist([dA dP],200);
        [h_low c_low] = hist([dA(I) dP(I)],200);
        [h_high c_high] = hist([dA(J) dP(J)],200);
        h_all = h_all/Ntrials;
        h_low = h_low/sum(I);
        h_high = h_high/sum(J);
        
        figure; 
        subplot(3,1,1); plot(c_low,h_low); title('Low, All, High'); ylabel('Probability')
        subplot(3,1,2); plot(c_all,h_all); ylabel('Probability')
        subplot(3,1,3); plot(c_high,h_high);xlabel('Log likelihood ratio'); ylabel('Probability')
        
        figure; hold on;
        plot(cumsum(h_low(end:-1:1,1)),cumsum(h_low(end:-1:1,2)),'r'); 
        plot(cumsum(h_all(end:-1:1,1)),cumsum(h_all(end:-1:1,2)),'b'); 
        plot(cumsum(h_high(end:-1:1,1)),cumsum(h_high(end:-1:1,2)),'g'); 
        plot([0 1], [0 1], 'k--')
        legend('Low','All','High')
        xlabel('False-alarm rate'); ylabel('Detection rate');
        axis([0 1 0 1]);
        set(gca,'xtick',[0 0.5 1]); set(gca,'ytick',[0 0.5 1]); box on;         
        title('Conditioned on target contrast')
end