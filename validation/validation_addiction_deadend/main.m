%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Main file to run the opioid model %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%format long e
% site on linear chain trick
% http://sherrytowers.com/2019/02/19/contagion-models-with-non-exponentially-distributed-sojourn-times-in-the-infectious-state/
figure()

%% Input parameters
sigavg = 0.00744; 
range = 0;

class_vec = [1 5 10];


nc = length(class_vec);

for classy = 1:nc
    
    classes = class_vec(classy);
    p = parameters(classes,range,sigavg);
    
    %% Define initial conditions
    
    inits = zeros(p.num_classes,1); 
    
    prescribe_init = 0.0001; % total percent prescribed
    for i = 1:p.k
        inits(i) = prescribe_init*p.N/p.k; % prescribed, P, divide equally between all prescribed classes
    end
    inits(p.k+2) = 0;%0.05*p.N;    % addicted
    inits(p.k+1) = p.N - sum(inits);% susceptible, S
    
    
    %% Run the model

   options = odeset('NonNegative',1:p.num_classes); 
   [t, pop]=ode45(@opioid,[0 p.endT],inits,options,p);
    
    
    %% Extract population vectors
    
    S = pop(:,p.k+1)./p.N*100;
    A = pop(:,p.k+2)./p.N*100;
    
    P = pop(:,1)./p.N*100; % total prescribed
    if p.k > 1
        for i = 2:p.k
            P = P + pop(:,i)./p.N*100;
        end
    end
    
    cumulP = pop(:,end-1);
    cumulA = pop(:,end);

    %% Plotting fun
 if classes == 1        
        subplot(2,2,1); plot(t,P,'k','LineWidth',2)
        ylabel('Percent Prescribed')
        axis([0 p.endT 0 10])
         hold on
        subplot(2,2,2); plot(t,A,'k','LineWidth',2)
        ylabel('Percent Addicted')
        axis([0 p.endT 0 10])
        hold on
        subplot(2,2,3); plot(t,S,'k','LineWidth',2)
        ylabel('Percent Susceptible')
        xlabel('Time (years)')
        axis([0 p.endT 85 100])
        hold on
        subplot(2,2,4); plot(t,cumulP,'k','LineWidth',2)
        ylabel('Cumulative Number Prescribed')
        xlabel('Time (years)')
        axis([0 p.endT 0 max(cumulP)])
        hold on

        hold on

       
     else
        subplot(2,2,1); plot(t,P,'--k')
         hold on
        subplot(2,2,2); plot(t,A,'--k')
        subplot(2,2,3); plot(t,S,'--k')
        subplot(2,2,4); plot(t,cumulP,'--k')
        
     end
    

    

end

 %  figure();
 % plot(t,S+P+A)