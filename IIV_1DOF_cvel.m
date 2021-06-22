clear fquantiles A
% This code is meant for running several ice speeds. Note that it relies on
% quantiles function from statistics and machine learning toolbox.

totalruns   =	6;             % 1 and >=2 give different visualisations.
tstep       =	0.0002;          %should be under 0.0003.
endtref     =	7.5;             % This is just a reference

% Ice properties
nofelements =	55;
maxelasticdelta     =	-0.002;
k2  =   72000;       
k1  =   k2 * 0.1 * 2;                   
c1  =   k1*3.18*0.65;                      
c2  =   5.2585e+06;    
c2tolerance=0.02*maxelasticdelta*k2;
c2exp=1.2;
spawnmax    =	0.012; %meters
spawnmin    =	0.000; %meters

c2comp= c2^c2exp;   %Corrects linear c2 input. 
c2_initialguess=	c2comp * 0.10549;     %median real c2 damping based on recording c2_his

%Structure properties
struc_k     =	1500000;      %2460000*1
struc_m     =	65;        %52400
struc_c     =	0.005 *2*(struc_k*struc_m)^0.5;

icevweightstr   = 11;
icevfromto  = [-0.0001^(1/icevweightstr) -0.05^(1/icevweightstr)]; %0.0005 - 0.5 
icev        = (icevfromto(1):((icevfromto(2)-icevfromto(1))/(totalruns-1)):icevfromto(2)).^icevweightstr;
endt        = -endtref./(icev-0.09);
ntsteps     = floor(endt/tstep);
sub1        =   min([k2,k1,c1/tstep,c2comp/tstep]);       %%balances matrix A to avoid ill-conditions matrices. use rcond(A) to determine a good value for this. Note that A2 might be more ill-conditioned than A.
t(1:(max(ntsteps)+1))  = (0:tstep:tstep*max(ntsteps));

% Initialize data storage
p_his=nan(max(ntsteps),totalruns);       %%Initialize data storage matrices
v_his=nan(max(ntsteps),totalruns);
f_his=nan(max(ntsteps),totalruns);
c2_his=NaN(max(ntsteps),totalruns);   %%use c2_his to determine mean real c2 e.g. for comparing to linear damping
break_his=NaN(max(ntsteps),totalruns);
initialstruc_p=0.0;
initialstruc_v=0;


parfor run_n = 1:totalruns %% Program repeat loop
    A=[c2_initialguess/tstep,0,0,     sub1,      0;      % another way of avoiding bad scaling would be to convert stuff to N/cm which would have a cubic effect on c2
     0,     k1,     0,      sub1,      sub1;   
     0,     -c1/tstep,    0,      0,      sub1;
     0,     0,      k2,     sub1,      0;
     sub1,     sub1,      sub1       , 0,    0];    
    % A key : 1st col variable is plastic velocity in timestep, 2nd col is
    % viscoelastic compression velocity, 3rd col elastic spring compression,
    % 4th col is the force, 5th col is needed to couple eq rows 2 and 3.
 
 
    struc_v=initialstruc_v;        %Initial structural velocity.
    struc_p=initialstruc_p;     %Initial displacement.
    elements    = [(spawnmax-spawnmin).*rand(nofelements,1)+spawnmin , zeros([nofelements 4])];
%     elements    = zeros(nofelements,5)
    v=icev(run_n);
    nextp_his=NaN;
    nextv_his=NaN;
    runntsteps=ntsteps(run_n);
    p_hisrun=nan(max(ntsteps),1);
    v_hisrun=nan(max(ntsteps),1);
    c2_hisrun=nan(max(ntsteps),nofelements);
    f_hisrun=nan(max(ntsteps),1);

    for it = 1:runntsteps         %%Time integration loop (per timestep in run)
        if ~isnan(nextp_his)
            p_hisrun(it)=nextp_his;
            v_hisrun(it)=nextv_his;        
        end
        
        fstep=zeros(1,nofelements);
        parinputelements=elements;
        for e = 1:nofelements       %%per element loop in timestep in run LOOP
            fcontact=1;
            if (-struc_p+elements(e,1)-elements(e,2)-elements(e,3))>0    % No contact condition, element will move towards equilibrium
                fcontact=0; 
            else
                elements(e,5)=1;      %(e,5) is an activation var. 1 = Contact has been made, element eq of motion is activated.
            end
            
            if elements(e,5) == 1
                A(1,1)=c2_initialguess/tstep;
                B=[0; -k1*elements(e,3); 0 ; 0; sub1*fcontact*(-elements(e,2)-elements(e,3)-struc_p+elements(e,1))];
                X=linsolve(A,B); %solves with c2_initialguess. There is now a X(1) solution that is used in the first while iteration
                c2nonlin=c2_initialguess;
                
                while abs((abs(sub1*X(4)) - abs((c2comp)*(X(1)/tstep)^c2exp))) > -c2tolerance   
                    c2nonlin=((c2comp)*abs((X(1)/tstep)^(c2exp-1))+0.5*c2nonlin)/1.5;     %a weight of 2 for the old value seems to strike a good balance avoiding overshoot, optimizing the speed of finding c2. 2 proved unreliable for large C2, causing overshoot and ill-conditioned matrices -> algorithm gets stuck, use a more conservative weight.

                    A(1,1)=c2nonlin/tstep;
                    X=linsolve(A,B);                    
                    %disp(['Timestepindex ',num2str(it),', Element ',num2str(e),', Force ',num2str(sub1*X(4)),' N, Tolerance ',num2str(abs((abs(sub1*X(4)) - abs((c2comp)*(X(1)/tstep)^3)))),' vs ', num2str(-c2tolerance), ', and X(1) is ', num2str(X(1))])
                    % Use the disp to diagnose what the right c2_tolerance
                    % is as well as for optimizing c2_nonlin guessing
                    % routine. use a breakpoint.
                    
                end
                c2_hisrun(it,e)=c2nonlin;      %use for recording c2
                if X(1)<0                   %condition: dont add any plastic expansion to plastic deformation variable!!!
                    elements(e,2)=elements(e,2)+X(1);
                end
                elements(e,3)=elements(e,3)+X(2); %X(2) is viscoelastic compression/expansion velocity
                if fcontact == 1
                    fstep(e)=-sub1*X(4);     %X(4) is force (same over all 3 element parts)
                else
                    fstep(e)=0;         %No tension forces on structure!
                end
                if X(3) < maxelasticdelta
                    fstep(e)=fstep(e)*((maxelasticdelta-elements(e,4))/(X(3)-elements(e,4)));
                    %break_his(it,run)=struc_p;
                    %elements(e,:)= [(spawnmax-spawnmin).*rand(1,1)+spawnmin+struc_p, zeros([1 4])];
                    elements(e,:)= [(spawnmax-spawnmin).*rand(1,1)+spawnmin+struc_p+(-v*tstep+(X(3)-elements(e,4)))*((maxelasticdelta-elements(e,4))/(X(3)-elements(e,4))) , zeros([1 4])]; %% -v*tstep*(subtimestep timestep portion before failure) compensates for the fact that 
                else
                    elements(e,4)=X(3);     %X(3) is elastic compression 
                end                
            end
            elements(e,1)=elements(e,1)+v*tstep; 
          
        end
        
        %struc_k=(1000*1500*(abs(1000*struc_p)+0.001).^0.85)./(abs(1000*struc_p)+0.001);   %%NONLINEAR STRUCTURE -- This is Secant Modulus.

        struc_a=(sum(fstep)-struc_k*struc_p-struc_v*struc_c)/struc_m;
        struc_v=struc_v+struc_a*tstep;
        struc_p=struc_p+struc_v*tstep;  % You can make this =0, if you want a rigid structure representation
 
        nextp_his=struc_p;
        nextv_his=struc_v;        
        f_hisrun(it)=sum(fstep,'omitnan');

    end
    p_his(:,run_n)=p_hisrun;
    v_his(:,run_n)=v_hisrun;
    f_his(:,run_n)=f_hisrun;
    c2_his(:,run_n)=mean(c2_hisrun,2,'omitnan');

    disp(run_n)
end
        
p_his(1,:)=initialstruc_p;
  
disp(['median c2 / c2comp = ', num2str((nanmedian(c2_his,'all'))/c2comp)]) %%Use this number value to calibrate c2 initialguess. Use this if code is slow or the ratio of A2/B linsolves is more than 5 times A/B
disp(['mean c2 / c2comp = ', num2str((nanmean(c2_his,'all'))/(c2comp))]) %use this value to determine relation between linear and cubic damping (aim for a value of 1 here)   

%% Plots 
throwawayfactor= 0.1    %Some of the beginning data should be thrown away so that means are representative
fmeans=0;
f_fig=figure('Name','Force on structure');

hold on
for run_n = 1:totalruns
    %f_fig=plot(t',movmean(sum(f(:,:,run),2),0.1/tstep),'Color','#0392cf');%Narrow
    %moving mean of force
    fmeans(run_n)=mean(f_his(round(throwawayfactor*ntsteps(run_n)):ntsteps(run_n),run_n),1,'omitnan');
    %f_fig=scatter(-icev(run),mean(sum(f_his(round(throwawayfactor*ntsteps(run)):ntsteps(run),:,run),2,'omitnan'),'omitnan'),'filled','k'); 
    fquantiles(run_n,:)=[quantile(f_his(round(throwawayfactor*ntsteps(run_n)):ntsteps(run_n),run_n),0.25) quantile(f_his(round(throwawayfactor*ntsteps(run_n)):ntsteps(run_n),run_n),0.75)];
end
f_fig=plot(-icev,-fmeans,'-ob','MarkerSize',2);
f_fig=plot(-icev,-fquantiles(:,1),'--b');
f_fig=plot(-icev,-fquantiles(:,2),'--b', 'HandleVisibility', 'off');
legend('Simulated Means','Simulated Quartiles')
xlabel('v_{ice} [m/s]')
ylabel('Global ice load [N]')

%%
f_fig2=figure('Name','Force on structure, the part of the runs that is used for the means');
hold on
for run_n = 1:totalruns
    %f_fig=plot(t',movmean(sum(f(:,:,run),2),0.1/tstep),'Color','#0392cf');%Narrow
    if mod(run_n,1)~=1 % For plotting only a portion of the graphs
        f_fig2=plot(t(round(throwawayfactor*ntsteps(run_n)):ntsteps(run_n))',movmean(f_his(round(throwawayfactor*ntsteps(run_n)):ntsteps(run_n),run_n),1/tstep)); 
    end
end
xlabel('t (s)')
ylabel('moving average (10 seconds) global ice load F (N)')


