%close all
% This code is general code for simulating dynamic interaction between ice
% and some structure. The _cvel version is used for calibration purposes.
% These programs are somewhat cut down from what is used in the thesis, and
% contain only the essentials.

% Simulation parameters
totalruns   =	1;                  % 1 and >=2 give different visualisations.
tstep       =	0.0002;             %should be under 0.0003. Implement Runge-Kutta or Heuns/Midpoint method to increase viable stepsize. Implement sub-timestep ice force for when element fails
endt        =	30;                 

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
iceinitialv =	-0.00             %icev at t=0 (m/s)
iceendv     =	-0.15               %icev at t=endt

c2comp= c2^c2exp; % To get e.g. cube root of C2 for exp 1/3, which is how C2 is used in Hayos code.
c2_initialguess=	c2comp *   0.09164;     %median real c2 damping based on recording c2_his

% Structure properties
struc_k     =	150000;     
struc_m     =	500;        
struc_c     =	0.025   *2*(struc_k*struc_m)^0.5;              %2*(struc_k*struc_m)^0.5/struc_m is equal to critical damping, so give as a fraction of critical damping.

initialstruc_p=0; % use eg. 0.05 to see how strong struc_c is
initialstruc_v=0;

t           =	(0:tstep:endt);
ntsteps     =   endt/tstep;
sub1        =   min([k2,k1,c1/tstep,c2comp/tstep]);       %%balances matrices A and A2 to avoid ill-conditions matrices. use rcond(A) to determine a good value for this. Note that A2 might be more ill-conditioned than A.
toyoyama    =   [0 -iceinitialv*1.2 ;endt -iceendv*1.2];

vincrement=(iceendv-iceinitialv)/(endt/tstep);

% Initialize data storage
%contact_his=nan(length(t),totalruns); %If you want to record contact
%history, this may be used
p_his=nan(length(t),totalruns);       %%Initialize data storage matrices
v_his=nan(length(t),totalruns);
f_his=nan(length(t),nofelements,totalruns);
c2_his=NaN(max(ntsteps),totalruns);   %%use c2_his to determine mean real c2 e.g. for comparing to linear damping
break_his=NaN(length(t),totalruns);
p_his(1,:)=initialstruc_p;

for run = 1:totalruns %% Program repeat loop, parfor also an option
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
    v=iceinitialv;
    nextp_his=NaN;
    nextv_his=NaN;
    c2_hisrun=nan(max(ntsteps),nofelements);
    
    for i = 1:ntsteps          %%Time integration loop (per timestep in run)
        if ~isnan(nextp_his)
            p_his(i,run)=nextp_his;
            v_his(i,run)=nextv_his;        
        end
        
        fstep=zeros(1,nofelements);
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
                
                while abs((abs(sub1*X(4)) - abs((c2comp)*(X(1)/tstep)^c2exp))) > -c2tolerance    %do i need all these abs() to calculate difference
                    c2nonlin=((c2comp)*abs((X(1)/tstep)^(c2exp-1))+0*c2nonlin)/1;     %a weight of 2 for the old value seems to strike a good balance avoiding overshoot, optimizing the speed of finding c2. 2 proved unreliable for large C2, causing overshoot and ill-conditioned matrices -> algorithm gets stuck, use a more conservative weight.
                    A(1,1)=c2nonlin/tstep;
                    X=linsolve(A,B);                    
                    %disp(['Timestepindex',num2str(i),', Element',num2str(e),', Force ',num2str(sub1*X(4)),', Tolerance ',num2str(abs(sub1*X(4)) - abs((c2comp)*(X(1)/tstep)^(c2exp))),' vs ', num2str(-c2tolerance), ', and X(1) is ', num2str(X(1))])
                    % Use the disp to diagnose what the right c2_tolerance
                    % is as well as for optimizing c2_nonlin guessing
                    % routine. use a breakpoint.
                end
                
                %disp(c2nonlin/c2comp)
                c2_hisrun(i,e)=c2nonlin;      %use for recording c2
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
                    break_his(i,run)=struc_p; %%break_his records when and where any ice elements break
                    %elements(e,:)= [(spawnmax-spawnmin).*rand(1,1)+spawnmin+struc_p , zeros([1 4])];
                    elements(e,:)= [(spawnmax-spawnmin).*rand(1,1)+spawnmin+struc_p+(-v*tstep+(X(3)-elements(e,4)))*((maxelasticdelta-elements(e,4))/(X(3)-elements(e,4))) , zeros([1 4])]; %% -v*tstep*(subtimestep timestep portion before failure) compensates for the fact that                     
                else
                    elements(e,4)=X(3);     %X(3) is elastic compression 
                end                
            end
            elements(e,1)=elements(e,1)+v*tstep; 
          
        end
        
        struc_a=(sum(fstep)-struc_k*struc_p-struc_v*struc_c)/struc_m;
        struc_v=struc_v+struc_a*tstep;
        struc_p=struc_p+struc_v*tstep;
        
        v=v+vincrement;
        
        %contact_his(i,run)=nnz(fstep)/nofelements; %record contact area
        %history
        f_his(i,:,run)=fstep;
        nextp_his=struc_p;
        nextv_his=struc_v;

    end
    c2_his(:,run)=mean(c2_hisrun,2,'omitnan');
    disp(run)
end

%% Plotting (p_his and f_his)
disp('Plotting starting')
f_fig=figure('Name','Force on structure');
hold on
for run = 1:totalruns
    f_fig=plot(t',sum(f_his(:,:,run),2),'Color','#0392cf'); %No force plot smoothing
end
f_fig=plot(t',movmean(mean(squeeze(sum(f_his,2)),2),2/tstep),'LineWidth',3,'Color','b');
xlabel('Time [s]');
ylabel('Crushing load [N]');

p_fig=figure('Name','Struc_p plot');
hold on
p_fig=plot(toyoyama(:,1),toyoyama(:,2), '--','LineWidth',1 ,'Color','#878787'); %TOYAMA
p_fig=plot(t',movmean(movmax(abs(v_his),0.6/tstep),0.6/tstep),'LineWidth',1,'Color','#00a9e6');

if totalruns == 1
    p_fig=scatter(t',p_his(:,1),20,sum(f_his(:,:,1),2),'fill');
    colorbar;
    title(colorbar,'F_i_c_e [N]','FontSize',10);
else
    p_fig=plot(t',p_his,'Color','k');    
end
xlabel('t [s]')
ylabel('u [m]')
legend('Toyama law 120%', 'Structure max velocity', 'Structure position')

disp(['median c2 / c2comp = ', num2str((nanmedian(c2_his,'all'))/c2comp)]) %%Use this number value to calibrate c2 initialguess. Use this if code is slow or the ratio of A2/B linsolves is more than 5 times A/B
%%  Contact history plot
% yyaxis right;
% plot(t,movmean(contact_his,100))
% ylim([0 1])