clear
close all

addpath utils/

plotfig = 0; % 1 to plot figures as the run happens (this is slower), 0 to not plot

if plotfig
    F1 = figure(1);
    set(F1, 'Position',  [100, 100, 1000, 1600])
    figrun = 1;
    stemc = 'b';
end


platform = "GPU";
%restart the simulation from the last event, 0 if no, 1 if yes. Parameters
%have to be consistent, doing a restart overwrites all the previous
%parameters set in the simulations so you may be running a different
%simulation than you think

restart = 0;



%% Input parameters

input_table = readtable("BP4input_mat.txt","Delimiter",':',ReadRowNames=true);
sample_points = readtable("BP4_sample_points.txt","Delimiter",' ',ReadRowNames=true);

sample_point_id = 8;

% onfault_data = readtable(['./jiang/' num2str(sample_point_id) '.csv'],"Delimiter",' ',ReadRowNames=true);
% 
% global_data = readtable('./jiang/others/jiang.csv',"Delimiter",' ',ReadRowNames=true);
% 
% 
% benchmark_time = onfault_data.t;

% General inputs


cs = input_table.Var1("cs");
rho = input_table.Var1("rho");
fract = input_table.Var1("fract"); % used in frac to scale accuracy, a larger value means more accurate and shower
frac = 1/2^(fract-1);
tolup = frac*input_table.Var1("tolup");
tollo = frac*input_table.Var1("tollo");
AL = input_table.Var1("AL");  % select law. 1: Ageing, else: Slip
ac = input_table.Var1("a");
aRSc = input_table.Var1("aRS"); %a value to create rate strengthening
bc = input_table.Var1("b");
Dc = input_table.Var1("Dc");
fr = input_table.Var1("fr");
si0 = input_table.Var1("si0");
nu = input_table.Var1("nu");
miniter = input_table.Var1("miniter"); %should not be changed
tf = input_table.Var1("tf"); % End time of simulation in years

G = cs^2*rho;
eta = G/(2*cs); % Radiation damping coefficient
Catalog = zeros(10000,13);

% Velocities

Vpl = input_table.Var1("Vp"); % Plate slip rate
V0  = input_table.Var1("Vi"); % Initial slip rate
Vr  = input_table.Var1("Vr"); % Reference slip rate
Vi  = 0.001; % Slip rate for fav zone


% Domain
W = input_table.Var1("Wf"); % Half width of the fault
H = input_table.Var1("H"); % half width of the VW domain
L = input_table.Var1("L"); % half length of the VW domain
h = input_table.Var1("h"); % buffer zone dimension

gridside = input_table.Var1("gridside");

% Nucleation zone

wi = input_table.Var1("wi"); % width and height of the nucleation zone
hi = input_table.Var1("hi"); % distance of the nucleation zone from VW boundary




% find minimum timestep, the aim is that shear wave speed would be
% resolved.
if frac <= 1
    mindt = frac/cs*2*gridside;
else
    mindt = 1/cs*2*gridside;
end
dt = mindt;
% mindt = 0.;
% dt = 1e-9;


%% Check important dimensions

Lb = G/(1-nu)*Dc/(bc*si0); %cohesive zone length, should be resolved

if gridside*3>Lb
    disp('Cohesive zone may be poorly resolved. Stopping')
    return
end

Linf = pi/4*Lb*(bc/(bc-ac))^2;

if gridside*10>Linf
    disp('Linf is poorly resolved. Stopping')
    return
end

if H/Linf < 1
    disp('VW domain smaller then Linf. Stopping')
    disp("R/Linf: "+num2str(H/Linf))
    return
end

disp("R/Linf: "+num2str(H/Linf))


%% Initialize geometry

domain_size = 1.2e5; % Total domain size in m (120km)
half_gridside = gridside/2;

domainpower = round(log2((domain_size/2)/half_gridside)); %round the result so that we have a nice power of 2

L_domain = half_gridside*2^domainpower;
W_domain = half_gridside*2^domainpower;

X = -W_domain:gridside:W_domain-1;
Y = -W_domain:gridside:W_domain-1;
% X = linspace(-W_domain, W_domain, 2^domainpower+1);
% X = X(1:end-1);
% Y = linspace(-W_domain, W_domain, 2^domainpower+1);
% Y = Y(1:end-1);

[x,y] = meshgrid(X,Y);
N = 2^domainpower;

sample_point_index = find(x == sample_points{sample_point_id,1}*1000 & y == sample_points{sample_point_id,2}*1000);



dLo = (abs(y) <= W);  % inside loading area 
dRW = (abs(x) < L) .* (abs(y) < H); %inside rate/velocity weakening area
dRS = ((abs(x) >= L+h) + ((abs(y) <= W_domain) .*(abs(y) >= H+h)))>=1;
dTR = (dRS-dRW) == 0;

dNU = ((x>(-L+hi)).* (x < (-L+wi+hi)) .* (y>(-H+hi)).* (y < (-H+wi+hi)))==1;

dMR = (abs(x)<=L+2*h) .* (abs(y)<=H+2*h); % area used to calculate the moment rates
NNodes = sum(sum(dMR));



dCR = dLo == 0; %creeping region
dFD = dNU == 0; %everything out of the nucleation zone

 
%% Initial values and matrices

% Calculate values for a in the transition zone
x_buffer = x(dTR);
y_buffer = y(dTR);

r = (max([abs(y(dTR))-H,abs(x(dTR))-L], [], 2))./h;

template = zeros(size(dLo));
a_buffer = ac + r*(aRSc-ac);
a = ac + template;
a(dRS) = aRSc;
a(dTR) = a_buffer;

% for ix = 1:length(x(1,:))
%     for ixx = 1:length(x(:,1))
%         if dTR(ix,ixx)
%         ra = (max([abs(y(ix,ixx))-H,abs(x(ix,ixx))-L]))/h;    
%         a(ix,ixx) = ac + ra*(aRSc-ac);   
%         end
%     end
% end


b = bc + template;

%Intial slip
dx = template; % slip along the x

% Initial slip rate
V = V0+template; % slip rate along the x
V(dNU) = Vi; % We impose higher slip rate in the favorable nucleation zone

% Initial state value
theta = Dc/V0+template;


% Initial tau variable
tau0 = si0.*a.*asinh((V0/(2*Vr))*exp((fr+bc*log(Vr/V0))./a))+eta*V0;
tau0(dNU) = si0*ac.*asinh((Vi/(2*Vr))*exp((fr+bc*log(Vr/V0))./ac))+eta*Vi;

%% Fourier
% Initialize Fourier coefficienents
F = fft2(dx);
Fs = 1/gridside;
frex = Fs*(-N/2:1:N/2-1)/N;
frey = Fs*(-N/2:1:N/2-1)/N;
[frex,frey] = meshgrid(frex,frey);
kvx = ifftshift(2*pi*frex); %!
kvy = ifftshift(2*pi*frey); %!



% Pre-compute constants outside the loop to save time:

antiplane = G;
inplane = G/(1-nu);
denom = sqrt((kvx).^2 + (kvy).^2);

akdhatx=   -inplane/2*(kvx).^2./denom;
akdhaty=   -antiplane/2*(kvy).^2./denom;


akdhatx(isnan(akdhatx)) = 0;
akdhaty(isnan(akdhaty)) = 0;
akdhat = (akdhatx + akdhaty);


%% Main Loop Starts
tryagain = 0;


counter = 0;
maxtolvio = 0;
cumutolvio = 0;


t = 0; % time in seconds
time = 0; % time in years
timep = 0;

eventrun = 1;
eventstart = 0;
err = 1.;
erra= 1.;

timep = 0;
time = 0;


it = 1;

Vcontour = template+NaN;
first_contour = 0;
t_fc = 0;


if restart
    load('BP4_QD_state');  
else
    itstart = 2; 
    save('BP4_QD_state')
    save('BP4_QD_Catalog','Catalog')
end

sliprates = zeros(1,70000);
max_sliprates = zeros(1,70000);
momentrates = zeros(1,70000);
slips = zeros(1,70000);
states = zeros(1,70000);
shears = zeros(1,70000);
times = zeros(1,70000);
dts = zeros(1,70000);
errs = zeros(1,70000);
erras = zeros(1,70000);


if platform == "GPU"
    dx = gpuArray(dx);
    V  = gpuArray(V);
    theta = gpuArray(theta);
    akdhat = gpuArray(akdhat);
    tau0 = gpuArray(tau0);
    a = gpuArray(a);
    b = gpuArray(b);
    dLo = gpuArray(dLo);
    dCR = gpuArray(dCR);
end

% NT = 10000;

while t <= tf*(365*24*60*60)
% while it <= NT
    dxp = dx;
    Vp = V;
    thetap = theta;

    if err < tollo && erra < 10*tollo
        dt = dt*1.1;
    end

    if dt < mindt
        dt = mindt;
    end

    for tryagaincount=0:9

        % If this is the first attempt
        if tryagaincount < 1

            Vg =  V;
            Vmg = V;
            dxg = dx + dt*V;

        elseif tryagaincount == 1
            % First iteration is treated as a special case where time-step 
            % is not decreased, in attempt to optain a better guess.

            dxg = dx;
            Vg = V;
            Vmg = Vm;

        else

            dxg = 0.5*(dx + dxp);
            Vg = 0.5*(Vp + V);
            Vmg = 0.5*(Vp + Vg);

        end

        F =  fft2(dxg);

        taur =  real(ifft2(akdhat.*F));

        tau = tau0+taur;


        % guess the state variable

        thetag = calc_theta(thetap, Vmg, dt, Dc, AL);


        V = calc_sliprate(Vg, tau, si0, a, eta, Vr, fr, b , thetag, Dc);
        V = V.*dLo + Vpl.*dCR;
        Vm = 0.5*(Vp + V);

        s = dt*Vm;
        dx = dxp + s;


        if tryagaincount < miniter
            err = 1;  
            erra = 1;
        else
            [err, erra] = calc_err(V, Vg);
        end

        if (tryagaincount < miniter || ((err > tolup || erra > 10*tolup) && dt > mindt))
            if tryagaincount>0
                % If accuracy metrics are violated then time-step is refined
                dt = dt/2;
            end

        else

           maxV = max(V, [], "all");

           if plotfig

               if maxV>1e-3
                  plotrat = 1; %plot less frequently at high speeds, increase
               else
                  plotrat = 2; %plot less frequently at low speeds, increase
               end

               if rem(it,plotrat*100) == 0 || it < 3

                   set(0,'CurrentFigure', F1)
                   subplot(3,2,1)
                   plotmat(log10(V))
                   title('log_{10} slip speed [m/s]')
                   caxis([-15 1])
                   subplot(3,2,2)
                   plotmat(dx)
                   title('Cumulative slip [m]')
                   subplot(3,2,3)
                   plotmat(tau/1000000)
                   title('Shear stress [MPa]')

                   if first_contour == 1
                       subplot(3,2,4)
                       [c,h] = contour(Vcontour, -10:5:30);
                       % clabel(c,h);
                       colorbar
                       title('Time to Vth > 0.1m/s')
                       daspect([1 1 1])
                   end
                   % if eventrun == 1
                   %     line(benchmark_time, benchmark_data.slip_2)
                   % end
                   % hold on
                   % scatter(t, dx(sample_point_index), 'black','filled')
                   % hold off

                   % hold on
                   % scatter(t, log10(V(sample_point_index)))
                   % hold off
                   if eventrun > 2
                       subplot(3,2,5:6)
                       stem(Catalog(1:eventrun-1,4),Catalog(1:eventrun-1,6),'filled',stemc)
                       ylabel('Magnitude')
                       xlabel('Years')
                       ylim([0, 9])
                   end
                   sgtitle(strcat('---total time:',string(seconds(t))))
                   drawnow
               end
           figrun = figrun + 1;
           end

           % identify that an even has started and gather information
           if maxV > 1e-3 && eventstart == 0
              eventstart = 1;

              TEMPD = dx;
              TEMPS = tau;
              ImaxV = find(V > 1e-3);
              meanx = mean(x(ImaxV));
              meany = mean(y(ImaxV));
              timep = time;
              timesecstart = t;
              disp('Event has started')
           end

           if eventrun == 1 

                Vth = (V >= 0.1) & isnan(Vcontour);
                if (maxV >= 0.1) && first_contour == 0
                    first_contour = 1;
                    t_fc = t;
                end
                Vcontour(Vth) = t-t_fc;

           end
           % indentify that event has ended and compute various properties
           if maxV <= 1e-3 && eventstart == 1 && t-timesecstart>10
              eventstart = 0;
              timesecend = t;
              time = timesecend/(365*24*60*60);
              Trup = timesecend - timesecstart;
              Slip = dx - TEMPD;
              DTAU = tau - TEMPS;
              Ievent = dx - TEMPD > Dc*0.5;
              DX = max(max(x(Ievent)))-min(min(x(Ievent)));
              DY = max(max(y(Ievent)))-min(min(y(Ievent)));
              Nevent = sum(sum(Ievent));
              Area = Nevent*gridside^2;
              MeanSlip = sum(sum(Slip(Ievent)))/Nevent;
              MeanStress = sum(sum(DTAU(Ievent)))/Nevent;
              Momen = G*MeanSlip*Area;
              Mag = (log10(Momen) - 9.05)/1.5;
              Catalog(eventrun,:) = [eventrun meanx meany time Momen Mag time-timep Area DX DY MeanStress MeanSlip Trup];
              eventrun = eventrun + 1;
              itstart = it;

              % clf
              save('BP4_QD_state')
              save('BP4_QD_Catalog','Catalog')
              disp('state saved, Catalog,saved')
              disp(strcat('Event nr. ',num2str(eventrun-1),' has ended, magnitude = ',num2str(Mag),' at time = ',num2str(time),'years, intereventtime = ',num2str(time-timep),' years'))
           end

           break

        end       

    end

    t = t + dt;

    theta = calc_theta(thetap, Vm, dt, Dc, AL);
    sliprates(it) = V(sample_point_index);
    max_sliprates(it) = log10(maxV);
    meansliprate = sum(sum(V.*dMR))/NNodes;
    momentrates(it) = G*meansliprate*NNodes*gridside^2;
    slips(it) = dx(sample_point_index);
    states(it) = theta(sample_point_index);
    shears(it) = tau(sample_point_index);
    times(it) = t;
    dts(it) = dt;
    errs(it) = err;
    erras(it) = erra;



    

    it = it+1;
end



save output/out_GPU_500.mat sliprates slips states shears times Vcontour errs erras dts
% save("slip_rates.mat","sliprates")
% save("times.mat","times")

% subplot(2,3,3)
% title('Shear stress comparison')
% 
% ylabel('Shear stress [MPa]')
% xlabel('Time (s)')
% 
% line(benchmark_time, onfault_data.shear_stress_2,'b-')
% hold on
% plot(times, shears/1e6,'black')
% hold off
% subplot(2,3,4)
% title('State comparison')
% ylabel('Log10 state')
% xlabel('Time (s)')
% 
% line(benchmark_time, onfault_data.state,'b-')
% hold on
% plot(times, log10(states),'black')
% hold off
% 
% subplot(2,3,5)
% title('Max slip rate comparison')
% ylabel('Log10 max slip rate')
% xlabel('Time (s)')
% 
% line(global_data.t, global_data.max_slip_rate)
% hold on
% plot(times, max_sliprates,'black')
% hold off
% 
% subplot(2,3,6)
% title('Moment rate comparison')
% ylabel('Moment rate')
% xlabel('Time (s)')
% 
% line(global_data.t, global_data.moment_rate)
% hold on
% plot(times, momentrates,'black')
% hold off

% while t <= tf*(365*24*60*60)
% 
%     it = it+1;
%     tryagaincount = 0;
%     enter = 1;
%     %% guess    
% 
%     Vg =  V;
%     Vmg = 0.5*(V + Vg);
%     dxg = dx + dt*Vmg;
% 
% 
% 
%     dxp = dx;
%     thetap = theta;
% 
%     Vp = V;
% 
%     while tryagain || enter
%         counter = counter + 1;
% 
%         % If this is the first attempt
%         if enter == 1
% 
%             F = fft2(dxg);
% 
% 
%         elseif tryagaincount > 1
% 
%             dxg = 0.5*(dx + dxp);
% 
%             F =  fft2(dxg);
% 
%             Vg = 0.5*(Vp + V);
%             Vmg = 0.5*(Vp + Vg);
%         elseif tryagaincount == 1
%             % First iteration is treated as a special case where time-step 
%             % is not decreased, in attempt to optain a better guess.
% 
%             dxg = dx;
% 
%             F =  fft2(dxg);
% 
%             Vg = V;
%             Vmg = 0.5*(Vp + Vg);
% 
%         end
% 
%         %if enter
% 
%         taur =  real(ifft2(akdhat.*F));
% 
%         %end
% 
%         enter = 0;
% 
%         maxV = max(max(Vg));
%         %guess the state variable
%         if AL == 1
%             thetag = thetap.*exp(-Vmg*dt/Dc) + Dc./Vmg.*(1 - exp(-Vmg*dt/Dc));
%         else
%             thetag = Dc./Vmg.*(Vmg.*thetap/Dc).^(exp(-Vmg*dt/Dc));
%         end
% 
% 
% 
%         tau = tau0 + taur;
% 
%         if maxV>0.00001
%          %use linearized method   
%          AA = tau./si0./a;
%          BB = eta./si0./a; %not special meaning just a combination that occurs in the linearized update
%          CC = 2*Vr./( exp((fr + b .* log(Vr * thetag / Dc)) ./ a) );
%          dV = -(Vg - CC.*sinh(AA-BB.*Vg))./(BB.*CC.*cosh(AA-BB.*Vg) + 1);
%          V = dV + Vg;
%         else
%          %pure explicity
%          V = 2 * Vr * sinh((tau - eta * Vg) ./ (si0) ./ a) ./ exp((fr + b .* log(Vr * thetag / Dc)) ./ a);      
%         end
% 
%         V = V.*dLo + Vpl.*dCR;
%         Vm = 0.5*(Vp + V);
% 
%         % Aging/Slipping law for state variable update
%         if AL == 1
%             theta = thetap.*exp(-Vm*dt/Dc) + Dc./Vm.*(1 - exp(-Vm*dt/Dc));
%         else
%             theta = Dc./Vm.*(Vm.*thetap/Dc).^(exp(-Vm*dt/Dc));
%         end
%         maxV = max(max(V));
%         dx = dxp + dt*Vm;
%             if tryagaincount < miniter
%             err = 1;  
%             erra = 1;
%             else
%             err = norm(V-Vg,1)/(norm(V,1));
%             [erra,ierra] = max(abs(V(:)-Vg(:)));
%             erra = erra/V(ierra);
%             end
%         if ((tryagaincount < miniter || ((err > tolup || erra > 10*tolup) && dt > mindt)) && tryagaincount < 10 )
%             tryagain = 1;
%             tryagaincount = tryagaincount + 1;
%             if tryagaincount==1
%                 % First iteration refines the prediction and does not reduce the time-step
%                 % dt = dt;
%             else
%                 % If accuracy metrics are violated then time-step is refined
%                 dt = dt/2;
%             end
% 
%         else
%            tryagain = 0;
%            tryagaincount = 0; 
%            if plotfig
%            if maxV>0.001 
%               plotrat = 1; %plot less frequently at high speeds, increase
%            else
%               plotrat = 1; %plot less frequently at low speeds, increase
%            end
%            if rem(it,plotrat*100) == 0 || it < 3
% 
%            hhh = figure(1);
%            subplot(3,2,1)
%            plotmat(log10(V))
%            title('log_{10} slip speed [m/s]')
%            daspect([1 1 1])
%            colorbar
%            caxis([-15 1])
%            subplot(3,2,2)
%            plotmat(dx)
%            title('Cumulative slip [m]')
%            daspect([1 1 1])
%            colorbar
%            subplot(3,2,3)
%            plotmat(tau/1000000)
%            daspect([1 1 1])
%            % caxis([tau0 - 4.0e6 tau0 + 4.0e6]/1000000)
%            title('Shear stress [MPa]')
%            colorbar
%            if eventrun > 2
%            subplot(3,2,4)
%            CUM = fliplr(cumsum(ones(1,eventrun-2)));
%             MAGS = sort(Catalog(2:eventrun-1,6));
%             NMle0 = sum(MAGS >= 0);
%             plot(MAGS,log10(CUM/NMle0))
%             hold on
%             plot(MAGS,log10(10.^(-1*MAGS)))
%             plot(MAGS,log10(10.^(-0.75*MAGS)),'-.')
%             plot(MAGS,log10(10.^(-1.25*MAGS)),'--')
%             ylabel('log_{10}(N/N(M>0))')
%             hold off
%            subplot(3,2,5:6)
%            stem(Catalog(1:eventrun-1,4),Catalog(1:eventrun-1,6),'filled',stemc)
%            ylabel('Magnitude')
%            xlabel('Days')
%            ylim([0 7])
%            end
%            sgtitle(strcat('---total time:',datestr(seconds(t),'yy-mm-dd HH-MM-SS.FFF')))
%            drawnow
%            end
%            figrun = figrun + 1;
%            end
% 
% 
% 
%            % identify that an even has started and gather information
%            if maxV > 0.01 && eventstart == 0
%               eventstart = 1;
% 
%               TEMPD = dx;
%               TEMPS = tau;
%               ImaxV = find(V >= 0.001);
%               meanx = mean(x(ImaxV));
%               meany = mean(y(ImaxV));
%               timep = time;
%               time = t/(24*60*60);
%               timesecstart = t;
%               disp('Event has started')
%            end
% 
%           if eventrun <= 2 %&& (t == 0 || t==10 || t == 20 || t==30)
%                 subplot(3,2,4)
%                 Vth = (V >= 0.001) & isnan(Vcontour);
%                 Vcontour(Vth) = t;
%            end
%            %indentify that event has ended and compute various properties
%            if maxV <= 0.001 && eventstart == 1
%               eventstart = 0;
%               timesecend = t;
%               Trup = timesecend - timesecstart;
%               Slip = dx - TEMPD;
%               DTAU = tau - TEMPS;
%               Ievent = dx - TEMPD > 3*Dc;
%               DX = max(max(x(Ievent)))-min(min(x(Ievent)));
%               DY = max(max(y(Ievent)))-min(min(y(Ievent)));
%               Nevent = sum(sum(Ievent));
%               Area = Nevent*gridside^2;
%               MeanSlip = sum(sum(Slip(Ievent)))/Nevent;
%               MeanStress = sum(sum(DTAU(Ievent)))/Nevent;
%               Momen = G*MeanSlip*Area;
%               Mag = (log10(Momen) - 9.05)/1.5;
%               Catalog(eventrun,:) = [eventrun meanx meany time Momen Mag time-timep Area DX DY MeanStress MeanSlip Trup];
%               eventrun = eventrun + 1;
%               toc
%               itstart = it;
%               clf
%               % disp(it/NT)
%               save('state')
%               save('Catalog','Catalog')
%               disp('state saved, Catalog,saved')
%               disp(strcat('Event nr. ',num2str(eventrun-1),' has ended, magnitude = ',num2str(Mag),' at time = ',num2str(time),'days, intereventtime = ',num2str(time-timep),' days'))
%            end
% 
%         end       
% 
%     end
% 
% 
%     sliprates(it) = V(64, 64);
%     times(it) = t;
% 
%     t = t + dt;
%     dtp = dt;
%     %dt = frac*L./maxV;% 1.0e9*frac*min(theta)/8]); 
% 
%     if err < tollo && erra < 10*tollo
%        dt = (2^(miniter-1) + 0.1)*dtp;
%     else
%        dt =  2^(miniter-1)*dtp;
%     end
% 
%     if dt < mindt
%        dt = mindt;
%     end
% 
% end
% 
