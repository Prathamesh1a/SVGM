clc
clear
%% INPUT PARAMETERS
filename= 'Time-History.xlsx';%Input time history at station 1
% interval=[0,5.135,10.75,14.305,18.91,25.965,34.99,39.015]; %Input time history segments H2
interval=[0,5.275,12.76,20.075,24.35,27.135,30.36,35.335,39.01]; %Input time history segments H1
D = [0,50,100,150,200,250];%Input distances
wk=[15,15,15,15,15,15];%Filter Frequency at each site
zk=[0.6,0.6,0.6,0.6,0.6,0.6];%Soil Damping at each site
alpha=0.2; %Incoherence coefficient
Vs=400; %(m/s)

%% Initializing parameters
TH_input=readmatrix(filename);
n=length(D);
d = zeros(n);
for ii=1:n
    for jj=1:n
        d(ii,jj)= (D(jj)-D(ii));
    end
end

for cc=1:1
for qq= 1:length(interval)-1
    start_time= interval(qq); %Start time of the segment
    end_time= interval(qq+1); %End time of the segment
    dt=diff(TH_input(1:2,1));%Time Step
    %
    index1= find(abs(TH_input(:,1)-start_time)<=1e-10);
    index2= find(abs(TH_input(:,1)-end_time)<=1e-10);

    TH= TH_input(index1:index2,:);
    TH(:,2)= tukeywin(length(TH),0.1).*TH(:,2);
    
    N=length(TH);
    f=(0:(N-1))/(N*dt); %frequency matrix
    df=diff(f(1:2));%frequency increment

    %% Fourier coefficients for input TH
    A0 = sum(TH(1:N,2))/N; %corresponding to wp=0

    Ap = zeros(N/2-1,1);%corresponding to pth frequency
    Bp = zeros(N/2-1,1);

    for p=1:N/2-1
        Ap(p,1)=2/N*sum(TH(1:N,2).*cos(2*pi*p/N/dt*TH(1:N,1)));
        Bp(p,1)=2/N*sum(TH(1:N,2).*sin(2*pi*p/N/dt*TH(1:N,1)));
    end
    An2=0; %corresponding to Nyquist frequency
    for ii=1:N
        An2 = An2 + (TH(ii,2)/N)*(-1)^(ii);
    end

    I=zeros(N/2+1,1); %Estimator for Auto-PSD at station-1 from equation 10
    I(1)=N*dt/4/pi*A0^2;
    I(2:N/2)=N*dt/4/pi*(Ap(1:N/2-1).^2+Bp(1:N/2-1).^2);
    I(N/2+1)=N*dt/4/pi*An2^2;

    %Matrices to store fourier coefficients of generated TH
    APK=zeros(N/2-1,(n-1));
    BPK=zeros(N/2-1,(n-1));
    AN2=zeros(1,(n-1));

    for ff=2:N/2+1
        %% Frequency Response Function

        %from equation 8
        H=zeros(n,2);
        for ii=1:n %for +ve frequency
            H(ii,1)=(wk(ii)^2+2*1i*zk(ii)*wk(ii)*2*pi*f(ff))/(wk(ii)^2-(2*pi*f(ff))^2+2*1i*zk(ii)*wk(ii)*2*pi*f(ff));
        end
        for ii=1:n %for -ve frequency
            H(ii,2)=(wk(ii)^2+2*1i*zk(ii)*wk(ii)*2*pi*(-f(ff)))/(wk(ii)^2-(2*pi*(-f(ff)))^2+2*1i*zk(ii)*wk(ii)*2*pi*(-f(ff)));
        end
        %% Coherancy Function
        gamma= zeros(n);
        for ii=1:n
            for jj=1:n
                theta_sr= atan2(imag(H(ii,1)*H(jj,2)),real(H(ii,1)*H(jj,2))); %phase angle for site response
                %considering incoherance and site effects from coherancy model of Kiureghian
                gamma(ii,jj)=exp(-(alpha*d(ii,jj)*2*pi*f(ff)/Vs)^2)*exp(1i*(theta_sr));
            end
        end

        %% Auto and Cross PSD
        G=zeros(n);
        G(1,1)=I(ff); %doubt

        for ii=2:n
            G(ii,ii)=G(1,1)*abs(H(ii,1))^2/abs(H(1,1))^2; %Auto PSD at unknown station from equation 12
        end

        for ii=1:n
            for jj=1:n
                if ii==jj
                    G(ii,ii)=G(ii,ii);
                elseif ii~=jj
                    G(ii,jj)=gamma(ii,jj)*sqrt(G(ii,ii)*G(jj,jj)); %Cross PSDs from equation 3
                end
            end
        end

        %% Covariance Matrix
        %From the formulation given in equation 2 from Konalki & Kiureghian 2012
        %details of the formulation attached in the excelfile named
        %COVARIANCE_MATRIX

        S=zeros(2*n);
        for ii=1:2*n
            for jj=ii:2*n
                if ii==jj %Diagonal elements
                    S(ii,jj)=G(ceil(jj/2),ceil(jj/2))*2*pi*df;

                elseif mod(ii,2)==1 && jj==ii+1
                    S(ii,jj)=0;

                elseif mod(ii,2)==1 && mod(jj,2)==1 %odd odd
                    S(ii,jj)=real(G(ceil(ii/2),ceil(jj/2)))*2*pi*df;

                elseif mod(ii,2)==0 && mod(jj,2)==0 %even even
                    S(ii,jj)=real(G(ceil(ii/2),ceil(jj/2)))*2*pi*df;

                elseif mod(ii,2)==1 && mod(jj,2)==0 %odd even
                    S(ii,jj)=-imag(G(ceil(ii/2),ceil(jj/2)))*2*pi*df;

                else
                    S(ii,jj)=imag(G(ceil(ii/2),ceil(jj/2)))*2*pi*df;
                end
                S(jj,ii)=S(ii,jj); %symmetric matrix
            end
        end
        %% Conditioned Mean and Covariance Matrix
        if ff<=N/2
            xp1= [Ap(ff-1);Bp(ff-1)] ;%Fourier coefficents vector for known site
        else
            xp1=[An2;0];%doubt
        end
        %generation of sub matrices
        Spp11= S(1:2,1:2);
        Spp12= S(1:2,3:2*n);
        Spp21= S(3:2*n,1:2);
        Spp22= S(3:2*n,3:2*n);

        M=Spp21/(Spp11)*xp1; %mean matrix equation 13
        s= Spp22-Spp21/Spp11*Spp12; %covariance matrix equation 14

        zp=zeros(2*(n-1),1); %generation of uncorrelated standard normal variables
        for ii=1:2*(n-1)
            zp(ii,1)=randn(1,1);
        end
        xp2=M+chol(s)'*zp; %formulation for fourier coefficients for unknown sites

        %% Seperation of Apk and Bpk terms and AN2k
        if ff<=N/2
            APK(ff-1,:)=reshape(xp2(1:2:end,1),n-1,[])';
            BPK(ff-1,:)=reshape(xp2(2:2:end,1),n-1,[])';
        else
            AN2=reshape(xp2(1:2:end,1),n-1,[])'; %doubt
        end
    end
    %% Applying extension to the time window
    len_taper=(round((index2-index1)*0.03)); % 3% overlap considered

    t1=index1-len_taper;
    t2=index2+len_taper;

    if t1<=0
        a= [TH;dt*(index2:t2-1)',zeros(abs(t2-index2),1)];
    elseif t2>= TH_input(end,1)/dt
        a= [dt*(t1-1:index1-2)',zeros(abs(index1-t1),1);TH];
    else
        a= [dt*(t1-1:index1-2)',zeros(abs(index1-t1),1);TH;dt*(index2:t2-1)',zeros(abs(t2-index2),1)];
    end

    A= [a,zeros(length(a),n-1)];
    %% Generation of Time Histories
    %using equation 1 summing up over all frequencies
    for kk=3:n+1
        for ii=1:length(A)
            for pp=1:N/2-1
                A(ii,kk)=APK(pp,kk-2)*cos(2*pi*f(pp+1)*A(ii,1))+BPK(pp,kk-2)*sin(2*pi*f(pp+1)*A(ii,1))+A(ii,kk);
            end
        end
    end

    A(:,3:n)=A(:,3:n)+A0; %summing up zeroth frequency

    for ii=1:length(A)
        A(ii,3:n+1)=A(ii,3:n+1)+AN2*(-1)^(ii-1); %summing up N/2th frequency
    end

    if qq>1

        start_index=find(abs(A_previous(:,1)-(A(1,1)))<1e-10); %Index of start time of second segment in first segment
        end_index=length(A_previous); %Index of end of first segment
        len_taper= end_index-start_index+1; %Length of overlap of two segments
        Y=zeros(end_index-start_index+1,min(size(A_previous))); %Variable to store the combined time histories
        Y(:,1)=A_previous(start_index,1):dt:A_previous(end,1);


        for ee=2:min(size(A_previous))
            hh=1;
            for ii=1:length(Y)
                Y(ii,ee)=cos(0.5*pi*((start_index-1+ii)-start_index)/(end_index-start_index))*A_previous(start_index-1+ii,ee)+(1-cos(0.5*pi*((start_index-1+ii)-start_index)/(end_index-start_index)))*A(hh,ee); %Weighted cosine function
                hh=hh+1;
            end
        end

        % Combine the weighted overlapping intervals of TH1 and TH2 to create a single time-history.
        A= [A_previous(1:start_index-1,:); Y; A(len_taper+1:end,:)];
    else
    end
    A_previous=A;
end

%% Post Processing of simulated time-histories

%Initalizing original time-history
A(:,2)=TH_input(:,2);

%Applying cosine tapper to set initial accleration to zero
initial_taper= round(length(A)*0.02); %considering 2% inital taper 

for ii=1:initial_taper
A(ii,3:end)=(1-cos(pi*(ii)/(2*initial_taper)))*A(ii,3:end);
end

%Applying high pass filter to enforce zeros residual velocity and
%displacemnt
fc = 0.1; % corner frequency
fs = 1/(A(2,1)-A(1,1)); % sampling frequency
order = 4; % filter order
Wn = fc / (fs / 2); % Calculate normalized frequency
[b, a] = butter(order, Wn, 'high');% Design Butterworth filter
for ii=3:min(size(A))
A(:,ii) = filtfilt(b, a, A(:,ii));% Apply filter to acceleration data
end

%Subtracting mean- Base line correction 
A(:,3:end) = detrend(A(:,3:end));
velocity = cumtrapz(dt,A(:,2:end)); %calculating velocity time history
disp     = cumtrapz(dt,velocity); %calculating displacement time history

%% Generating the plot and sheet


R=zeros(1201,n);

for ii=2:n+1
    R(:,ii-1)=responsespectrum(A(1:end,ii),5);
end
col_letter = char('A' + mod(cc-1, 20));
for ii = 1:n
    eval(sprintf('RS_%d = R(:,%d);', ii,ii));
    
    output_filename = sprintf('RS_%s_%d.xlsx',filename,ii); 
    writematrix(eval(sprintf('RS_%d', ii)), output_filename, 'Range', sprintf('%s1:%s1200', col_letter, col_letter));
end
for ii = 1:n
    eval(sprintf('acc_%d = A(:,%d+1);', ii,ii));
    
    output_filename = sprintf('acc__%s_%d.xlsx',filename,ii); 
    writematrix(eval(sprintf('acc_%d', ii)), output_filename, 'Range', sprintf('%s1:%s%d', col_letter,col_letter, length(A)));
end
for ii = 1:n
    eval(sprintf('disp_%d = disp(:,%d);', ii,ii));
    
    output_filename = sprintf('disp_%s_%d.xlsx',filename,ii); 
    writematrix(eval(sprintf('disp_%d', ii)), output_filename, 'Range', sprintf('%s1:%s%d', col_letter, col_letter,length(A)));
end
for ii = 1:min(size(A))-1
    eval(sprintf('a_%d = [A(:,1) A(:,%d+1)];', ii,ii));
    
    output_filename = sprintf('%d_%s_%d.txt',cc,filename,ii); 
    writematrix(eval(sprintf('a_%d', ii)),output_filename, 'delimiter', ' ');
end
for ii=2:n
    subplot((n-1),2,ii-1)
    plot(TH_input(1:length(A),1),TH_input(1:length(A),2),A(:,1),A(:,ii+1));
    xlabel('Time (s)')
    ylabel('Acceleration (g)')
end
for ii=2:n
    subplot((n-1),2,n+ii-2)
    plot(TH_input(1:length(A),1),disp(:,1),TH_input(1:length(A),1),disp(:,ii));
    xlabel('Time (s)')
    ylabel('displacment (m)')
end

end

function [A] = responsespectrum(accel, ee)
close all
% ee - damping in % - 5 is recommended
% y - gamma in newmark's method - 0.5 is recommended
% b - beta in newmark's method - 0.25 is recommended
% td - time till which you want graph to be plotted
Tn = 6; % time period till which u want response spectrum
y = 0.5;
b = 0.25;
uo = 0;
vo = 0;
m = 1;
z = ee/100;
dt = 0.005;
na = length(accel);
nl = 2*na;
T = (0.005:0.005:Tn)';
accel = [accel;zeros(nl-na,1)];
p = -m*accel;
A = zeros(length(T),1);% acclelration response spectrum - total accelaration
V = zeros(length(T),1);% velocity response spectrum - relative velocity
D = zeros(length(T),1);% displacement response spectrum - relative displacement
for j = 1:length(T)
    
    fn = 1/T(j);
    wn = 2*pi*fn;
    k = m*wn^2;
    c = 2*m*wn*z;
    
    u = zeros(nl,1);
    v = zeros(nl,1);
    ac = zeros(nl,1);
    
    u(1) = uo;
    v(1) = vo;
    ac(1) = (p(1)-c*vo-k*uo)/m;
    
    kf = k + y*c/(b*dt) + m/(b*dt^2);
    a = m/(b*dt) + y*c/b;
    b2 = m/(2*b) + dt*(y/(2*b) - 1)*c;
    
    for i = 1:nl-1
        p1=p(i);
        p2=p(i+1);
        dpf = (p2 - p1) + a*v(i) + b2*ac(i);
        du = dpf/kf;
        dv = y/(b*dt)*du - (y/b)*v(i) + dt*(1 - y/(2*b))*ac(i);
        da = du/(b*dt^2) - v(i)/(b*dt) - ac(i)/(2*b);
        u(i+1) = u(i) + du;
        v(i+1) = v(i) + dv;
        ac(i+1) = ac(i) + da;
    end
    
    asd = ac + accel;
    A(j) = max(abs(asd));
    V(j) = max(abs(v));
    D(j) = max(abs(u));
end
A = [max(abs(accel(:)));A];
end