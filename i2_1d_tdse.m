close all
clear all 
%constants & conversions
h_bar=1; i1=sqrt(-1);
AU_TO_EV=27.2113961;
NM_TO_AU=45.56335;
AU_TO_CM=2.194746e5;
AU_TO_ANG=0.529177249;
AU_TO_FS=0.0241889;
AMU_TO_AU=1.8228885e3;
WCM2_TO_AU=3.5094452*1e16;
conv_r= AU_TO_ANG;
conv_en=AU_TO_CM;
conv_t= AU_TO_FS;

% Turn on and off plotting of the wavepacket wavefunction
plot_flag=1;

% set up grids (in Bohr)
N=2^12; % grid points
Rmin = 2;  Rmax = 66; % box size
R = linspace(Rmin,Rmax,N);
dR = R(2)-R(1);

% mass of I2
mass = (126.904473/2)*AMU_TO_AU;

% set up k grid
dk = 2*pi/(N*dR);
k =[(0:N/2-1),(-N/2:-1)]*dk;

% set up kinetic energy and k-space propagator
kin= (0.5*(k.^2)/mass).';

% X (ground) state of I2
De_X= 12440/AU_TO_CM;
be_X= 1.875*AU_TO_ANG;
Re_X= 2.656/AU_TO_ANG;
V1= De_X*(1-exp(-be_X*(R-Re_X))).^2;

% B (excited) state of I2
De_B=  5169/AU_TO_CM;
be_B= 1.696*AU_TO_ANG;
Re_B= 3.025/AU_TO_ANG;
Te_B= 15771/AU_TO_CM;
V2= De_B*(1-exp(-be_B*(R-Re_B))).^2 + Te_B;

% C (or B'' excited) state of I2
De_C=.0543;
alpha_C= -3.6060;
be_C=  1.6041;
Re_C=3.5161;
Te_C= 522.5292/AU_TO_CM;
V3= De_C*(1-alpha_C*exp(-be_C*(R-Re_C))).^2 + Te_C;

%% get initial psi, use eigenstate of ground state via Fourier Grid Method
kk=pi/dR;
% kinetic energy
[nx, ny]=meshgrid(1:N);
T=-2*kk^2*(-1).^(nx-ny)./(pi*(nx-ny)).^2;
T(1:N+1:end)=0; % zeros inf on diagonal
T=-(0.5*h_bar^2/mass)*(T-eye(N)*kk^2/3);

[evec1,evl1]=eig(T+diag(V1)); % eigen-states of X-state
norm1= 1/sqrt(dot(evec1(:,1), evec1(:,1)));
psi_GS= norm1*evec1(:,1);
E_GS=evl1(1,1);


update_freq= 1; % update_freq*t_sample is the interval at which plots updated
%%
% laser pulse parameters
I0     =  1e11;            % peak power in W/cm^2
FWHM   =  40.0/AU_TO_FS;   % pulse FWHM duration
lambda =  520;             % pulse wavelength in nm
w0     =  NM_TO_AU/lambda; % pulse frequency
t0     =  3*FWHM;          % pulse envelope center

% Pulse intensities in atomic units
E0=sqrt(I0/WCM2_TO_AU);    % the peak field amplitude in a.u.

% Transition dipole in atomic units
mu12 = 0.181;

% time grid propagation parameters
dt       = 1/AU_TO_FS;
t_sample = 15/AU_TO_FS;
t_end    = t0 + 1900/AU_TO_FS;
n_steps  = round(t_sample/dt);
n_sample = round(t_end/t_sample);
t_end    = n_sample*t_sample;

t        = 0:dt:t_end; % time grid

% The laser pulse
E        = E0*exp(-2*log(2)/FWHM^2*(t-t0).^2).*exp(-1i*w0*(t-t0));

%% propagate !

% set up propagators
u_kin= exp(-1i*dt*kin);
u_kin_half= exp(-1i*0.5*dt*kin);
u_pot1= exp(-1i*dt*(V1.'));
u_pot2= exp(-1i*dt*(V2.'));
u_pot3= exp(-1i*dt*(V3.'));

psi_1=    psi_GS;          % The GS unperturbed wavefunction
psi_2=    zeros(N, 1);     % First order wavefunction
psi_3=    zeros(N, 1);     % First order wavefunction

% Plot potential curves and initial wavefunction
if (plot_flag)
    h1=figure(1);
    set(h1,'position',[50, 50,500, 650])%
    subplot(1,1,1);
    clf;
    pot_plot=plot( R*conv_r, [V1*conv_en; V2*conv_en ; V3*conv_en]);
    axis manual;
    axis([conv_r*min(R) conv_r*max(R) 0 1.25*conv_en*V2(N)]);
    hold on;
    wave_scale= 0.2*De_X*conv_en/(max(psi_GS)-min(psi_GS));
    
    psi_plot=plot( R*conv_r, 0*wave_scale*psi_GS+ E_GS*conv_en, ...
        R*conv_r, 0*psi_GS+ (E_GS+w0)*conv_en,...
        R*conv_r, 0*psi_GS+ (E_GS+w0)*conv_en);
    xlim([2 8]);
end

in=1; %
tot_steps=0;
t=0;
n_sample=0;
max_steps= round(t_end/dt);
n=0;
while (tot_steps<max_steps)
    % propagate wavefunctions 1/2 step forward
    U= u_kin_half;
    psi_2=    ifft(U.*fft(psi_2));
    psi_3=    ifft(U.*fft(psi_3));
    
    % Propagate until next sample time
    for step=1:n_steps
        % potential step
        psi_1=    psi_1*exp(-1i*dt*E_GS);
        psi_2=    u_pot2.*psi_2;
        psi_3=    u_pot3.*psi_3;
        %psi_1=    u_pot1.*psi_1;
        % kinetic step
        if (step==n_steps)
            % apply only 1/2 kinetic step when wavefunction is sampled
            U= u_kin_half;
        else
            U= u_kin;
        end
        psi_3=    ifft(U.*fft(psi_3));
        psi_2=    ifft(U.*fft(psi_2));
        %psi_1=    ifft(U.*fft(psi_1));
        % Add contribution due to field interactions
        W1= E(tot_steps+1)*mu12;
        psi_2=    psi_2    + 1i*dt*(W1.*psi_1);
        psi_3=    psi_3    + 1i*dt*(W1.*psi_1);
        %psi_1=    psi_1    + 1i*dt*(W1.*(psi_1));
        t=t+dt;
        tot_steps=tot_steps+1;
        
        
    end
    n_sample=n_sample+1;
    pop= dot(psi_2,psi_2);
    pop_2(n_sample)= dot(psi_2,psi_2);
    pop_3(n_sample)= dot(psi_3,psi_3);
    
    n=n+1;
    PSI3(:,n)=abs(psi_3).^2;
    
    PSI2(:,n)=abs(psi_2).^2;
    PSI1(:,n)=abs(psi_1).^2;
    time(n)=t.*conv_t;
    if (pop)
        r_1(n_sample)= real(dot(psi_2, (R.').*psi_2)/pop_2(n_sample));
        r_1_sqr(n_sample)= real(dot(psi_2, ((R.^2).').*psi_2)/pop_2(n_sample));
        
        r_2(n_sample)= real(dot(psi_3, (R.').*psi_3)/pop_3(n_sample));
        r_2_sqr(n_sample)= real(dot(psi_3, ((R.^2).').*psi_3)/pop_3(n_sample));
        
    end
    if (plot_flag & mod(tot_steps/n_steps, update_freq)==0)
        set(psi_plot(1), 'YData', 1*wave_scale*abs(psi_1)+ ...
            (E_GS)*conv_en);
        set(psi_plot(2), 'YData', 8*wave_scale*abs(psi_2)+ ...
            (E_GS+w0)*conv_en);
        set(psi_plot(3), 'YData', 8*wave_scale*abs(psi_3)+ ...
            (E_GS+w0)*conv_en);
        
        s= sprintf('\\psi_2: Time= %.0f fs', conv_t*t);
        title(s,'fontsize',6);
        drawnow;
        
        f(in) = getframe(h1);
        in = in+1;
        
    end
end

%%
%Plot pulses and populations
figure(2);
clf;
subplot(2,2,1);
plot((0:dt:t_end)*conv_t, abs(E).^2);
ylabel('Field intensity');
subplot(2,2,2);
plot((t_sample:t_sample:t_end)*conv_t, pop_2,(t_sample:t_sample:t_end)*conv_t, pop_3);
ylabel('Excited state population');
subplot(2,2,3);
plot((t_sample:t_sample:t_end)*conv_t, r_1*conv_r,(t_sample:t_sample:t_end)*conv_t, r_2*conv_r);
ylabel('<bond length>');
subplot(2,2,4);
plot((t_sample:t_sample:t_end)*conv_t, sqrt(r_1_sqr - r_1.^2)*conv_r, ...
    (t_sample:t_sample:t_end)*conv_t, sqrt(r_2_sqr - r_2.^2)*conv_r);
ylabel('<dbond>');
xlabel('time (fs)');



figure(100)
imagesc((bsxfun(@times,PSI1,1-pop_2-pop_3)+bsxfun(@times,PSI2,1./pop_2))+bsxfun(@times,PSI3,1./pop_3));colorbar
%caxis([-10 -1])