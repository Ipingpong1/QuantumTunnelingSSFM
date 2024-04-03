N = 512; % Number of Fourier modes
dt = .0005; % Size of time step
tfinal = 100; % Final time
M = round(tfinal/dt); % Total number of time steps
J = 100; % Steps between output
L = 20; % Space period
q = 1;
v = 1;
h = L/N; % Space step
n = (-N/2:1:N/2-1); % Indices
ts = (0:1:M); % Indices

x = n*h; % Grid points
u0 = exp( -((x+6).^2) ) .* exp(11i.*x);
u0 = normalize(u0);
u = u0; % Initial Condition
k = 2*n.*pi/L; % Wavenumbers.
V = 1*(x.^2)/50*0;
val = ts.*0 + 1;


t = 0;
for m = 1:1:M % Start time loop
    
    u = exp(dt*1j*V).*u; % Solve non-constant part of LSE
    c = fftshift(fft(u)); % Take Fourier transform
    c = exp(-dt*1j*k.^2).*c; % Advance in Fourier Space
    u = ifft(fftshift(c)); % Return to Physical Space
    
    V(320:340)=1500;
    V(1) = 10000;
    V(N) = 10000;
    u(1) = 0;
    u(N) = 0;

    
    hold on
    plot(ts, val)
    disp(trapz(x, abs(u).^2))

    hold off
    pause(.001)
    clf()
 
end
