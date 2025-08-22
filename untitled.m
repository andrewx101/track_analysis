% Parameters
M      = 20;   % # of tracers
L      = 2000;    % box size

% given M, L
rho = M / L^2;             
E_R = 1/(2*sqrt(rho));     
xi_min = 0.3*E_R;          
xi_max = 3*E_R;   

fprintf('Mean NN distance ≃ %.1f px.\n', E_R);
fprintf('Recommend ξ in [%.1f, %.1f] px for nonzero χ_cross.\n', xi_min, xi_max);

% pick one ξ in that interval
xi = E_R;    % e.g. set ξ = mean NN ≃ 224 px

xi=xi/10;

% determine Nx for p grid‐points/ξ
p      = 8;  
Nx_raw = ceil(p*L/xi);
Nx     = 2^nextpow2(Nx_raw);

fprintf('Using ξ = %.1f px → Nx = %d (Δx = %.2f px)\n', xi, Nx, L/Nx);

D0     = 10;     % mean diffusivity
sigmaD = 10;   % std of diffusivity
dim    = 2;     % 1/2/3
dt     = 0.01;  % time step
Nsteps = 5000;  % total number of steps
nLag   = 100;   % lag in steps for overlap: Delta_t = nLag*dt


%--- Build correlated Gaussian field G(r) in arbitrary dim --------------%
h = L/Nx;
m = ceil(3*xi/h);
coords = (-m:m)*h;    % vector of length 2m+1

switch dim
    case 1
        % 1D kernel
        K = exp(-coords.^2/(2*xi^2));
        K = K / sum(K);
        W = randn(Nx,1);
        G = conv(W, K, 'same');

    case 2
        % 2D kernel
        [X,Y] = meshgrid(coords, coords);
        K = exp(-(X.^2 + Y.^2)/(2*xi^2));
        K = K / sum(K(:));
        W = randn(Nx,Nx);
        G = conv2(W, K, 'same');

    case 3
        % 3D kernel
        [X,Y,Z] = ndgrid(coords, coords, coords);
        K = exp(-(X.^2 + Y.^2 + Z.^2)/(2*xi^2));
        K = K / sum(K(:));
        W = randn(Nx,Nx,Nx);
        G = convn(W, K, 'same');

    otherwise
        error('dim must be 1, 2 or 3');
end

% Standardize G to zero mean, unit variance
G = (G - mean(G(:))) / std(G(:));

%--- Exponentiate to get a strictly positive, log‐normal D_field ------%
sigma2 = log(1 + (sigmaD/D0)^2);
sigmaLN = sqrt(sigma2);
muLN    = log(D0) - 0.5*sigma2;
Df = exp(muLN + sigmaLN * G);   % Df is now >0 everywhere

%--- Build grid coordinates for interpolation ------------------------%
switch dim
    case 1
        xgrid = (0:Nx-1)*h;
        Fi=griddedInterpolant({xgrid},Df,'linear','nearest');

    case 2
        xgrid = (0:Nx-1)*h;
        ygrid = xgrid;
        Fi = griddedInterpolant({xgrid,ygrid}, Df, 'linear', 'nearest');

    case 3
        xgrid = (0:Nx-1)*h;
        ygrid = xgrid;
        zgrid = xgrid;
        Fi=griddedInterpolant({xgrid,ygrid,zgrid},Df,'linear','nearest');

end
% Fi(xi, yi) now returns Df at all (xi,yi), defaulting to 'nearest' if out‐of‐bounds.

%--- Simulate M tracers ------------------------------------------------
% initialize
Xraw = zeros(M,dim,Nsteps+1);
X    = zeros(M,dim,Nsteps+1);
Xraw(:,:,1) = L*rand(M,dim);
X(:,:,1)    = Xraw(:,:,1);

for t=1:Nsteps
    % 1) gather positions
    Xi_wrapped = X(:,:,t);          % M×dim
    % 2) interpolate vectorized
    switch dim
        case 1
            Di=Fi(Xi_wrapped(:,1));
        case 2
            Di = Fi(Xi_wrapped(:,1), Xi_wrapped(:,2));
        case 3
            Di=Fi(Xi_wrapped(:,1),Xi_wrapped(:,2),Xi_wrapped(:,3));
    end
    % 3) random step for all
    dW = randn(M,dim)*sqrt(dt);
    dXraw      = sqrt(2*Di).*dW;          % true small step
    Xraw(:,:,t+1) = Xraw(:,:,t) + dXraw;   % unwrapped pos
    X(:,:,t+1)    = mod(Xraw(:,:,t+1), L); % for interpolation next stepend
end

close all

tr=[];
for i=1:M
    
    Xpos=squeeze(Xraw(i,1,:));
    Ypos=squeeze(Xraw(i,2,:));
    step_id=1:length(Xpos);
    tr=[tr;[Xpos(:),Ypos(:),step_id(:),ones(length(step_id),1).*i]];

    
    plot(Xpos,Ypos,'.-')

    


    hold on;
end
hold off

