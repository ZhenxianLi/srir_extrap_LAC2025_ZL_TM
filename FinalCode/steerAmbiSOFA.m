function sofaOutFile = steerAmbiSOFA(sofaInFile, elDeg, azDeg, beta, sigma, sofaOutFile)
% steerAmbiSOFA  ── Beam-steer (focus) a Higher-Order Ambisonic SRIR/BRIR
%                 stored in SOFA format and write the result as a new SOFA.
%
%   sofaOutFile = steerAmbiSOFA(sofaInFile, elDeg, azDeg, beta, sigma)
%   sofaOutFile = steerAmbiSOFA(sofaInFile, elDeg, azDeg, beta, sigma, outName)
%
%   INPUT ARGUMENTS
%       sofaInFile   : Path to the source SOFA file (ACN/N3D channel order).
%                      The file can contain either SRIRs or arbitrary HOA
%                      block signals in Data.IR  (size: M×C×Ns).
%       elDeg        : Elevation angle of the target look-direction
%                      in degrees.   +90 = zenith, 0 = horizontal plane,
%                      -90 = nadir.   (ITU convention.)
%       azDeg        : Azimuth angle of the target look-direction
%                      in degrees.   0 = front, +90 = left, increasing CCW
%                      around the listener.   (ITU convention.)
%       beta         : Linear cross-fade factor, 0…1.
%                      0   → pure original sound field,
%                      1   → only the beam-steered (focused) version.
%       sigma        : Beam-width parameter (>0).  Rough guideline:
%                      σ ≈ 1   → narrow beam (highly directional),
%                      σ = 2-3 → wider, softer beam.
%       sofaOutFile  : (Optional) Desired output filename.  If omitted, a
%                      name in the form *_steer_<el>_<az>_beta<β>.sofa is
%                      generated alongside the input file.
%
%   OUTPUT
%       sofaOutFile  : The actual filename that has been written to disk.
%
%   ALGORITHM SUMMARY
%       1) Load the SOFA file and extract Data.IR (M measurements ×
%          C HOA channels × Ns samples).
%       2) Build a real ACN/N3D spherical-harmonic steering vector Y that
%          points to (az,el).   Apply a per-order Gaussian window
%          w_n = exp(-n/σ) to control beam width (σ).
%       3) For every measurement (row m) perform:
%          • Decode to a plane-wave signal: s(t) = X · Yᵀ
%          • Encode the plane wave back to HOA: X_steered = s(t) · Y
%          • Cross-fade: X_out = (1–β)·X  +  β·X_steered
%       4) Store X_out back into S.Data.IR and save as a new SOFA file.
%
%   REQUIREMENTS
%       • SOFA API for MATLAB/Octave (run SOFAstart once before calling).
%       • Input HOA channels must follow ACN order and N3D normalization.
%
%   2025-08-01   Example code by Zhenxian Li, commented & refined by ChatGPT-o3
% --------------------------------------------------------------------------
% --- 0)  Validate/complete optional arguments ------------------------------
if nargin < 6 || isempty(sofaOutFile)
    [pth,nm,~] = fileparts(sofaInFile);
    sofaOutFile = fullfile( ...
        pth, ...
        sprintf('%s_steer_%d_%d_beta%.2f.sofa', nm, round(elDeg), round(azDeg), beta));
end
if nargin < 5 || isempty(sigma), sigma = 1;   end
if nargin < 4 || isempty(beta),  beta  = 0.5; end

% --- 1) Load the SOFA file -------------------------------------------------
S  = SOFAload(sofaInFile);        % Make sure SOFAstart has been executed
IR = double(S.Data.IR);           % Cast to double for numerical accuracy
[M, C, Ns] = size(IR);            % M = measurement count,
                                  % C = HOA channels, Ns = samples/IR length
N = round(sqrt(C) - 1);           % HOA order (because C = (N+1)^2 in 3-D)

% --- 2) Create steering vector Y with optional Gaussian per-order taper ----
az = deg2rad(azDeg);
el = deg2rad(elDeg);
Y  = getRealSH_N3D(N, az, el);    % 1×C steering vector (unit-energy)

% Build Gaussian weights: w_n = exp(–n/σ) then expand to channel vector
w_per_order = exp( -(0:N).' / sigma );   % column vector, (N+1)×1
w_vec      = expandWeights_ACN(w_per_order);  % 1×C, repeated per order
Yw         = Y .* w_vec;                 % apply taper
Yw         = Yw / norm(Yw);              % re-normalize to unit energy

% --- 3) For every measurement, decode – re-encode – cross-fade ------------
Yw_col = Yw.';                            % C×1, for matrix multiply
for m = 1:M
    X  = squeeze(IR(m,:,:)).';           % Ns×C   (transpose to sample rows)
    
    s  = X * Yw_col;                     % Decode: Ns×1 plane-wave signal
    Xsteered = s * Yw;                   % Encode: Ns×C steered HOA signal
    
    % Simple linear cross-fade.
    % (Replace by equal-power fade if you prefer: g1=cos(pi/2*beta); etc.)
    g1 = 1 - beta;
    g2 = beta;
    
    Xout = g1 * X + g2 * Xsteered;       % Ns×C
    
    % Write back.  IR is M×C×Ns, but Xout' is C×Ns, so permute dimensions:
    IR(m,:,:) = permute(Xout.', [1 3 2]);   % (C×Ns) → (1×C×Ns)
end
S.Data.IR = single(IR);                   % Cast back to single if desired

% --- 4) Save the result as a new SOFA file --------------------------------
SOFAsave(sofaOutFile, S);
fprintf('[steerAmbiSOFA]  Finished writing →  %s\n', sofaOutFile);
end
% ==========================================================================

% ---------- Helper function: real ACN/N3D spherical harmonics -------------
function Y = getRealSH_N3D(N, az, el)
% getRealSH_N3D  Compute real, N3D-normalized SH coefficients (ACN order)
%
%   Y = getRealSH_N3D(N, az, el)
%       N    : maximum SH order
%       az   : azimuth  (rad)
%       el   : elevation (rad)
%       Y    : 1×(N+1)^2 vector of real SHs  (unit-energy)
%
theta = pi/2 - el;                % colatitude
nCh   = (N+1)^2;
Y     = zeros(1, nCh);            % pre-allocate

for n = 0:N
    % Associated Legendre polynomials P_n^m(cos θ), m = 0…n
    P = legendre(n, cos(theta), 'sch');   % 'sch' gives Schmidt-semi-norm.
    
    for m = -n:n
        idx = n*(n+1) + m + 1;    % ACN channel index (0-based → 1-based)
        Pnm = P(abs(m)+1);        % Retrieve correct order |m|
        
        % N3D normalization factor K_n^m
        K = sqrt((2*n+1)/(4*pi) * factorial(n-abs(m))/factorial(n+abs(m)));
        
        if m > 0
            Y(idx) = sqrt(2) * K * Pnm * cos( m * az );
        elseif m < 0
            Y(idx) = sqrt(2) * K * Pnm * sin( abs(m) * az );
        else % m == 0
            Y(idx) =         K * Pnm;
        end
    end
end
% Y is already N3D-normalized; ensure unit-energy for safety
Y = Y / norm(Y);
end
% --------------------------------------------------------------------------

% ---------- Helper function: expand per-order weights to ACN vector -------
function wvec = expandWeights_ACN(w_order)
% expandWeights_ACN  Convert an (N+1)×1 vector of per-order weights
%                    into a 1×(N+1)^2 vector aligned to ACN channel order.
%
N = numel(w_order) - 1;           % highest order
wvec = zeros(1, (N+1)^2);
idx  = 1;
for n = 0:N
    numM = 2*n + 1;               % number of SHs for this order
    wvec(idx : idx+numM-1) = w_order(n+1);
    idx = idx + numM;
end
end
