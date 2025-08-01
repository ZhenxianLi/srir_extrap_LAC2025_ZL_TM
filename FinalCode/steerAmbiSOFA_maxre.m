%% Narrow, very focused beam, broadband
%out1 = steerAmbiSOFA_maxre(inFile,  10, -30, 'Beta',0.8);

%% Wide beam (σ = 2) below 500 Hz, max-rE above
%out2 = steerAmbiSOFA_maxre(inFile, -15,  90, ...
              'Beta',0.5, 'FreqMode',true, 'Fsplit',500, ...
              'Sigma',2, 'Weighting','gaussian');

%% GPU-accelerated processing with equal-power fade
%out3 = steerAmbiSOFA_maxre(inFile,  0,   0, ...
              'Beta',0.6, 'FadeLaw','equalpower', 'UseGPU',true);


function sofaOut = steerAmbiSOFA_maxre(sofaIn, elDeg, azDeg, varargin)
% steerAmbiSOFA_maxre  |  Advanced HOA beam-steering with max-rE weighting
%
%   sofaOut = steerAmbiSOFA_maxre(sofaIn, elDeg, azDeg, ...
%                'Beta',0.6, 'Weighting','maxre', 'FadeLaw','equalpower',
%                'FreqMode',true, 'Fsplit',600, 'Sigma',1.5,
%                'UseGPU',false, 'OutFile','mySteered.sofa')
%
% INPUTS
%   sofaIn     – Path to ACN/N3D SOFA file containing HOA IRs or blocks
%   elDeg      – Elevation   (deg, +90=zenith)
%   azDeg      – Azimuth     (deg, 0=front, +90=left, CCW)
%
% NAME–VALUE PAIRS  (defaults shown)
%   'Beta'        , 0.5        – Cross-fade ratio (0–1)
%   'Weighting'   , 'maxre'    – {'maxre'|'gaussian'|'custom'|'none'}
%   'Sigma'       , 1          – Beam width for Gaussian (if selected)
%   'CustomWeights', []        – Vector w_n (overrides both above)
%   'FadeLaw'     , 'linear'   – {'linear'|'equalpower'}
%   'FreqMode'    , false      – true = STFT per-band processing
%   'Fsplit'      , 500        – Crossover [Hz] if FreqMode==true
%   'STFTlen'     , 1024       – Window length (samples)
%   'UseGPU'      , false      – Send heavy matrices to GPU
%   'OutFile'     , ''         – Filename of output SOFA
%
% OUTPUT
%   sofaOut   – Path to the written SOFA file
%
% ----------------------------------------------------------------------
% 2025-08-01  |  Z. Li & ChatGPT-o3  |  Prototype for research use
% ----------------------------------------------------------------------

%% ---------- 0)  Parse options -----------------------------------------
P = inputParser;
P.addParameter('Beta',              0.5,        @(x) isnumeric(x)&&isscalar(x));
P.addParameter('Weighting',         'maxre',    @(s) ischar(s)||isstring(s));
P.addParameter('Sigma',             1,          @(x) isnumeric(x)&&isscalar(x));
P.addParameter('CustomWeights',     [],         @(x) isnumeric(x));
P.addParameter('FadeLaw',           'linear',   @(s) any(strcmpi(s,{'linear','equalpower'})));
P.addParameter('FreqMode',          false,      @islogical);
P.addParameter('Fsplit',            500,        @(x) isnumeric(x)&&isscalar(x));
P.addParameter('STFTlen',           1024,       @(x) isnumeric(x)&&isscalar(x));
P.addParameter('UseGPU',            false,      @islogical);
P.addParameter('OutFile',           '',         @(s) ischar(s)||isstring(s));
P.parse(varargin{:});
opts = P.Results;

%% ---------- 1)  Load SOFA ---------------------------------------------
S   = SOFAload(sofaIn);                         % requires SOFAstart
IR  = double(S.Data.IR);                        % M × C × Ns
[M,C,Ns] = size(IR);
fs  = S.Data.SamplingRate;
N   = round(sqrt(C)-1);                         % HOA order (3-D)

%% ---------- 2)  Build steering vector Y (real, ACN/N3D) ---------------
az = deg2rad( azDeg );
el = deg2rad( elDeg );
Y  = getRealSH_N3D(N, az, el);                  % 1 × C

%% ---------- 3)  Derive per-order weights w_n --------------------------
switch lower(opts.Weighting)
    case 'maxre'
        w_n = getMaxRE_perOrder(N);             % 1 × (N+1)
    case 'gaussian'
        w_n = exp(-(0:N)/opts.Sigma);
    case 'custom'
        assert(~isempty(opts.CustomWeights),...
           'CustomWeights must be provided when Weighting=''custom''');
        w_n = opts.CustomWeights(:).';
        assert(numel(w_n)==N+1,'CustomWeights length must be N+1');
    otherwise % 'none'
        w_n = ones(1,N+1);
end
w_vec = expandWeights_ACN(w_n);                 % 1 × C
Yw    = Y .* w_vec;                             % steering + taper
Yw    = Yw / norm(Yw);                          % unit energy
Ycol  = Yw.';                                   % column for decoding

%% ---------- 4)  Cross-fade gains --------------------------------------
switch lower(opts.FadeLaw)
    case 'equalpower'
        g1 = cos(pi/2*opts.Beta);
        g2 = sin(pi/2*opts.Beta);
    otherwise % 'linear'
        g1 = 1-opts.Beta;  g2 = opts.Beta;
end

%% ---------- 5)  GPU? ---------------------------------------------------
if opts.UseGPU && gpuDeviceCount>0
    IR  = gpuArray(IR);
    Ycol= gpuArray(Ycol);
    Yw  = gpuArray(Yw);
end

%% ---------- 6)  Main processing path ----------------------------------
if ~opts.FreqMode          % ------- 6a)  Sample-domain (fast) ---------
    for m = 1:M
        X  = squeeze(IR(m,:,:)).';              % Ns×C
        s  = X * Ycol;                          % decode
        Xs = s * Yw;                            % re-encode
        Xo = g1 * X + g2 * Xs;                  % cross-fade
        IR(m,:,:) = permute(Xo.',[1 3 2]);      % C×Ns → 1×C×Ns
    end
    
else                       % ------- 6b)  STFT per-band steering -------
    L  = opts.STFTlen;     hop = L/2;
    win = hann(L,'periodic');
    [B,~,~] = stft(randn(Ns,1),fs,'Window',win,'OverlapLength',L-hop,'FFTLength',L);
    nBins = size(B,1);  f = (0:nBins-1).' * fs/L;
    hiBand = f >= opts.Fsplit;
    
    for m = 1:M
        Xf = squeeze(IR(m,:,:)).';              % Ns×C
        Sf = stft(Xf,'Window',win,'OverlapLength',L-hop,'FFTLength',L);
        % Sf:  nBins × nFrames × C
        Sdec = pagemtimes( permute(Sf,[3 1 2]), Ycol );      % 1 × nBins × nFrames
        Senc = pagemtimes( Sdec, Yw );                      % 1 × nBins × nFrames × C
        % Frequency-selective weighting:
        %   low band → leave untouched
        %   high band → apply cross-fade
        mask  = reshape(hiBand,[nBins 1 1]);
        Xf_hi = g1 * Sf .* mask + g2 * Senc .* mask;
        Xf_lo = Sf .* (~mask);
        Xf_out= Xf_lo + Xf_hi;
        
        Xo = istft(Xf_out,'Window',win,'OverlapLength',L-hop,'FFTLength',L, ...
                    'ConjugateSymmetric',true,'OutputLength',Ns);
        IR(m,:,:) = permute(Xo.',[1 3 2]);      % C×Ns → 1×C×Ns
    end
end

%% ---------- 7)  Cast + Save -------------------------------------------
if isa(IR,'gpuArray'), IR = gather(IR); end
S.Data.IR = single(IR);                         % back to single
if isempty(opts.OutFile)
    [pth,nm,~] = fileparts(sofaIn);
    opts.OutFile = fullfile(pth, ...
       sprintf('%s_MAXRE_%d_%d_beta%.2f.sofa', nm, round(elDeg), round(azDeg), opts.Beta));
end
SOFAsave(opts.OutFile, S);
fprintf('[steerAmbiSOFA_maxre]  Saved →  %s\n', opts.OutFile);
sofaOut = opts.OutFile;
end
% ================== Helper: closed-form max-rE weights ===================
function w_n = getMaxRE_perOrder(N)
try                     % prefer toolbox version if it exists
    w_n = getMaxREweights(N).';
catch                    % fallback to closed-form
    n = 0:N;
    w_n = sin((n+1)*pi/(2*N+2)) ./ ((n+1)*sin(pi/(2*N+2)));
end
end
% ================== Helper: real N3D SH (same as before) ================
function Y = getRealSH_N3D(N, az, el)
theta = pi/2 - el;
Y = zeros(1,(N+1)^2);
for n = 0:N
    P = legendre(n, cos(theta), 'sch');
    for m = -n:n
        idx = n*(n+1)+m+1;
        Pnm = P(abs(m)+1);
        K   = sqrt((2*n+1)/(4*pi) * factorial(n-abs(m))/factorial(n+abs(m)));
        if m>0,   Y(idx)=sqrt(2)*K*Pnm*cos(m*az);
        elseif m<0, Y(idx)=sqrt(2)*K*Pnm*sin(abs(m)*az);
        else,       Y(idx)=K*Pnm;
        end
    end
end
Y = Y/norm(Y);
end
% ================== Helper: expand per-order weights to ACN vector =======
function wvec = expandWeights_ACN(w_order)
N = numel(w_order)-1;  idx=1;  wvec = zeros(1,(N+1)^2);
for n=0:N, mCnt=2*n+1; wvec(idx:idx+mCnt-1)=w_order(n+1); idx=idx+mCnt; end
end
