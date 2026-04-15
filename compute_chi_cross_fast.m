function [chicrossraw, chicrossunb, chising, Qsum, varraw, varunb, biashat] = ...
    compute_chi_cross_fast(S, DeltaT, deltalist)
% COMPUTE_CHI_CROSS_FAST High-performance estimator for inter-particle dynamical correlation.
% Utilizes FFT exact-support variance. Designed to be called inside a Parfor loop over lags.

if nargin < 3
    error('Usage: compute_chi_cross_fast(S, DeltaT, deltalist)');
end

nDelta = numel(deltalist);
nFov = numel(S);

% --- 1. Precompute Grids & Single-Track Metrics per FOV ---
fovData = cell(nFov, 1);
totalPairs = 0;

for f = 1:nFov
    Tf = S{f};
    if isempty(Tf), continue; end

    [trks, nTracks] = local_extract_tracks(Tf);
    if nTracks < 2, continue; end

    [dr2_mat, min_t, dt_frame] = local_build_fov_grid(trks, nTracks, DeltaT);
    [f_chising, f_Qsum, f_nObsQ, f_nTrack] = local_compute_track_metrics(dr2_mat, deltalist);

    pairs = nchoosek(1:nTracks, 2);
    nP = size(pairs, 1);

    fovData{f} = struct('dr2_mat', dr2_mat, 'pairs', pairs, 'nP', nP, ...
        'chising', f_chising, 'Qsum', f_Qsum, ...
        'nObsQ', f_nObsQ, 'nTrack', f_nTrack);
    totalPairs = totalPairs + nP;
end

if totalPairs == 0
    chicrossraw = NaN(nDelta, 1); chicrossunb = NaN(nDelta, 1); chising = NaN(nDelta, 1);
    Qsum = NaN(nDelta, 1); varraw = NaN(nDelta, 1); varunb = NaN(nDelta, 1); biashat = NaN(nDelta, 1);
    return;
end

% --- 2. Flatten pairs for iteration ---
pairList = zeros(totalPairs, 3);
idx = 1;
for f = 1:nFov
    if isempty(fovData{f}), continue; end
    nP = fovData{f}.nP;
    pairList(idx:idx+nP-1, 1) = f;
    pairList(idx:idx+nP-1, 2:3) = fovData{f}.pairs;
    idx = idx + nP;
end

pair_raw    = NaN(totalPairs, nDelta);
pair_unb    = NaN(totalPairs, nDelta);
pair_bias   = NaN(totalPairs, nDelta);
pair_varraw = NaN(totalPairs, nDelta);
pair_varunb = NaN(totalPairs, nDelta);
pair_valid  = false(totalPairs, 1);

% --- 3. Compute exact pair-level overlaps (Sequential) ---
for p = 1:totalPairs
    f = pairList(p, 1);
    i = pairList(p, 2);
    j = pairList(p, 3);

    dr2_i = fovData{f}.dr2_mat(i, :);
    dr2_j = fovData{f}.dr2_mat(j, :);

    valid_mask = ~isnan(dr2_i) & ~isnan(dr2_j);
    N = sum(valid_mask);
    if N < 2, continue; end

    pair_valid(p) = true;
    dr2_i_v = dr2_i(valid_mask);
    dr2_j_v = dr2_j(valid_mask);
    scale = N / (N - 1);

    p_raw = zeros(1, nDelta); p_unb = zeros(1, nDelta); p_bias = zeros(1, nDelta);
    p_vraw = zeros(1, nDelta); p_vunb = zeros(1, nDelta);

    % >>> NEW: Precompute the fixed missing-data mask FFT outside the delta loop <<<
    L = length(valid_mask);
    is_fast_fft = L < 50000;
    if is_fast_fft
        Mf = zeros(L, 1);
        Mf(valid_mask) = 1;
        mhall = round(conv(Mf, flip(Mf)));
        posidx = L+1 : 2*L-1;
        mh = mhall(posidx);
        valid_h = mh > 0;
        mh_v = mh(valid_h);
        idx_v = posidx(valid_h);
    else
        % Variables for generic fallback
        mh_v = []; idx_v = []; Mf = [];
    end
    % >>> END NEW <<<

    for d = 1:nDelta
        delta = deltalist(d);
        qic = exp(-dr2_i_v / (2 * delta^2));
        qjc = exp(-dr2_j_v / (2 * delta^2));

        z_val = (qic - mean(qic)) .* (qjc - mean(qjc));
        rawij = mean(z_val);
        unbij = scale * rawij;

        p_raw(d)  = rawij;
        p_unb(d)  = unbij;
        p_bias(d) = unbij - rawij;

        % Pass the precomputed mask arrays to the updated variance function
        vrawij = local_exact_var_z_optimized(z_val, valid_mask, is_fast_fft, Mf, mh_v, idx_v);
        p_vraw(d) = vrawij;
        p_vunb(d) = scale^2 * vrawij;
    end


    pair_raw(p, :)    = p_raw;
    pair_unb(p, :)    = p_unb;
    pair_bias(p, :)   = p_bias;
    pair_varraw(p, :) = p_vraw;
    pair_varunb(p, :) = p_vunb;
end

% --- 4. Re-aggregate to FOV level ---
fov_raw    = NaN(nFov, nDelta);
fov_unb    = NaN(nFov, nDelta);
fov_bias   = NaN(nFov, nDelta);
fov_varraw = NaN(nFov, nDelta);
fov_varunb = NaN(nFov, nDelta);
fov_nPair  = zeros(nFov, 1);

for f = 1:nFov
    if isempty(fovData{f}), continue; end
    p_idx = find(pairList(:, 1) == f & pair_valid);
    nP_f = length(p_idx);
    fov_nPair(f) = nP_f;

    if nP_f > 0
        fov_raw(f, :)    = mean(pair_raw(p_idx, :), 1);
        fov_unb(f, :)    = mean(pair_unb(p_idx, :), 1);
        fov_bias(f, :)   = mean(pair_bias(p_idx, :), 1);
        fov_varraw(f, :) = sum(pair_varraw(p_idx, :), 1) / nP_f^2;
        fov_varunb(f, :) = sum(pair_varunb(p_idx, :), 1) / nP_f^2;
    end
end

% --- 5. Global Weighted Aggregations ---
[chicrossraw, chicrossunb, biashat, varraw, varunb] = local_aggregate_fovs(fov_raw, fov_unb, fov_bias, fov_varraw, fov_varunb, fov_nPair, nDelta);

fov_chising = NaN(nFov, nDelta); fov_Qsum = NaN(nFov, nDelta);
fov_nObs = zeros(nFov, nDelta); fov_nTrk = zeros(nFov, nDelta);
for f = 1:nFov
    if ~isempty(fovData{f})
        fov_chising(f, :) = fovData{f}.chising;
        fov_Qsum(f, :)    = fovData{f}.Qsum;
        fov_nObs(f, :)    = fovData{f}.nObsQ;
        fov_nTrk(f, :)    = fovData{f}.nTrack;
    end
end
[chising, Qsum] = local_aggregate_track_metrics(fov_chising, fov_Qsum, fov_nObs, fov_nTrk, nDelta);

end

% =========================================================================
% HELPER FUNCTIONS
% =========================================================================

function [trks, nTracks] = local_extract_tracks(Tf)
if iscell(Tf)
    nTracks = numel(Tf);
    trks = struct('t', cell(nTracks,1), 'r', cell(nTracks,1));
    for k = 1:nTracks
        trk = local_standardize_track(Tf{k});
        trks(k).t = trk.t; trks(k).r = trk.r;
    end
elseif isstruct(Tf)
    nTracks = numel(Tf);
    trks = struct('t', cell(nTracks,1), 'r', cell(nTracks,1));
    for k = 1:nTracks
        trk = local_standardize_track(Tf(k));
        trks(k).t = trk.t; trks(k).r = trk.r;
    end
else
    nTracks = 0; trks = [];
end
end

function trkStd = local_standardize_track(trk)
tvec = trk.t;
if isfield(trk, 'r')
    rmat = trk.r;
elseif isfield(trk, 'pos')
    rmat = trk.pos;
elseif isfield(trk, 'x') && isfield(trk, 'y')
    rmat = [trk.x, trk.y];
end
if size(rmat,1) ~= numel(tvec), rmat = rmat.'; end
[tu, ia] = unique(tvec, 'stable');
[~, ord] = sort(tu);
trkStd = struct('t', tu(ord), 'r', rmat(ia(ord), :));
end

function [dr2_mat, min_t, dt_frame] = local_build_fov_grid(trks, nTracks, DeltaT)
all_t = cell2mat({trks.t}');
tu = unique(all_t);
dt_all = diff(tu);
dt_frame = min(dt_all(dt_all > 1e-4));
if isempty(dt_frame), dt_frame = 1; end
min_t = min(tu);

L_grid = round((max(tu) - min_t) / dt_frame) + 1;
dr2_mat = NaN(nTracks, L_grid);

for k = 1:nTracks
    tvec = trks(k).t;
    rmat = trks(k).r;

    [Lia, LocB] = ismember(tvec + DeltaT, tvec);
    valid_dr = find(Lia);
    if isempty(valid_dr), continue; end

    dr2_val = sum((rmat(LocB(valid_dr), :) - rmat(valid_dr, :)).^2, 2);
    t_valid = tvec(valid_dr);

    idx = round((t_valid - min_t) / dt_frame) + 1;
    dr2_mat(k, idx) = dr2_val;
end
end

function [f_chising, f_Qsum, f_nObsQ, f_nTrack] = local_compute_track_metrics(dr2_mat, deltalist)
nDelta = numel(deltalist);
nTracks = size(dr2_mat, 1);
f_chising = NaN(1, nDelta); f_Qsum = NaN(1, nDelta);
f_nObsQ = zeros(1, nDelta); f_nTrack = zeros(1, nDelta);

for d = 1:nDelta
    delta = deltalist(d);
    qAll = [];
    trackVars = NaN(nTracks, 1);
    for k = 1:nTracks
        dr_k = dr2_mat(k, :);
        dr_k = dr_k(~isnan(dr_k));
        if isempty(dr_k), continue; end
        qk = exp(-dr_k / (2*delta^2));
        qAll = [qAll; qk(:)];
        trackVars(k) = mean((qk - mean(qk)).^2);
    end
    if ~isempty(qAll)
        f_chising(d) = mean(trackVars(~isnan(trackVars)));
        f_Qsum(d) = mean(qAll);
        f_nObsQ(d) = numel(qAll);
        f_nTrack(d) = sum(~isnan(trackVars));
    end
end
end

function v_out = local_exact_var_z_optimized(z_val, valid_mask, is_fast_fft, Mf, mh_v, idx_v)
    % EXACT-SUPPORT COVARIANCE: FFT convolution threshold methodology
    % Uses pre-calculated missing-data mask convolutions to save CPU time.
    N = sum(valid_mask);
    mz = mean(z_val);
    v0 = mean((z_val - mz).^2) / N;
    L = length(valid_mask);
    vcorr = 0;
    
    if is_fast_fft 
        Xf = zeros(L, 1); 
        Xf(valid_mask) = z_val;
        
        % Pad to next power of 2 for maximum FFT speed
        N_fft = 2^nextpow2(2*L - 1);
        
        % Explicit FFT for faster execution than built-in conv
        Xf_fft = fft(Xf, N_fft);
        Mf_fft = fft(Mf, N_fft);
        Xf_rev_fft = fft(flip(Xf), N_fft);
        Mf_rev_fft = fft(flip(Mf), N_fft);
        
        % Compute only the variables that depend on z_val
        sumAx_full = ifft(Mf_fft .* Xf_rev_fft);
        sumBx_full = ifft(Xf_fft .* Mf_rev_fft);
        gxx_full   = ifft(Xf_fft .* Xf_rev_fft);
        
        if any(mh_v > 0)
            % Extract the valid positive-offset slices
            meanAx = sumAx_full(idx_v) ./ mh_v;
            meanBx = sumBx_full(idx_v) ./ mh_v;
            Gh = gxx_full(idx_v) ./ mh_v - meanAx .* meanBx;
            
            vcorr = sum(2 * mh_v / N^2 .* Gh);
        end
    else 
        % Generic fallback path for ultra-sparse segments
        steps = find(valid_mask);
        offsets = unique(bsxfun(@minus, steps, steps')); 
        offsets = offsets(offsets > 0);
        for h = offsets(:)'
            [tf, loc] = ismember(steps + h, steps);
            idxA = find(tf);
            idxB = loc(tf);
            m = length(idxA);
            if m > 0
                A = z_val(idxA); B = z_val(idxB);
                ma = mean(A); mb = mean(B);
                Gh = (A - ma)' * (B - mb) / m;
                vcorr = vcorr + 2 * m / N^2 * Gh;
            end
        end
    end
    v_out = max(real(v0 + vcorr), 0);
end

function [out_raw, out_unb, out_bias, out_varraw, out_varunb] = local_aggregate_fovs(fov_raw, fov_unb, fov_bias, fov_varraw, fov_varunb, fov_nPair, nDelta)
    out_raw    = NaN(nDelta, 1); 
    out_unb    = NaN(nDelta, 1); 
    out_bias   = NaN(nDelta, 1);
    out_varraw = NaN(nDelta, 1); 
    out_varunb = NaN(nDelta, 1);
    
    valid_fovs = fov_nPair > 0; % Full length mask (e.g., length 3 or 4)
    if ~any(valid_fovs)
        return; 
    end
    
    wPair = fov_nPair(valid_fovs);
    wPair_norm = wPair / sum(wPair);
    
    for d = 1:nDelta
        r_col = fov_raw(valid_fovs, d);
        valid = ~isnan(r_col); % Sub-mask (shorter length)
        
        if any(valid)
            % Create a fresh full-length mask of mostly false
            current_valid_mask = false(size(valid_fovs));
            % Inject the true values back exactly where they belong!
            current_valid_mask(valid_fovs) = valid;
            
            w = wPair_norm(valid) / sum(wPair_norm(valid));
            
            % Now we safely index using a mask that is guaranteed to be length nFov
            out_raw(d)    = sum(w .* r_col(valid));
            out_unb(d)    = sum(w .* fov_unb(current_valid_mask, d));
            out_bias(d)   = sum(w .* fov_bias(current_valid_mask, d));
            out_varraw(d) = sum(w.^2 .* fov_varraw(current_valid_mask, d));
            out_varunb(d) = sum(w.^2 .* fov_varunb(current_valid_mask, d));
        end
    end
end

function [out_chising, out_Qsum] = local_aggregate_track_metrics(fov_chising, fov_Qsum, fov_nObs, fov_nTrk, nDelta)
out_chising = NaN(nDelta, 1); out_Qsum = NaN(nDelta, 1);
for d = 1:nDelta
    val_trk = fov_nTrk(:, d) > 0 & ~isnan(fov_chising(:, d));
    if any(val_trk)
        wT = fov_nTrk(val_trk, d) / sum(fov_nTrk(val_trk, d));
        out_chising(d) = sum(wT .* fov_chising(val_trk, d));
    end
    val_obs = fov_nObs(:, d) > 0 & ~isnan(fov_Qsum(:, d));
    if any(val_obs)
        wO = fov_nObs(val_obs, d) / sum(fov_nObs(val_obs, d));
        out_Qsum(d) = sum(wO .* fov_Qsum(val_obs, d));
    end
end
end