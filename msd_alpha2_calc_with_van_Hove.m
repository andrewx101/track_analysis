function out = msd_alpha2_calc(tr_list, dt, res, dim)
% msd_alpha2_calc
%
% Combined version of msd_calc and alpha2_calc using one shared exact-support
% displacement preparation step.
%
% INPUT
% tr_list : cell array of FOVs; each FOV contains trajectory cells.
% Each trajectory matrix must have columns
% [x1, ..., xd, step_id, track_id]
% dt : physical time step between successive frame indices
% res : spatial resolution multiplier
% dim : spatial dimension d
%
% OUTPUT
% Union of the original outputs from msd_calc and alpha2_calc.
% Field names are kept as close as possible to the originals.

d = dim;
K = (d + 2) / d;
nFov = numel(tr_list);

%% Flatten trajectories into a shared track list
tracks = struct( ...
    'fov', {}, 'track_id', {}, 'steps', {}, 'pos', {}, ...
    'ni', {}, 'ell', {}, ...
    'lag_frames_msd', {}, 'lag_msd', {}, 'msd', {}, 'var_hac', {}, 'q_msd', {}, ...
    'lag_frames_alpha2', {}, 'lag_alpha2', {}, ...
    'X', {}, 'Y', {}, ...
    'varX_full', {}, 'varY_full', {}, 'covXY_full', {}, ...
    'varX_corr', {}, 'varY_corr', {}, 'covXY_corr', {}, ...
    'alpha2_raw', {}, 'alpha2_norm', {}, 'alpha2_bc', {}, ...
    'var_alpha2', {}, 'bias_corr_raw', {}, 'bias_corr_bc', {}, ...
    'q_alpha2', {});

for f = 1:nFov
    trf = tr_list{f};
    for j = 1:numel(trf)
        t = trf{j};
        if isempty(t)
            continue;
        end
        if size(t,2) ~= d + 2
            error('Each trajectory must have d position columns + step_id + track_id.');
        end

        t = sortrows(t, d+1);
        steps = t(:, d+1);
        pos = t(:, 1:d) * res;
        tr_id = t(1, d+2);

        if numel(steps) < 2
            continue;
        end

        s.fov = f;
        s.track_id = tr_id;
        s.steps = steps(:);
        s.pos = pos;
        s.ni = steps(end) - steps(1);
        s.ell = numel(steps);

        s.lag_frames_msd = [];
        s.lag_msd = [];
        s.msd = [];
        s.var_hac = [];
        s.q_msd = [];

        s.lag_frames_alpha2 = [];
        s.lag_alpha2 = [];
        s.X = [];
        s.Y = [];
        s.varX_full = [];
        s.varY_full = [];
        s.covXY_full = [];
        s.varX_corr = [];
        s.varY_corr = [];
        s.covXY_corr = [];
        s.alpha2_raw = [];
        s.alpha2_norm = [];
        s.alpha2_bc = [];
        s.var_alpha2 = [];
        s.bias_corr_raw = [];
        s.bias_corr_bc = [];
        s.q_alpha2 = [];

        tracks(end+1) = s; %#ok<AGROW>
    end
end

nTracks = numel(tracks);
if nTracks == 0
    error('No valid trajectories found.');
end

all_ni = [tracks.ni]';
total_ni_all = sum(all_ni);
maxLag = max(all_ni);
if ~isfinite(maxLag) || maxLag < 1
    error('Not enough span across tracks to compute any lag.');
end

%% Per-track outputs
msd_i = cell(nTracks, 1);
alpha2_i = cell(nTracks, 1);
van_Hove_dx = cell(maxLag, 1);

for i = 1:nTracks
    steps = tracks(i).steps;
    pos = tracks(i).pos;
    tr_id = tracks(i).track_id;
    fov_i = tracks(i).fov;
    ni_i = tracks(i).ni;
    ell_i = tracks(i).ell;

    lag_frames_msd = [];
    lag_time_msd = [];
    msd_vals = [];
    var_vals = [];
    q_vals_msd = [];

    lag_frames_alpha2 = [];
    lag_time_alpha2 = [];
    Xvals = [];
    Yvals = [];
    varXfullVals = [];
    varYfullVals = [];
    covXYfullVals = [];
    varXcorrVals = [];
    varYcorrVals = [];
    covXYcorrVals = [];
    aRawVals = [];
    aNormVals = [];
    aBcVals = [];
    varAvals = [];
    bCorrRawVals = [];
    bCorrBcVals = [];
    qVals_alpha2 = [];

    maxLag_i = min(maxLag, ni_i);

    for n = 1:maxLag_i
        [r2, start_steps, dx] = local_sqdisp_series_exactsupport(steps, pos, n);
        q = numel(r2);

        if q == 0
            continue;
        end

        x = r2(:);

        %% MSD branch
        lag_frames_msd(end+1,1) = n; %#ok<AGROW>
        lag_time_msd(end+1,1) = n * dt; %#ok<AGROW>
        msd_vals(end+1,1) = mean(x); %#ok<AGROW>
        var_vals(end+1,1) = exact_var_of_mean_from_support(x, start_steps); %#ok<AGROW>
        q_vals_msd(end+1,1) = q; %#ok<AGROW>

        %% alpha2 branch
        if q < 2
            continue;
        end

        y = x.^2;
        X = mean(x);
        Y = mean(y);

        if ~(isfinite(X) && isfinite(Y)) || X <= 0
            continue;
        end

        Vdec = exact_vcov_of_means_exactsupport([x, y], start_steps);

        varX_full = Vdec.full(1,1);
        varY_full = Vdec.full(2,2);
        covXY_full = Vdec.full(1,2);

        varX_corr = Vdec.corr(1,1);
        varY_corr = Vdec.corr(2,2);
        covXY_corr = Vdec.corr(1,2);

        if any(~isfinite([varX_full, varY_full, covXY_full, ...
                          varX_corr, varY_corr, covXY_corr]))
            continue;
        end

        aRaw = Y / (K * X^2) - 1;

        cNorm = (q * d + 2) / (q * d);
        aNorm = cNorm * (aRaw + 1) - 1;

        bCorrRaw = (3 * Y / (K * X^4)) * varX_corr ...
                 - (2 / (K * X^3)) * covXY_corr;
        bCorrBc = cNorm * bCorrRaw;
        aBc = aNorm - bCorrBc;

        fx = -2 * Y / (K * X^3);
        fy = 1 / (K * X^2);
        varRaw = fx^2 * varX_full + fy^2 * varY_full + 2 * fx * fy * covXY_full;
        varA = cNorm^2 * max(real(varRaw), 0);

        lag_frames_alpha2(end+1,1) = n; %#ok<AGROW>
        lag_time_alpha2(end+1,1) = n * dt; %#ok<AGROW>
        Xvals(end+1,1) = X; %#ok<AGROW>
        Yvals(end+1,1) = Y; %#ok<AGROW>
        varXfullVals(end+1,1) = varX_full; %#ok<AGROW>
        varYfullVals(end+1,1) = varY_full; %#ok<AGROW>
        covXYfullVals(end+1,1) = covXY_full; %#ok<AGROW>
        varXcorrVals(end+1,1) = varX_corr; %#ok<AGROW>
        varYcorrVals(end+1,1) = varY_corr; %#ok<AGROW>
        covXYcorrVals(end+1,1) = covXY_corr; %#ok<AGROW>
        aRawVals(end+1,1) = aRaw; %#ok<AGROW>
        aNormVals(end+1,1) = aNorm; %#ok<AGROW>
        aBcVals(end+1,1) = aBc; %#ok<AGROW>
        varAvals(end+1,1) = varA; %#ok<AGROW>
        bCorrRawVals(end+1,1) = bCorrRaw; %#ok<AGROW>
        bCorrBcVals(end+1,1) = bCorrBc; %#ok<AGROW>
        qVals_alpha2(end+1,1) = q; %#ok<AGROW>
    end

    %% Store MSD track data
    tracks(i).lag_frames_msd = lag_frames_msd;
    tracks(i).lag_msd = lag_time_msd;
    tracks(i).msd = msd_vals;
    tracks(i).var_hac = var_vals;
    tracks(i).q_msd = q_vals_msd;

    msd_i{i} = [lag_frames_msd, lag_time_msd, msd_vals, var_vals, q_vals_msd, ...
                repmat(tr_id, numel(lag_frames_msd), 1), ...
                repmat(fov_i, numel(lag_frames_msd), 1), ...
                repmat(ni_i, numel(lag_frames_msd), 1), ...
                repmat(ell_i, numel(lag_frames_msd), 1)];

    %% Store alpha2 track data
    tracks(i).lag_frames_alpha2 = lag_frames_alpha2;
    tracks(i).lag_alpha2 = lag_time_alpha2;
    tracks(i).X = Xvals;
    tracks(i).Y = Yvals;
    tracks(i).varX_full = varXfullVals;
    tracks(i).varY_full = varYfullVals;
    tracks(i).covXY_full = covXYfullVals;
    tracks(i).varX_corr = varXcorrVals;
    tracks(i).varY_corr = varYcorrVals;
    tracks(i).covXY_corr = covXYcorrVals;
    tracks(i).alpha2_raw = aRawVals;
    tracks(i).alpha2_norm = aNormVals;
    tracks(i).alpha2_bc = aBcVals;
    tracks(i).var_alpha2 = varAvals;
    tracks(i).bias_corr_raw = bCorrRawVals;
    tracks(i).bias_corr_bc = bCorrBcVals;
    tracks(i).q_alpha2 = qVals_alpha2;

    alpha2_i{i} = [lag_frames_alpha2, lag_time_alpha2, ...
                   aRawVals, aBcVals, varAvals, ...
                   qVals_alpha2, ...
                   repmat(tr_id, numel(lag_frames_alpha2), 1), ...
                   repmat(fov_i, numel(lag_frames_alpha2), 1), ...
                   repmat(ni_i, numel(lag_frames_alpha2), 1), ...
                   Xvals, Yvals, ...
                   varXfullVals, varYfullVals, covXYfullVals, ...
                   varXcorrVals, varYcorrVals, covXYcorrVals, ...
                   bCorrRawVals, bCorrBcVals, aNormVals];
end

%% MSD aggregation
lag_frames_all = (1:maxLag)';
lag_time_all = lag_frames_all * dt;

msd_naive = NaN(maxLag,1);
var_naive = NaN(maxLag,1);

for n = 1:maxLag
    mi = NaN(nTracks,1);
    vi = NaN(nTracks,1);
    qi = NaN(nTracks,1);
    elig = false(nTracks,1);

    for i = 1:nTracks
        idx = find(tracks(i).lag_frames_msd == n, 1, 'first');
        if ~isempty(idx)
            elig(i) = true;
            mi(i) = tracks(i).msd(idx);
            vi(i) = tracks(i).var_hac(idx);
            qi(i) = tracks(i).q_msd(idx);
        end
    end

    if ~any(elig)
        continue;
    end

    mi_e = mi(elig);
    vi_e = vi(elig);
    qi_e = qi(elig);

    sb = sum(qi_e);
    if isfinite(sb) && sb > 0
        b = qi_e / sb;
        msd_naive(n) = sum(b .* mi_e);

        if all(isfinite(vi_e))
            var_naive(n) = sum((b.^2) .* vi_e);
            var_naive(n) = max(real(var_naive(n)), 0);
        end
    end
end

%% alpha2 reference total q at lag 1
total_q1_all = 0;
for i = 1:nTracks
    idx1 = find(tracks(i).lag_frames_alpha2 == 1, 1, 'first');
    if ~isempty(idx1)
        total_q1_all = total_q1_all + tracks(i).q_alpha2(idx1);
    end
end
if total_q1_all <= 0
    total_q1_all = NaN;
end

%% alpha2 aggregation
alpha2_sp_raw = NaN(maxLag,1);
alpha2_sp_norm = NaN(maxLag,1);
alpha2_sp_bc = NaN(maxLag,1);
var_sp = NaN(maxLag,1);

alpha2_mp_raw = NaN(maxLag,1);
alpha2_mp_norm = NaN(maxLag,1);
alpha2_mp_ftbc = NaN(maxLag,1);
alpha2_mp_final = NaN(maxLag,1);
var_mp_final = NaN(maxLag,1);

n_elig_tracks = zeros(maxLag,1);
q_total = NaN(maxLag,1);

theta_ni = NaN(maxLag,1);
theta_q = NaN(maxLag,1);

fov_w = NaN(maxLag, nFov);
fov_w_q = NaN(maxLag, nFov);

ni_all = [tracks.ni]';
fi_all = [tracks.fov]';

for n = 1:maxLag
    elig = false(nTracks,1);

    aRaw = NaN(nTracks,1);
    aNorm = NaN(nTracks,1);
    aBc = NaN(nTracks,1);
    vA = NaN(nTracks,1);
    X = NaN(nTracks,1);
    Y = NaN(nTracks,1);
    vXfull = NaN(nTracks,1);
    vYfull = NaN(nTracks,1);
    cXYfull = NaN(nTracks,1);
    vXcorr = NaN(nTracks,1);
    vYcorr = NaN(nTracks,1);
    cXYcorr = NaN(nTracks,1);
    qv = NaN(nTracks,1);

    for i = 1:nTracks
        idx = find(tracks(i).lag_frames_alpha2 == n, 1, 'first');
        if ~isempty(idx)
            elig(i) = true;
            aRaw(i) = tracks(i).alpha2_raw(idx);
            aNorm(i) = tracks(i).alpha2_norm(idx);
            aBc(i) = tracks(i).alpha2_bc(idx);
            vA(i) = tracks(i).var_alpha2(idx);
            X(i) = tracks(i).X(idx);
            Y(i) = tracks(i).Y(idx);
            vXfull(i) = tracks(i).varX_full(idx);
            vYfull(i) = tracks(i).varY_full(idx);
            cXYfull(i) = tracks(i).covXY_full(idx);
            vXcorr(i) = tracks(i).varX_corr(idx);
            vYcorr(i) = tracks(i).varY_corr(idx);
            cXYcorr(i) = tracks(i).covXY_corr(idx);
            qv(i) = tracks(i).q_alpha2(idx);
        end
    end

    if ~any(elig)
        continue;
    end

    ni_e = ni_all(elig);
    fi_e = fi_all(elig);
    aRaw_e = aRaw(elig);
    aNorm_e = aNorm(elig);
    aBc_e = aBc(elig);
    vA_e = vA(elig);
    X_e = X(elig);
    Y_e = Y(elig);
    vXf_e = vXfull(elig);
    vYf_e = vYfull(elig);
    cXYf_e = cXYfull(elig);
    vXc_e = vXcorr(elig);
    vYc_e = vYcorr(elig);
    cXYc_e = cXYcorr(elig);
    q_e = qv(elig);

    n_elig_tracks(n) = numel(ni_e);
    q_total(n) = sum(q_e);

    sw_sp = sum(ni_e);
    if isfinite(sw_sp) && sw_sp > 0
        w_sp = ni_e / sw_sp;

        alpha2_sp_raw(n) = sum(w_sp .* aRaw_e);
        alpha2_sp_norm(n) = sum(w_sp .* aNorm_e);
        alpha2_sp_bc(n) = sum(w_sp .* aBc_e);

        if all(isfinite(vA_e))
            var_sp(n) = sum((w_sp.^2) .* vA_e);
            var_sp(n) = max(real(var_sp(n)), 0);
        end
    end

    [a_mp_bc, var_mp_bc, a_mp_norm, a_mp_raw] = weighted_mp_normcorr( ...
        X_e, Y_e, vXf_e, vYf_e, cXYf_e, vXc_e, cXYc_e, q_e, K, d);

    alpha2_mp_raw(n) = a_mp_raw;
    alpha2_mp_norm(n) = a_mp_norm;
    alpha2_mp_ftbc(n) = a_mp_bc;

    if numel(q_e) >= 2 && isfinite(a_mp_bc)
        jk_vals = NaN(numel(q_e),1);

        for j = 1:numel(q_e)
            keep = true(numel(q_e),1);
            keep(j) = false;

            [a_jk, ~] = weighted_mp_normcorr( ...
                X_e(keep), Y_e(keep), ...
                vXf_e(keep), vYf_e(keep), cXYf_e(keep), ...
                vXc_e(keep), cXYc_e(keep), q_e(keep), K, d);

            jk_vals(j) = a_jk;
        end

        if all(isfinite(jk_vals))
            m = numel(jk_vals);
            jk_mean = mean(jk_vals);
            alpha2_mp_final(n) = m * a_mp_bc - (m - 1) * jk_mean;
            var_mp_final(n) = (m - 1) / m * sum((jk_vals - jk_mean).^2);
            var_mp_final(n) = max(real(var_mp_final(n)), 0);
        else
            alpha2_mp_final(n) = a_mp_bc;
            var_mp_final(n) = var_mp_bc;
        end
    else
        alpha2_mp_final(n) = a_mp_bc;
        var_mp_final(n) = var_mp_bc;
    end

    if isfinite(total_ni_all) && total_ni_all > 0
        theta_ni(n) = sum(ni_e) / total_ni_all;
    end

    if isfinite(total_q1_all) && total_q1_all > 0
        theta_q(n) = sum(q_e) / total_q1_all;
    end

    if isfinite(sum(ni_e)) && sum(ni_e) > 0
        for f = 1:nFov
            fov_w(n,f) = sum(ni_e(fi_e == f)) / sum(ni_e);
        end
    end

    if isfinite(sum(q_e)) && sum(q_e) > 0
        for f = 1:nFov
            fov_w_q(n,f) = sum(q_e(fi_e == f)) / sum(q_e);
        end
    end
end

%% Build van Hove output
van_Hove = cell(maxLag, 2);
for n = 1:maxLag
    if ~isempty(van_Hove_dx{n})
        [N_vh, edges_vh] = histcounts(van_Hove_dx{n}, 100, 'Normalization', 'pdf');
        van_Hove{n,1} = N_vh;
        van_Hove{n,2} = edges_vh;
    end
end

%% Shared track info
track_info = table((1:nTracks)', [tracks.fov]', [tracks.track_id]', ...
                   [tracks.ni]', [tracks.ell]', ...
                   'VariableNames', {'track_index','fov','track_id', ...
                                     'duration_frames','observed_positions'});

%% Output
out = struct();

% msd_calc outputs
out.msd_i = msd_i;
out.msd_naive = [lag_time_all, msd_naive, var_naive];
out.track_info = track_info;
out.variance_definition = [ ...
    'Per-track variance: exact-support covariance-style estimator for the sample mean ', ...
    'of squared displacements using the actual supported lag-n start-time offsets derived from the data. ', ...
    'Naive pooled variance: weighted sum of per-track variances with weights proportional to q_i(n), assuming independence across tracks.' ];

% alpha2_calc outputs
out.alpha2_i = alpha2_i;
out.alpha2_sp = [lag_time_all, alpha2_sp_raw, alpha2_sp_bc, var_sp];
out.alpha2_sp_norm = [lag_time_all, alpha2_sp_norm];
out.alpha2_mp = [lag_time_all, alpha2_mp_raw, alpha2_mp_final, var_mp_final, n_elig_tracks];
out.alpha2_mp_norm = [lag_time_all, alpha2_mp_norm];
out.alpha2_mp_ftbc = [lag_time_all, alpha2_mp_ftbc];
out.theta = [lag_time_all, theta_ni];
out.theta_q = [lag_time_all, theta_q];
out.fov_w = fov_w;
out.fov_w_q = fov_w_q;
out.q_total = [lag_time_all, q_total];
out.K = K;
out.sampling = 'exact-support';
out.sp_weighting = 'duration-weighted';
out.mp_weighting = 'q-weighted pooled moments';
out.bc_definition = ['exact finite-q normalization plus correlation-only ', ...
                     'positive-offset exact-support bias term'];
out.alpha2_i_columns = {'lag_frame','lag','alpha2_raw_i','alpha2_bc_i','var_alpha2_i', ...
                        'q_i','track_id','fov_id','n_i','X_i','Y_i', ...
                        'varX_full_i','varY_full_i','covXY_full_i', ...
                        'varX_corr_i','varY_corr_i','covXY_corr_i', ...
                        'bias_corr_raw_i','bias_corr_bc_i','alpha2_norm_i'};
out.van_Hove = van_Hove;

end

function [a_bc, var_a, a_norm, a_raw] = weighted_mp_normcorr( ...
    X, Y, vX_full, vY_full, cXY_full, vX_corr, cXY_corr, q, K, d)

a_bc = NaN;
var_a = NaN;
a_norm = NaN;
a_raw = NaN;

if isempty(X) || isempty(Y) || isempty(q)
    return;
end

if any(~isfinite(X)) || any(~isfinite(Y)) || ...
   any(~isfinite(vX_full)) || any(~isfinite(vY_full)) || any(~isfinite(cXY_full)) || ...
   any(~isfinite(vX_corr)) || any(~isfinite(cXY_corr)) || any(~isfinite(q))
    return;
end

Q = sum(q);
if ~(isfinite(Q) && Q > 0)
    return;
end

w = q(:) / Q;

Xw = sum(w .* X(:));
Yw = sum(w .* Y(:));

if ~(isfinite(Xw) && isfinite(Yw)) || Xw <= 0
    return;
end

a_raw = Yw / (K * Xw^2) - 1;

cNorm = (Q * d + 2) / (Q * d);
a_norm = cNorm * (a_raw + 1) - 1;

VarXw_corr = sum((w.^2) .* vX_corr(:));
CovXYw_corr = sum((w.^2) .* cXY_corr(:));

bias_corr_raw = (3 * Yw / (K * Xw^4)) * VarXw_corr ...
              - (2 / (K * Xw^3)) * CovXYw_corr;

a_bc = a_norm - cNorm * bias_corr_raw;

VarXw_full = sum((w.^2) .* vX_full(:));
VarYw_full = sum((w.^2) .* vY_full(:));
CovXYw_full = sum((w.^2) .* cXY_full(:));

fx = -2 * Yw / (K * Xw^3);
fy = 1 / (K * Xw^2);
var_raw = fx^2 * VarXw_full + fy^2 * VarYw_full + 2 * fx * fy * CovXYw_full;
var_a = cNorm^2 * max(real(var_raw), 0);

end

function S = exact_vcov_of_means_exactsupport(Z, steps)

Z = Z(:,:);
steps = steps(:);
q = size(Z,1);
p = size(Z,2);

if numel(steps) ~= q
    error('exact_vcov_of_means_exactsupport:SizeMismatch', ...
          'steps must have one entry per row of Z.');
end

if q <= 1
    S = struct('full', NaN(p,p), 'zero', NaN(p,p), 'corr', NaN(p,p));
    return;
end

G0 = cov_centered_block(Z, Z);
V0 = G0 / q;
Vcorr = zeros(p,p);

offsets = unique(steps(:) - steps(:)');
offsets = offsets(offsets > 0);

for k = 1:numel(offsets)
    h = offsets(k);
    [idxA, idxB] = supported_pairs_by_offset(steps, h);
    m = numel(idxA);
    if m == 0
        continue;
    end

    A = Z(idxA, :);
    B = Z(idxB, :);
    G = cov_centered_block(A, B);
    term = (m / q^2) * (G + G');
    Vcorr = Vcorr + term;
end

Vfull = V0 + Vcorr;

S = struct();
S.full = (Vfull + Vfull') / 2;
S.zero = (V0 + V0') / 2;
S.corr = (Vcorr + Vcorr') / 2;

end

function G = cov_centered_block(A, B)
ma = mean(A, 1);
mb = mean(B, 1);
G = (A' * B) / size(A,1) - (ma' * mb);
end

function [r2, start_steps, dx] = local_sqdisp_series_exactsupport(steps, pos, n)

steps = steps(:);
target_steps = steps + n;
[tf, loc] = ismember(target_steps, steps);

idx1 = find(tf);
if isempty(idx1)
    r2 = [];
    start_steps = [];
    dx = [];
    return;
end

idx2 = loc(tf);
start_steps = steps(idx1);
[start_steps, ord] = sort(start_steps);
idx1 = idx1(ord);
idx2 = idx2(ord);

dr = pos(idx2, :) - pos(idx1, :);
dx = dr(:,1);
r2 = sum(dr.^2, 2);

end

function [idxA, idxB] = supported_pairs_by_offset(steps, h)
steps = steps(:);
[tf, loc] = ismember(steps + h, steps);
idxA = find(tf);
idxB = loc(tf);
end

function vmean = exact_var_of_mean_from_support(x, start_steps)

x = x(:);
start_steps = start_steps(:);

q = numel(x);
if q <= 1
    vmean = NaN;
    return;
end

if numel(start_steps) ~= q
    error('x and start_steps must have the same length.');
end

mhat = mean(x);
z = x - mhat;

D = bsxfun(@minus, start_steps.', start_steps);
C = z * z.';

lags = unique(D(:));
accum = 0;

for k = 1:numel(lags)
    h = lags(k);
    mask = (D == h);
    cnt = nnz(mask);

    if cnt == 0
        continue;
    end

    gamma_h = sum(C(mask)) / cnt;
    accum = accum + cnt * gamma_h;
end

vmean = accum / (q^2);
vmean = max(real(vmean), 0);

end