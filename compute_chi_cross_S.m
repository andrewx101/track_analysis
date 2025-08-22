function [chi_cross,chi_sing,Qsum] = compute_chi_cross_S(S, DeltaT, delta)
            % S        : pre-processed cell of FOV structs (see earlier code)
            % DeltaT   : frame-lag
            % delta    : tolerance
            %
            % This version subtracts each track’s own mean instead of a global mean.

            nFov = numel(S);
            % — First pass: extract each track’s q_i(t) series, and its mean —
            TrackSeries = cell(nFov,1);
            for f = 1:nFov
                T   = S{f};
                K   = numel(T);
                ts  = cell(K,1);    % times
                qs  = cell(K,1);    % q-values
                for k = 1:K
                    tvec = T(k).t;
                    x    = T(k).x;
                    y    = T(k).y;

                    [Lia, LocB] = ismember(tvec+DeltaT, tvec);
                    ok = find(Lia);
                    if isempty(ok), continue; end

                    ti   = tvec(ok);
                    dx   = x(LocB(ok)) - x(ok);
                    dy   = y(LocB(ok)) - y(ok);
                    qval = exp( - (dx.^2 + dy.^2) / (2*delta^2) );
                    % qval=exp(-(dx.^2)/(2*delta^2));

                    ts{k} = ti;
                    qs{k} = qval;
                end

                % Compute per-track means
                qbar_i = nan(K,1);
                for k = 1:K
                    if isempty(qs{k})
                        qbar_i(k) = NaN;
                    else
                        qbar_i(k) = mean(qs{k});
                    end
                end

                TrackSeries{f}.tq   = ts;
                TrackSeries{f}.qq   = qs;
                TrackSeries{f}.qb   = qbar_i;
            end

            % — Compute chi_sing and Qsum from per-track stats —
            chi_sing = 0;
            Qsum     = 0;
            nTracks=0;
            for f = 1:nFov
                qs = TrackSeries{f}.qq;
                qb = TrackSeries{f}.qb;
                for k = 1:numel(qb)
                    qv = qs{k};
                    if isempty(qv), continue; end
                    mu1 = qb(k);           % E[q_i]
                    mu2 = mean(qv.^2);     % E[q_i^2]
                    chi_sing = chi_sing + (mu2 - mu1^2);
                    Qsum     = Qsum     + mu1;
                    nTracks=nTracks + 1;
                end
            end
            if nTracks>0
                Qsum=Qsum/nTracks;
                chi_sing=chi_sing/nTracks;
            else
                chi_sing=NaN;
                Qsum=NaN;
            end


            % — Second pass: optimized cross-covariance in O(#events) per FOV —
            num = 0;
            den = 0;
            for f = 1:nFov
                ts = TrackSeries{f}.tq;   % cell(K,1) of time-vectors
                qs = TrackSeries{f}.qq;   % cell(K,1) of q-values
                qb = TrackSeries{f}.qb;   % (K×1) per-track means
                K  = numel(ts);
                if K<2, continue; end

                % 1) Preallocate and flatten all (t,qdev) pairs in this FOV
                nE = sum(cellfun(@numel, ts));
                if nE<1, continue; end
                T_all    = zeros(nE,1);
                Qdev_all = zeros(nE,1);
                idx = 1;
                for k = 1:K
                    tks = ts{k};
                    if isempty(tks), continue; end
                    nk = numel(tks);
                    T_all(idx:idx+nk-1)    = tks;
                    Qdev_all(idx:idx+nk-1) = qs{k} - qb(k);
                    idx = idx + nk;
                end

                % 2) Group by time
                [~,~,ic] = unique(T_all);              % ic(i)=bin-index of T_all(i)
                cnts     = accumarray(ic, 1);          % number of tracks at each unique time
                s1       = accumarray(ic, Qdev_all);   % sum of qdev
                s2       = accumarray(ic, Qdev_all.^2);% sum of qdev^2

                % 3) Only times with ≥2 tracks contribute
                mask = cnts >= 2;
                c    = cnts(mask);
                s1   = s1(mask);
                s2   = s2(mask);

                % 4) Use ½[(∑qdev)^2 − ∑(qdev^2)] = ∑_{i<j} qdev_i qdev_j
                num = num + sum( (s1.^2 - s2) * 0.5 );
                den = den + sum( c .* (c - 1) * 0.5 );
            end

            chi_cross = num / den;
        end