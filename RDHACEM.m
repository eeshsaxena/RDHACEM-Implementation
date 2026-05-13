% ==========================================================================
% RDHACEM.m
% Reversible Data Hiding with Automatic Contrast Enhancement for Medical Image
%
% Paper: Gao G., Tong S., Xia Z., Wu B., Xu L., Zhao Z.
%        Signal Processing, Vol. 178, 107817, January 2021.
%        DOI: 10.1016/j.sigpro.2020.107817
%
% Algorithm overview:
%   Separates medical image into ROI and NROI.
%   Automatically stretches ROI histogram to [0,255] to enlarge embedding
%   capacity and enhance contrast simultaneously.
%   Embeds secret bits into stretched ROI via histogram shifting at peak.
%   No brightness preservation (over-enhancement occurs at high capacity).
%
% Run:  RDHACEM     (full demo — 4 medical images, 4 embedding capacities)
% ==========================================================================
function RDHACEM()
    clc; close all;
    fprintf('=== RDHACEM: Automatic Contrast Enhancement for Medical Images ===\n');
    fprintf('    Gao et al., Signal Processing 178 (2021) 107817\n\n');

    imgs  = generate_medical_images();
    names = {'Brain01','Brain02','chest','xray'};
    cap_vals = [5000 10000 20000 50000];

    % Experiment 1: PSNR vs capacity
    fprintf('\n--- Exp 1: PSNR (dB) ---\n');
    fprintf('%-10s', 'Image');
    for c = cap_vals, fprintf('%12d', c); end; fprintf('\n');
    fprintf('%s\n', repmat('-',1,58));
    for k = 1:numel(names)
        fprintf('%-10s', names{k});
        for cap = cap_vals
            rng(42); pay = randi([0 1],1,cap,'uint8');
            [I_emb,~] = rdhacem_embed(imgs{k}, pay);
            fprintf('%12.2f', compute_psnr(imgs{k}, I_emb));
        end
        fprintf('\n');
    end

    % Experiment 2: Contrast improvement (ΔSD in ROI)
    fprintf('\n--- Exp 2: ΔSD in ROI (contrast gain) ---\n');
    fprintf('%-10s', 'Image');
    for c = cap_vals, fprintf('%12d', c); end; fprintf('\n');
    fprintf('%s\n', repmat('-',1,58));
    for k = 1:numel(names)
        fprintf('%-10s', names{k});
        for cap = cap_vals
            rng(42); pay = randi([0 1],1,cap,'uint8');
            I = imgs{k};
            [I_emb, meta] = rdhacem_embed(I, pay);
            roi = meta.roi_mask;
            ds = std(double(I_emb(roi))) - std(double(I(roi)));
            fprintf('%12.2f', ds);
        end
        fprintf('\n');
    end

    % Experiment 3: Brightness change (shows over-enhancement)
    fprintf('\n--- Exp 3: Brightness change |B - B_emb| (shows over-enhancement) ---\n');
    fprintf('%-10s', 'Image');
    for c = cap_vals, fprintf('%12d', c); end; fprintf('\n');
    fprintf('%s\n', repmat('-',1,58));
    for k = 1:numel(names)
        fprintf('%-10s', names{k});
        for cap = cap_vals
            rng(42); pay = randi([0 1],1,cap,'uint8');
            I = imgs{k};
            [I_emb, meta] = rdhacem_embed(I, pay);
            roi = meta.roi_mask;
            dB = abs(mean(double(I_emb(roi))) - mean(double(I(roi))));
            fprintf('%12.4f', dB);
        end
        fprintf('\n');
    end

    % Experiment 4: Reversibility
    fprintf('\n--- Exp 4: Reversibility (20000 bits) ---\n');
    for k = 1:numel(names)
        rng(42); pay = randi([0 1],1,20000,'uint8');
        [I_emb, meta] = rdhacem_embed(imgs{k}, pay);
        [I_rec, D_ext] = rdhacem_extract(I_emb, meta);
        ok   = isequal(imgs{k}, I_rec);
        n_ext = min(numel(D_ext), numel(pay));
        errs  = sum(D_ext(1:n_ext) ~= pay(1:n_ext));
        fprintf('  %-10s → Reversible: %s | Bit errors: %d\n', names{k}, string(ok), errs);
    end
    fprintf('\nDone.\n');
end

% ==========================================================================
%  STAGE 1: ROI / NROI SEGMENTATION  (Sec. 2, from [13])
% ==========================================================================
function roi_mask = get_roi_mask(I)
    thr      = graythresh(I) * 255;
    roi_mask = double(I) > thr;
    if sum(roi_mask(:)) < 0.05 * numel(I)
        [r,c] = size(I);
        mr = round(r*0.2); mc = round(c*0.2);
        roi_mask = false(r,c);
        roi_mask(mr:r-mr, mc:c-mc) = true;
    end
end

% ==========================================================================
%  STAGE 2: AUTOMATIC CONTRAST STRETCHING  (Sec. 2)
% ==========================================================================
function [I_str, I_MIN, I_MAX] = auto_contrast_stretch(I, roi_mask)
% Automatically stretches ROI pixel range [I_MIN, I_MAX] to [0, 255].
% This is the key step that distinguishes RDHACEM — no user-specified range.
    roi_pix = double(I(roi_mask));
    I_MIN   = min(roi_pix);
    I_MAX   = max(roi_pix);
    I_str   = I;
    if I_MAX > I_MIN
        stretched = round(255 * (roi_pix - I_MIN) / (I_MAX - I_MIN));
        I_str(roi_mask) = uint8(max(0,min(255,stretched)));
    end
end

% ==========================================================================
%  STAGE 3: EMBEDDING VIA HISTOGRAM SHIFTING (Sec. 2)
% ==========================================================================
function [I_emb_roi, P, Z, hs_dir, n_emb] = embed_hs(I_str, roi_mask, payload)
% Histogram shifting: find peak P and nearest zero bin Z in ROI.
% Shift bins between P and Z; embed bits at P.
    roi_pix = double(I_str(roi_mask));
    counts  = histcounts(roi_pix, 0:256);

    [~, pk_idx] = max(counts);
    P = pk_idx - 1;

    % Find nearest zero bin (prefer right)
    Z = -1;  hs_dir = 1;
    for z = P+1:255
        if counts(z+1) == 0, Z = z; hs_dir = 1; break; end
    end
    if Z < 0
        for z = P-1:-1:0
            if counts(z+1) == 0, Z = z; hs_dir = 0; break; end
        end
    end
    if Z < 0, I_emb_roi = I_str; P=0; Z=0; hs_dir=1; n_emb=0; return; end

    flat    = double(I_str(:));
    roi_idx = find(roi_mask(:));
    pay_ptr = 1;  n_emb = 0;

    for ii = 1:numel(roi_idx)
        idx = roi_idx(ii);  p = flat(idx);
        if hs_dir == 1      % right: shift (P,Z) right, embed at P
            if p > P && p < Z
                flat(idx) = p + 1;
            elseif p == P && pay_ptr <= numel(payload)
                flat(idx) = p + payload(pay_ptr);
                pay_ptr = pay_ptr + 1;
                n_emb   = n_emb + 1;
            end
        else                % left: shift (Z,P) left, embed at P
            if p < P && p > Z
                flat(idx) = p - 1;
            elseif p == P && pay_ptr <= numel(payload)
                flat(idx) = p - payload(pay_ptr);
                pay_ptr = pay_ptr + 1;
                n_emb   = n_emb + 1;
            end
        end
    end
    I_emb_roi = uint8(reshape(flat, size(I_str)));
end

% ==========================================================================
%  STAGE 4: NROI EMBEDDING (histogram shifting on NROI)
% ==========================================================================
function [I_out, P_N, Z_N, ndir] = embed_nroi(I, nroi_mask, payload)
    flat   = double(I(:));
    counts = histcounts(flat(nroi_mask), 0:256);
    [~, pk_idx] = max(counts);
    P_N = pk_idx - 1;
    Z_N = -1; ndir = 1;
    for z = P_N+1:255
        if counts(z+1)==0, Z_N=z; ndir=1; break; end
    end
    if Z_N < 0
        for z = P_N-1:-1:0
            if counts(z+1)==0, Z_N=z; ndir=0; break; end
        end
    end
    if Z_N < 0, I_out=I; P_N=0; Z_N=0; ndir=1; return; end
    flat_out = flat;
    nroi_idx = find(nroi_mask(:));
    pay_ptr  = 1;
    for ii = 1:numel(nroi_idx)
        idx=nroi_idx(ii); p=flat_out(idx);
        if ndir==1
            if p>P_N && p<Z_N, flat_out(idx)=p+1;
            elseif p==P_N && pay_ptr<=numel(payload)
                flat_out(idx)=p+payload(pay_ptr); pay_ptr=pay_ptr+1; end
        else
            if p<P_N && p>Z_N, flat_out(idx)=p-1;
            elseif p==P_N && pay_ptr<=numel(payload)
                flat_out(idx)=p-payload(pay_ptr); pay_ptr=pay_ptr+1; end
        end
    end
    I_out = uint8(reshape(flat_out, size(I)));
end

% ==========================================================================
%  MAIN EMBEDDING PIPELINE
% ==========================================================================
function [I_emb, meta] = rdhacem_embed(I, payload)
    roi_mask = get_roi_mask(I);
    [I_str, I_MIN, I_MAX] = auto_contrast_stretch(I, roi_mask);
    [I_emb_roi, P, Z, hs_dir, n_emb_roi] = embed_hs(I_str, roi_mask, payload);

    I_emb = I_emb_roi;
    P_N=0; Z_N=0; ndir=1;
    if numel(payload) > n_emb_roi
        [I_emb, P_N, Z_N, ndir] = embed_nroi(I_emb_roi, ~roi_mask, payload(n_emb_roi+1:end));
    end

    meta = struct('roi_mask',roi_mask, 'I_MIN',I_MIN, 'I_MAX',I_MAX, ...
        'P',P, 'Z',Z, 'hs_dir',hs_dir, 'n_emb_roi',n_emb_roi, ...
        'P_N',P_N, 'Z_N',Z_N, 'ndir',ndir);
end

% ==========================================================================
%  EXTRACTION AND RECOVERY
% ==========================================================================
function [I_rec, D_ext] = rdhacem_extract(I_emb, meta)
    roi_mask = meta.roi_mask;
    P=meta.P; Z=meta.Z; hs_dir=meta.hs_dir;
    I_MIN=meta.I_MIN; I_MAX=meta.I_MAX;
    n_emb_roi=meta.n_emb_roi;
    P_N=meta.P_N; Z_N=meta.Z_N; ndir=meta.ndir;

    % Step A: reverse NROI
    flat = double(I_emb(:));
    D_nroi = [];
    nroi_idx = find(~roi_mask(:));
    if Z_N ~= P_N
        for ii = 1:numel(nroi_idx)
            idx=nroi_idx(ii); p=flat(idx);
            if ndir==1
                if p==P_N+1, D_nroi(end+1)=1; flat(idx)=P_N;
                elseif p==P_N, D_nroi(end+1)=0;
                elseif p>P_N && p<=Z_N, flat(idx)=p-1; end
            else
                if p==P_N-1, D_nroi(end+1)=1; flat(idx)=P_N;
                elseif p==P_N, D_nroi(end+1)=0;
                elseif p<P_N && p>=Z_N, flat(idx)=p+1; end
            end
        end
    end
    I_step1 = uint8(reshape(flat, size(I_emb)));

    % Step B: extract from ROI + recover histogram shift
    flat2   = double(I_step1(:));
    D_roi   = zeros(1,n_emb_roi,'uint8');
    bit_ptr = 1;
    roi_idx = find(roi_mask(:));
    for ii = 1:numel(roi_idx)
        idx=roi_idx(ii); p=flat2(idx);
        if hs_dir==1
            if p==P+1, D_roi(bit_ptr)=1; flat2(idx)=P; bit_ptr=bit_ptr+1;
            elseif p==P, D_roi(bit_ptr)=0; bit_ptr=bit_ptr+1;
            elseif p>P && p<=Z, flat2(idx)=p-1; end
        else
            if p==P-1, D_roi(bit_ptr)=1; flat2(idx)=P; bit_ptr=bit_ptr+1;
            elseif p==P, D_roi(bit_ptr)=0; bit_ptr=bit_ptr+1;
            elseif p<P && p>=Z, flat2(idx)=p+1; end
        end
    end
    D_ext = [D_roi(1:bit_ptr-1), D_nroi];

    % Step C: inverse auto contrast stretch
    if I_MAX > I_MIN
        roi_pix = flat2(roi_mask(:));
        orig_pix = round((I_MAX-I_MIN)*roi_pix/255 + I_MIN);
        flat2(roi_mask(:)) = max(I_MIN,min(I_MAX,orig_pix));
    end
    I_rec = uint8(reshape(flat2, size(I_emb)));
end

% ==========================================================================
%  METRICS
% ==========================================================================
function p = compute_psnr(I, I_emb)
    mse = mean((double(I(:))-double(I_emb(:))).^2);
    if mse==0, p=Inf; else, p=10*log10(255^2/mse); end
end

% ==========================================================================
%  SYNTHETIC MEDICAL IMAGE GENERATOR
% ==========================================================================
function imgs = generate_medical_images()
    imgs = cell(4,1); sz=512;
    rng(1); I=uint8(ones(sz)*20); cx=sz/2; cy=sz/2;
    for r=1:sz; for c=1:sz
        d=sqrt((r-cx)^2+(c-cy)^2)/(sz*0.35);
        if d<1, I(r,c)=uint8(min(255,80+round(120*exp(-d*2))+randi(20))); end
    end; end
    imgs{1}=I;
    rng(2); I=imgs{1};
    for r=round(sz*0.35):round(sz*0.65); for c=round(sz*0.4):round(sz*0.6)
        d=sqrt((r-cx)^2+(c-cy)^2)/(sz*0.12);
        if d<1, I(r,c)=uint8(max(0,double(I(r,c))-round(60*exp(-d*2)))); end
    end; end
    imgs{2}=I;
    rng(3); I=uint8(zeros(sz));
    for r=1:sz; for c=1:sz
        I(r,c)=uint8(40+randi(20));
        if mod(c,round(sz/8))<round(sz/32), I(r,c)=uint8(min(255,double(I(r,c))+120+randi(30))); end
    end; end
    imgs{3}=I;
    rng(4); I=uint8(ones(sz)*80);
    for r=1:sz; for c=1:sz
        if abs(r-sz/2)<sz/6
            bf=max(0,1-abs(c-sz/2)/(sz*0.3));
            I(r,c)=uint8(min(255,80+round(150*bf)+randi(15)));
        end
    end; end
    imgs{4}=I;
end
