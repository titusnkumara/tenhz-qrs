"""Numba kernels for MinEl and TiTen (the element-by-element loops).

Copied unchanged from tools/tenhz_fast.py. Compiled without fastmath, so each
float operation rounds exactly as in the reference code.
"""
import math

import numpy as np
from numba import njit


# ===================================================================== kernels
@njit(cache=True)
def _bisect_left(a, x, lo, hi):
    while lo < hi:
        mid = (lo + hi) // 2
        if a[mid] < x:
            lo = mid + 1
        else:
            hi = mid
    return lo


@njit(cache=True)
def _bisect_right(a, x, lo, hi):
    while lo < hi:
        mid = (lo + hi) // 2
        if x < a[mid]:
            hi = mid
        else:
            lo = mid + 1
    return lo


@njit(cache=True)
def _sort_small(v, n):
    for i in range(1, n):            # insertion sort, n <= ~30
        x = v[i]
        j = i - 1
        while j >= 0 and v[j] > x:
            v[j + 1] = v[j]
            j -= 1
        v[j + 1] = x


@njit(cache=True)
def _median_small(src, s, e, buf):
    """core.median_small / np.median of src[s:e] (float64)."""
    n = e - s
    for i in range(n):
        buf[i] = src[s + i]
    _sort_small(buf, n)
    if n % 2:
        return buf[(n + 1) // 2 - 1]
    return (buf[n // 2 - 1] + buf[n // 2]) / 2


# ---------------------------------------------------------- 10 Hz stream
@njit(cache=True)
def _k_minmax_stream(x, edges):
    """Block min/max (first occurrence) merged in sample order, as
    core.minmax_stream (min first where one sample is both)."""
    nb = len(edges) - 1
    val = np.empty(2 * nb)
    idx = np.empty(2 * nb, np.int64)
    m = 0
    for k in range(nb):
        s = edges[k]
        e = edges[k + 1]
        lo = x[s]; hi = x[s]; il = s; ih = s
        for i in range(s + 1, e):
            v = x[i]
            if v < lo:
                lo = v; il = i
            elif v > hi:
                hi = v; ih = i
        if il < ih:
            idx[m] = il; val[m] = lo; idx[m + 1] = ih; val[m + 1] = hi; m += 2
        elif ih < il:
            idx[m] = ih; val[m] = hi; idx[m + 1] = il; val[m + 1] = lo; m += 2
        else:
            idx[m] = il; val[m] = lo; m += 1
    return val[:m], idx[:m]


# ------------------------------------------------------------------- MinEl
@njit(cache=True)
def _k_region_max(xs, st, en):
    nr = len(st)
    m = np.empty(nr)
    cl = np.empty(nr, np.int64)
    for r in range(nr):
        s = st[r]
        best = xs[s]; pos = s
        for i in range(s + 1, en[r] + 1):
            if xs[i] > best:
                best = xs[i]; pos = i
        m[r] = best; cl[r] = pos
    return m, cl


@njit(cache=True)
def _k_amp_ref(x, cl, half, perc):
    n = len(x)
    ref = np.empty(len(cl))
    buf = np.empty(2 * half + 1)
    for i in range(len(cl)):
        c = cl[i]
        a = max(0, c - half)
        b = min(n - 1, c + half)
        L = b - a + 1
        for q in range(L):
            buf[q] = x[a + q]
        _sort_small(buf, L)
        k = min(L, max(1, math.ceil(perc / 100 * L)))
        ref[i] = buf[k - 1]
    return ref


@njit(cache=True)
def _k_minel_main(cl, ca, base, tau, refr, rr_tol, low_f, high_f):
    nc = len(cl)
    locs = np.empty(nc, np.int64)
    amps = np.empty(nc)
    nl = 0
    rr = np.empty(10)
    nrr = 0
    buf = np.empty(10)
    for q in range(nc):
        t = cl[q]
        if nl >= 2 and nrr >= 1:
            mrr = _median_small(rr, 0, nrr, buf)
            dev = abs(tau[t] - (tau[locs[nl - 1]] + mrr)) / mrr
        else:
            dev = np.nan
        mult = low_f if (math.isnan(dev) or dev <= rr_tol) else high_f
        if nl > 0 and (tau[t] - tau[locs[nl - 1]]) < refr:
            continue
        if ca[q] > mult * base[q]:
            locs[nl] = t; amps[nl] = ca[q]; nl += 1
            if nl >= 2:
                r = tau[locs[nl - 1]] - tau[locs[nl - 2]]
                if nrr < 10:
                    rr[nrr] = r; nrr += 1
                else:
                    for z in range(9):
                        rr[z] = rr[z + 1]
                    rr[9] = r
    return locs[:nl], amps[:nl]


@njit(cache=True)
def _sb_flag(T, k, mult, buf):
    """core._sb_flag_at for 1-based beat k."""
    if k < 3:
        return False
    a = max(0, k - 12)
    nd = (k - 1) - a - 1                      # len(np.diff(T[a:k-1]))
    s = max(0, nd - 10)
    m = 0
    for q in range(s, nd):
        buf[m] = T[a + q + 1] - T[a + q]
        m += 1
    _sort_small(buf, m)
    if m % 2:
        med = buf[(m + 1) // 2 - 1]
    else:
        med = (buf[m // 2 - 1] + buf[m // 2]) / 2
    return T[k - 1] - T[k - 2] > mult * med


@njit(cache=True)
def _k_minel_searchback(locs, amps, tau, cl_all, ca_all, base_all, rel, tc, refr, sb_mult, sb_frac):
    nT = len(locs)
    cap = nT + len(cl_all) + 1
    L = np.empty(cap, np.int64); A = np.empty(cap); T = np.empty(cap); F = np.zeros(cap, np.bool_)
    buf = np.empty(12)
    for i in range(nT):
        L[i] = locs[i]; A[i] = amps[i]; T[i] = tau[locs[i]]
    for k in range(3, nT + 1):
        F[k - 1] = _sb_flag(T, k, sb_mult, buf)
    nc = len(tc)
    k = 3
    while True:
        q = k - 1
        while q < nT and not F[q]:
            q += 1
        if q >= nT:
            break
        k = q + 1
        lo = _bisect_left(tc, T[k - 2] + refr, 0, nc)
        hi = _bisect_right(tc, T[k - 1] - refr, 0, nc) - 1
        inserted = False
        if lo <= hi:
            j = lo
            bv = rel[lo]
            for z in range(lo + 1, hi + 1):
                if rel[z] > bv:
                    bv = rel[z]; j = z
            if ca_all[j] > sb_frac * base_all[j]:
                for z in range(nT, k - 1, -1):
                    L[z] = L[z - 1]; A[z] = A[z - 1]; T[z] = T[z - 1]; F[z] = F[z - 1]
                L[k - 1] = cl_all[j]; A[k - 1] = ca_all[j]; T[k - 1] = tc[j]; F[k - 1] = False
                nT += 1
                for q2 in range(k, min(k + 11, nT) + 1):
                    F[q2 - 1] = _sb_flag(T, q2, sb_mult, buf)
                inserted = True
        if not inserted:
            k += 1
        if k > nT:
            break
    return L[:nT].copy(), A[:nT].copy()


# ------------------------------------------------------------------- TiTen
@njit(cache=True)
def _k_peaks_time(n, locs, order, T, refr, w):
    blocked = np.zeros(n, np.uint8)
    keep = np.zeros(len(locs), np.bool_)
    for kk in range(len(order)):
        k = order[kk]
        L = locs[k]
        if blocked[L]:
            continue
        keep[k] = True
        lo = max(0, L - w); hi = min(n - 1, L + w)
        r0 = _bisect_right(T, T[L] - refr, lo, hi + 1)
        r1 = _bisect_left(T, T[L] + refr, lo, hi + 1)
        for r in range(r0, r1):
            blocked[r] = 1
    return keep


@njit(cache=True)
def _k_peaks_dist(n, locs, order, d):
    blocked = np.zeros(n, np.uint8)
    keep = np.zeros(len(locs), np.bool_)
    for kk in range(len(order)):
        k = order[kk]
        L = locs[k]
        if not blocked[L]:
            keep[k] = True
            for r in range(max(0, L - d + 1), min(n, L + d)):
                blocked[r] = 1
    return keep


@njit(cache=True)
def _raw_peak(E, n, k0):
    kk = min(k0, n - 1)
    l0 = max(0, kk - 2)
    m = E[l0]; pos = l0
    for i in range(l0 + 1, kk + 1):
        if E[i] > m:
            m = E[i]; pos = i
    return m, pos


@njit(cache=True)
def _mround1(x):
    return np.sign(x) * np.floor(np.abs(x) + 0.5)


@njit(cache=True)
def _k_titen_main(pk, pl, I, E, TAU, use_time, sig_alpha, rr_start, rr_low, rr_high, sb_mult,
                  thr_frac, thr_sig, thr_noise, thr_sig_r, thr_noise_r, wd):
    n = len(E); ni = len(I)
    t_fire = -np.inf; t_last = -np.inf     # watchdog: last time it fired / last reference time
    wmax = 0.0; wmax_r = 0.0               # largest candidate peaks since the last accepted beat or watchdog firing
    sig_memory = 1 - sig_alpha
    noise_alpha = 0.125; sb_alpha = 0.25
    noise_memory = 1 - noise_alpha; sb_memory = 1 - sb_alpha
    sig_lev = thr_sig; noise_lev = thr_noise
    sig_lev_r = thr_sig_r; noise_lev_r = thr_noise_r
    cap = 2 * len(pk) + 2
    qi = np.empty(cap, np.int64); nq = 0
    qr = np.empty(cap, np.int64); qa = np.empty(cap); nr = 0
    sel_rr = 0.0
    sb_from = -1; sb_to = -1; sb_v = -np.inf; sb_k = -1
    for q in range(len(pk)):
        pki = pk[q]; pli = pl[q]
        if wd > 0 and use_time:
            # Recovery from a deadlock. If no beat has been accepted for `wd` s, the signal level cannot be higher than the
            # largest candidate peak seen in that time: clamp it (e.g. a large initial transient otherwise keeps the
            # threshold above every real beat for the rest of the record). It never fires while beats are being found.
            t_acc = TAU[qi[nq - 1]] if nq > 0 else TAU[0]
            t_ref = max(t_acc, t_fire)
            if t_ref != t_last:
                wmax = 0.0; wmax_r = 0.0; t_last = t_ref
            if TAU[pli] - t_ref >= wd:
                if wmax > 0 and sig_lev > wmax:
                    sig_lev = wmax
                if wmax_r > 0 and sig_lev_r > wmax_r:
                    sig_lev_r = wmax_r
                t_fire = TAU[pli]
                thr_sig = noise_lev + thr_frac * (sig_lev - noise_lev); thr_noise = 0.5 * thr_sig
                thr_sig_r = noise_lev_r + thr_frac * (sig_lev_r - noise_lev_r); thr_noise_r = 0.5 * thr_sig_r
        bc = nq
        if bc >= rr_start:
            s = max(1, bc - 8)
            acc = 0.0
            for j in range(s, bc):
                acc += TAU[qi[j]] - TAU[qi[j - 1]]
            mrr = acc / (bc - s)
            cur = TAU[qi[bc - 1]] - TAU[qi[bc - 2]]
            if cur <= rr_low * mrr or cur >= rr_high * mrr:
                thr_sig *= 0.5; thr_sig_r *= 0.5
            else:
                sel_rr = mrr
        gsb = sb_mult * sel_rr
        if not use_time:
            gsb = _mround1(gsb)
        if sel_rr > 0 and bc > 0 and (TAU[pli] - TAU[qi[nq - 1]]) >= gsb:
            a0 = qi[nq - 1] + 2; b0 = pli - 2
            if a0 < b0 and a0 >= 0 and b0 <= ni - 1:
                if sb_from != a0 or b0 < sb_to:
                    sb_from = a0; sb_to = a0 - 1; sb_v = -np.inf; sb_k = -1
                for j in range(sb_to + 1, b0 + 1):
                    if I[j] > sb_v:
                        sb_v = I[j]; sb_k = j
                sb_to = b0
                pv = sb_v; k0 = sb_k
                if pv > thr_noise:
                    qi[nq] = k0; nq += 1
                    ar, kr = _raw_peak(E, n, k0)
                    if ar > thr_noise_r:
                        qr[nr] = kr; qa[nr] = ar; nr += 1
                        sig_lev_r = sb_alpha * ar + sb_memory * sig_lev_r
                    sig_lev = sb_alpha * pv + sb_memory * sig_lev
        ar, kr = _raw_peak(E, n, pli)
        if pki >= thr_sig:
            skip = False
            if nq >= 3:
                bw = 1
                for k in range(pli - 1, max(0, pli - 10) - 1, -1):
                    if I[k] < 0.5 * thr_sig:
                        break
                    bw += 1
                fw = 1
                for k in range(pli + 1, min(ni - 1, pli + 10) + 1):
                    if I[k] < 0.5 * thr_sig:
                        break
                    fw += 1
                if bw > 3 and fw > 3:
                    skip = True
                    noise_lev_r = noise_alpha * ar + noise_memory * noise_lev_r
                    noise_lev = noise_alpha * pki + noise_memory * noise_lev
            if not skip:
                qi[nq] = pli; nq += 1
                if ar >= thr_sig_r:
                    qr[nr] = kr; qa[nr] = ar; nr += 1
                    sig_lev_r = sig_alpha * ar + sig_memory * sig_lev_r
                sig_lev = sig_alpha * pki + sig_memory * sig_lev
        else:
            noise_lev_r = noise_alpha * ar + noise_memory * noise_lev_r
            noise_lev = noise_alpha * pki + noise_memory * noise_lev
        if noise_lev != 0 or sig_lev != 0:
            thr_sig = noise_lev + thr_frac * (sig_lev - noise_lev)
            thr_noise = 0.5 * thr_sig
        if noise_lev_r != 0 or sig_lev_r != 0:
            thr_sig_r = noise_lev_r + thr_frac * (sig_lev_r - noise_lev_r)
            thr_noise_r = 0.5 * thr_sig_r
        if wd > 0:
            if pki > wmax:
                wmax = pki
            if ar > wmax_r:
                wmax_r = ar
    return qr[:nr].copy(), qa[:nr].copy()


# ------------------------------------------------------------------ fusion
@njit(cache=True)
def _k_pair_flags(I, J, order, na, nb):
    pa = np.zeros(na, np.bool_); pb = np.zeros(nb, np.bool_)
    for z in range(len(order)):
        k = order[z]
        if not pa[I[k]] and not pb[J[k]]:
            pa[I[k]] = True; pb[J[k]] = True
    return pa, pb


@njit(cache=True)
def _k_local_median(rr, k):
    med = np.empty(len(k))
    buf = np.empty(11)
    for z in range(len(k)):
        i = k[z]
        med[z] = _median_small(rr, max(0, i - 6), min(len(rr), i + 5), buf)
    return med


@njit(cache=True)
def _k_viterbi(c, sc, logv, lo, hi, lam, start_pen):
    N = len(c)
    best = np.full(N, -np.inf)
    prev = np.full(N, -1, np.int64)
    far = 0.0; far_i = -1; ptr = 0; e = 0
    for i in range(N):
        while ptr < i and c[i] - c[ptr] > hi:
            if best[ptr] > far:
                far = best[ptr]; far_i = ptr
            ptr += 1
        b = far - start_pen; pj = far_i
        for j in range(ptr, i):
            d = c[i] - c[j]
            if d < lo:
                e += 1
                continue
            v = best[j] - lam * abs(logv[e])
            e += 1
            if v > b:
                b = v; pj = j
        best[i] = sc[i] + b; prev[i] = pj
    i = 0
    for z in range(1, N):            # np.argmax: first maximum
        if best[z] > best[i]:
            i = z
    out = np.empty(N, np.int64); m = 0
    while i >= 0:
        out[m] = c[i]; m += 1; i = prev[i]
    return np.sort(out[:m])
