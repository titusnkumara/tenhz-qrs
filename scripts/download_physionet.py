#!/usr/bin/env python3
"""
export_physionet.py - export PhysioNet ECG databases for the 10 Hz QRS pipeline.

Per record it writes, in the format load_db.m expects:
    <db>-mat/<rec>.txt       first lead, one value per line
    <db>-mat/<rec>-ann.txt   beat sample indices (1-based by default)

Raw WFDB files are mirrored under <data>/physionet_raw/<db>/. Downloads skip complete
files and resume partial ones, outputs are written atomically, and a per-database
manifest records what finished, so re-running the same command continues from
where it stopped.

Usage:
    python scripts/download_physionet.py                       # all databases, 8 workers
    python scripts/download_physionet.py --dbs ltafdb --workers 8
    python scripts/download_physionet.py --dbs qtdb --overwrite   # force re-export

Output goes to <repo>/data (the folder above scripts/), or to $TENHZ_DATA if set,
or to --out-root. Point it at an existing data folder to reuse finished records.
"""
import argparse
import contextlib
import io
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

import subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)                               # this script lives in <ROOT>/scripts/
DATA = os.environ.get("TENHZ_DATA") or os.path.join(ROOT, "data")   # same rule as Q10.data_root


def _ensure_dependencies():
    """Install the packages in requirements.txt into this Python if any is missing."""
    try:
        import numpy  # noqa: F401
        import wfdb   # noqa: F401
    except ImportError:
        req = os.path.join(ROOT, "requirements.txt")
        print(f"Installing dependencies from {req} ...", flush=True)
        subprocess.check_call([sys.executable, "-m", "pip", "install", "-r", req])


_ensure_dependencies()
import numpy as np  # noqa: E402
import wfdb         # noqa: E402

# Standard WFDB beat symbols. Everything else (rhythm '+', noise '~', ST 's',
# T-wave 'T', comments '"', waveform boundaries '(' ')' 'p' 't', ...) is dropped.
BEAT_SYMBOLS = list("NLRBAaJSVrFejnE/fQ?")

# name: (PhysioNet dir, annotator, expected fs, annotation quality)
DATABASES = {
    "mitdb": ("mitdb", "atr", 360, "reference"),
    "nstdb": ("nstdb", "atr", 360, "reference"),
    "qtdb":  ("qtdb",  "atr", 250, "reference"),   # .atr only, never automatic .pu/.qrs
    "edb":   ("edb",   "atr", 250, "reference"),
    "nsrdb": ("nsrdb", "atr", 128, "reference"),
    "svdb":  ("svdb",  "atr", 128, "reference"),
    "twadb": ("twadb", "qrs", 500, "automatic"),   # no reviewed beat annotations
    "ltafdb": ("ltafdb", "atr", 128, "reference"), # held-out hard database (PocketECG + technician review)
}

# Databases whose non-beat annotations are also written, to <db>-meta/<rec>.csv
# (outside <db>-mat on purpose: the loaders treat every .txt in <db>-mat as a record).
# For LTAFDB these hold the rhythm labels ("(AFIB", "(N", ...) used to split AF from non-AF.
META_DBS = {"ltafdb"}

# Scheduling order by record length (longest first): nsrdb ~24 h, edb 2 h,
# mitdb/nstdb/svdb ~30 min, qtdb 15 min, twadb 2 min
SIZE_RANK = {"ltafdb": 0, "nsrdb": 0, "edb": 1, "mitdb": 2, "nstdb": 2, "svdb": 2, "qtdb": 3, "twadb": 4}


# ----------------------------------------------------------------------------
# Direct file fetching. wfdb.dl_database/get_record_list send several index and
# access-check requests per record (content/<db>/, RECORDS, ANNOTATORS, one per file);
# with parallel workers that is enough to get 502 Bad Gateway from PhysioNet. Here each
# record costs 3-4 requests (header, signal, annotation), paced and retried with long waits.
import random
import urllib.error
import urllib.request

BASE_URL = "https://physionet.org/files"
VERSION = "1.0.0"          # all databases used here are version 1.0.0 on PhysioNet
UA = "tenhz-qrs-export/1.0 (research use; contact: MARCS Institute, Western Sydney University)"
RETRY_CODES = {429, 500, 502, 503, 504}
_last_request = [0.0]


def fetch(url, dest, args, must_exist=True):
    """Download url to dest (atomic, resumable via dest+'.part'). Returns False if the
    server says 404 and must_exist is False; raises after args.net_retries failed attempts."""
    if os.path.exists(dest):
        return True
    part = dest + ".part"
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    last = None
    for k in range(args.net_retries + 1):
        wait = 0.0
        try:
            gap = args.delay - (time.time() - _last_request[0])
            if gap > 0:
                time.sleep(gap)
            _last_request[0] = time.time()
            have = os.path.getsize(part) if os.path.exists(part) else 0
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            if have:
                req.add_header("Range", f"bytes={have}-")
            with urllib.request.urlopen(req, timeout=60) as r:
                mode = "ab" if (have and r.status == 206) else "wb"   # 200: server ignored Range
                with open(part, mode) as f:
                    while True:
                        chunk = r.read(1 << 20)
                        if not chunk:
                            break
                        f.write(chunk)
            os.replace(part, dest)
            return True
        except urllib.error.HTTPError as e:
            if e.code == 404:
                if must_exist:
                    raise FileNotFoundError(f"404 {url}")
                return False
            if e.code == 416 and os.path.exists(part):       # .part is already complete or stale
                os.remove(part)
            last = f"HTTP {e.code} {url}"
            if e.code not in RETRY_CODES and e.code != 416:
                raise RuntimeError(last)
            ra = e.headers.get("Retry-After") if e.headers else None
            wait = float(ra) if ra and ra.isdigit() else 0.0
        except (urllib.error.URLError, TimeoutError, ConnectionError, OSError) as e:
            last = f"{type(e).__name__}: {e} {url}"
        if k < args.net_retries:
            time.sleep(max(wait, min(args.max_wait, args.backoff * 2 ** k)) * random.uniform(0.8, 1.2))
    raise RuntimeError(f"gave up after {args.net_retries + 1} attempts: {last}")


def signal_files(hea_path):
    """Signal file names listed in a WFDB header."""
    names = []
    with open(hea_path) as f:
        lines = [ln for ln in f.read().splitlines() if ln.strip() and not ln.startswith("#")]
    nsig = int(lines[0].split()[1])
    for ln in lines[1:1 + nsig]:
        if ln.split()[0] not in names:
            names.append(ln.split()[0])
    return names


def fetch_record(db, rec, ann_ext, args):
    """Fetch header, signal file(s) and annotation of one record into raw_root/<db>.
    Returns False if the annotation file does not exist on PhysioNet. A <rec>.fetched
    marker means everything was downloaded completely by this function; raw files without
    it (e.g. left by an interrupted wfdb run) are discarded and fetched again."""
    raw_dir = os.path.join(args.raw_root, db)
    marker = os.path.join(raw_dir, f"{rec}.fetched")
    base_url = f"{args.base_url}/{DATABASES[db][0]}/{VERSION}"
    if os.path.exists(marker) and os.path.exists(os.path.join(raw_dir, f"{rec}.hea")):
        return os.path.exists(os.path.join(raw_dir, f"{rec}.{ann_ext}"))
    os.makedirs(raw_dir, exist_ok=True)
    for f in os.listdir(raw_dir):                          # legacy / unverified files of this record
        if f.startswith(rec + ".") and not f.endswith(".part") and f != f"{rec}.fetched":
            os.remove(os.path.join(raw_dir, f))
    fetch(f"{base_url}/{rec}.hea", os.path.join(raw_dir, f"{rec}.hea"), args)
    for name in signal_files(os.path.join(raw_dir, f"{rec}.hea")):
        fetch(f"{base_url}/{name}", os.path.join(raw_dir, name), args)
    has_ann = fetch(f"{base_url}/{rec}.{ann_ext}", os.path.join(raw_dir, f"{rec}.{ann_ext}"), args,
                    must_exist=False)
    open(marker, "w").close()
    return has_ann


# ----------------------------------------------------------------------------
def get_records(db, args):
    path = os.path.join(args.raw_root, db, "RECORDS")
    try:
        fetch(f"{args.base_url}/{DATABASES[db][0]}/{VERSION}/RECORDS", path, args)
    except Exception as e:
        raise RuntimeError(f"{db}: could not fetch RECORDS ({e})")
    recs = [ln.strip() for ln in open(path).read().splitlines() if ln.strip()]
    if db == "nstdb":   # only the 118/119 noise-stress records, not bw/em/ma
        recs = [r for r in recs if r.startswith(("118e", "119e"))]
    return recs


def write_txt_atomic(path, arr, fmt):
    tmp = path[:-4] + ".tmp"
    np.savetxt(tmp, arr, fmt=fmt)
    os.replace(tmp, path)


def out_dir_for(db, args):
    return os.path.join(args.out_root, f"{db}-mat")


def manifest_path(db, args):
    return os.path.join(out_dir_for(db, args), "_manifest.json")


def load_manifest(db, args):
    p = manifest_path(db, args)
    if os.path.exists(p):
        with open(p) as f:
            return json.load(f)
    return {}


def save_manifest(db, args, manifest):
    p = manifest_path(db, args)
    tmp = p + ".tmp"
    with open(tmp, "w") as f:
        json.dump(manifest, f, indent=1, sort_keys=True)
    os.replace(tmp, p)


def is_done(db, rec, manifest, args):
    entry = manifest.get(rec)
    if not entry or entry.get("status") != "OK":
        return False
    d = out_dir_for(db, args)
    return (os.path.exists(os.path.join(d, f"{rec}.txt")) and
            os.path.exists(os.path.join(d, f"{rec}-ann.txt")))


# ----------------------------------------------------------------------------
def export_record(db, rec, args):
    pn_dir, ann_ext, fs_expected, quality = DATABASES[db]
    raw_dir = os.path.join(args.raw_root, db)
    out_dir = out_dir_for(db, args)
    os.makedirs(raw_dir, exist_ok=True)
    os.makedirs(out_dir, exist_ok=True)

    last_err = None
    for attempt in range(args.retries + 1):
        try:
            # Paced, resumable download of this record's header, signal and annotation
            if not fetch_record(db, rec, ann_ext, args):
                return {"record": rec, "status": "NO_ANNOTATION",
                        "error": f"no .{ann_ext} file on PhysioNet"}
            base = os.path.join(raw_dir, rec)

            hdr = wfdb.rdheader(base)
            try:
                sig, fields = wfdb.rdsamp(base, channels=[0])
                if len(sig) != hdr.sig_len:
                    raise ValueError(f"signal has {len(sig)} samples, header says {hdr.sig_len}")
            except Exception:
                for nm in signal_files(base + ".hea") + [rec + ".fetched"]:   # force a clean re-download
                    p = os.path.join(raw_dir, nm)
                    if os.path.exists(p):
                        os.remove(p)
                raise
            ecg = sig[:, 0].astype(float)
            n_nan = int(np.isnan(ecg).sum())
            ecg = np.nan_to_num(ecg, nan=0.0)          # MATLAB load can't parse 'nan'
            fs = float(fields["fs"])

            ann = wfdb.rdann(base, ann_ext)
            sym = np.asarray(ann.symbol, dtype=object)
            is_beat = np.isin(sym, BEAT_SYMBOLS)
            qrs = np.asarray(ann.sample, dtype=np.int64)[is_beat]
            inside = (qrs >= 0) & (qrs < len(ecg))
            n_outside = int((~inside).sum())
            qrs = qrs[inside]
            n_short_rr = int((np.diff(qrs) < 0.1 * fs).sum())   # reported, NOT removed

            warnings = []
            if fs != fs_expected:
                warnings.append(f"fs={fs:g}, expected {fs_expected}")
            if n_nan:
                warnings.append(f"{n_nan} NaN samples set to 0")
            if n_outside:
                warnings.append(f"{n_outside} annotations outside signal dropped")
            if n_short_rr:
                warnings.append(f"{n_short_rr} beat pairs < 100 ms apart (kept)")
            if len(qrs) == 0:
                warnings.append("no beat annotations")

            if args.one_based:
                qrs = qrs + 1
            write_txt_atomic(os.path.join(out_dir, f"{rec}.txt"), ecg, args.fmt)
            write_txt_atomic(os.path.join(out_dir, f"{rec}-ann.txt"), qrs, "%d")
            if db in META_DBS:
                meta_dir = os.path.join(args.out_root, f"{db}-meta")
                os.makedirs(meta_dir, exist_ok=True)
                aux = np.asarray(ann.aux_note, dtype=object)
                allsamp = np.asarray(ann.sample, dtype=np.int64)
                tmp = os.path.join(meta_dir, f"{rec}.tmp")
                with open(tmp, "w") as fh:
                    fh.write("sample_1based,symbol,aux_note\n")
                    for i in np.flatnonzero(~is_beat):
                        fh.write(f"{int(allsamp[i]) + 1},{sym[i]},{str(aux[i]).rstrip(chr(0))}\n")
                os.replace(tmp, os.path.join(meta_dir, f"{rec}.csv"))

            return {"record": rec, "status": "OK", "fs": fs,
                    "n_samples": int(len(ecg)), "hours": len(ecg) / fs / 3600,
                    "n_beats": int(len(qrs)), "n_nonbeat_dropped": int((~is_beat).sum()),
                    "lead": fields["sig_name"][0], "annotator": ann_ext,
                    "annotation_quality": quality, "one_based": args.one_based,
                    "warnings": warnings}
        except Exception as e:
            last_err = f"{type(e).__name__}: {e}"
            if attempt < args.retries:
                time.sleep(args.backoff * 2 ** attempt)
    return {"record": rec, "status": "ERROR", "error": last_err}


# ----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dbs", nargs="+", default=list(DATABASES), choices=list(DATABASES))
    ap.add_argument("--workers", type=int, default=6, help="parallel records (PhysioNet throttles heavy parallel use)")
    ap.add_argument("--retries", type=int, default=2, help="re-processing attempts per record")
    ap.add_argument("--net-retries", type=int, default=8, help="download attempts per file")
    ap.add_argument("--backoff", type=float, default=5.0, help="seconds, doubled per retry")
    ap.add_argument("--max-wait", type=float, default=120.0, help="longest wait between retries, seconds")
    ap.add_argument("--delay", type=float, default=0.5, help="minimum seconds between requests per worker")
    ap.add_argument("--base-url", default=BASE_URL, help="PhysioNet files URL (change only for a mirror)")
    ap.add_argument("--overwrite", action="store_true", help="re-export finished records")
    ap.add_argument("--out-root", default=DATA, help=f"where <db>-mat folders go (default: {DATA})")
    ap.add_argument("--raw-root", default=None, help="raw WFDB files (default: <out-root>/physionet_raw)")
    ap.add_argument("--fmt", default="%.6f", help="signal number format")
    ap.add_argument("--zero-based", dest="one_based", action="store_false")
    args = ap.parse_args()
    args.raw_root = args.raw_root or os.path.join(args.out_root, "physionet_raw")
    print(f"Exporting to {args.out_root} (raw files in {args.raw_root})")

    manifests, todo = {}, []
    for db in args.dbs:
        os.makedirs(out_dir_for(db, args), exist_ok=True)
        manifests[db] = load_manifest(db, args)
        try:
            recs = get_records(db, args)
        except RuntimeError as e:
            sys.exit(f"{e}\nPhysioNet did not answer. Run the same command again later: finished records are skipped.")
        pending = [r for r in recs if args.overwrite or not is_done(db, r, manifests[db], args)]
        print(f"{db:6s}: {len(recs):3d} records, {len(recs) - len(pending):3d} already done, "
              f"{len(pending):3d} to export")
        if DATABASES[db][3] != "reference":
            print(f"        note: {db} uses automatic .{DATABASES[db][1]} annotations, "
                  f"not reviewed reference beats")
        todo += [(db, r) for r in pending]

    # Longest records first, so the 24 h NSRDB files don't end up as a slow
    # single-worker tail after everything else has finished
    todo.sort(key=lambda t: SIZE_RANK.get(t[0], 9))

    print(f"\nExporting {len(todo)} records with {args.workers} worker processes...\n")
    t0 = time.time()
    pool = ProcessPoolExecutor(max_workers=args.workers)
    try:
            futures = {pool.submit(export_record, db, r, args): (db, r) for db, r in todo}
            for i, fut in enumerate(as_completed(futures), 1):
                db, rec = futures[fut]
                res = fut.result()
                manifests[db][rec] = res
                save_manifest(db, args, manifests[db])      # main thread only: no lock needed
                tag = f"[{i}/{len(todo)}] {db}/{rec}"
                if res["status"] == "OK":
                    w = ("  ! " + "; ".join(res["warnings"])) if res["warnings"] else ""
                    print(f"OK  {tag}: {res['hours']:.2f} h, {res['n_beats']} beats, "
                          f"lead {res['lead']}{w}")
                else:
                    # Remove files left by older exporters (e.g. QTDB ecgpuwave labels),
                    # so load_db.m can never pick up a record that did not export cleanly
                    removed = []
                    for suffix in (".txt", "-ann.txt"):
                        p = os.path.join(out_dir_for(db, args), f"{rec}{suffix}")
                        if os.path.exists(p):
                            os.remove(p)
                            removed.append(os.path.basename(p))
                    note = f" (removed stale {', '.join(removed)})" if removed else ""
                    print(f"{res['status']:3s} {tag}: {res.get('error')}{note}")
    except KeyboardInterrupt:
        # stop the workers right away: otherwise they linger after Ctrl-C and flood the terminal
        print("\nInterrupted: stopping workers. Run the same command again to resume.", flush=True)
        procs = list((getattr(pool, "_processes", None) or {}).values())   # before shutdown clears them
        pool.shutdown(wait=False, cancel_futures=True)
        for p in procs:
            p.terminate()
        os._exit(130)
    pool.shutdown()

    print(f"\nDone in {(time.time() - t0) / 60:.1f} min\n")
    failed = 0
    for db in args.dbs:
        m = manifests[db]
        ok = [v for v in m.values() if v.get("status") == "OK"]
        bad = sorted(k for k, v in m.items() if v.get("status") != "OK")
        failed += sum(1 for k in bad if m[k]["status"] == "ERROR")
        print(f"{db:6s}: {len(ok):3d} OK, {sum(v['hours'] for v in ok):7.1f} h, "
              f"{sum(v['n_beats'] for v in ok):8d} beats"
              + (f" | not exported: {', '.join(bad)}" if bad else ""))

    if failed:
        print("\nSome records failed. Re-run the same command to retry them.")
    print("MATLAB rebuilds its caches automatically when the exported records change.")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:           # Ctrl-C anywhere (also before the workers start)
        print("\nInterrupted. Run the same command again to resume.", flush=True)
        os._exit(130)
