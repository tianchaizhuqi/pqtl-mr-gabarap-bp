"""Download 1000G VCFs with 2 parallel workers + gzip verify."""
import requests, os, subprocess, time
from concurrent.futures import ThreadPoolExecutor, as_completed

BASE = "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502"
FMT  = "ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
DEST = "D:/多组学MR数据/ld_reference"
CHRS = ["8", "12", "14"]

def verify(path):
    r = subprocess.run(["gzip", "-t", path], capture_output=True)
    return r.returncode == 0

def dl(ch):
    out = os.path.join(DEST, FMT.format(ch))
    if os.path.exists(out) and verify(out):
        return f"chr{ch}: EXISTS OK ({os.path.getsize(out)/1e6:.0f}MB)"
    if os.path.exists(out): os.remove(out)

    url = f"{BASE}/{FMT.format(ch)}"
    for attempt in range(3):
        try:
            t0 = time.time()
            r = requests.get(url, stream=True, timeout=(30, 7200))
            r.raise_for_status()
            with open(out, "wb") as f:
                for chunk in r.iter_content(chunk_size=4*1024*1024):
                    f.write(chunk)
            elapsed = time.time() - t0
            sz = os.path.getsize(out)/1e6
            if verify(out):
                return f"chr{ch}: OK ({sz:.0f}MB, {elapsed/60:.1f}min)"
            os.remove(out)
            return f"chr{ch}: CORRUPT attempt {attempt+1}"
        except Exception as e:
            if os.path.exists(out): os.remove(out)
            if attempt < 2: time.sleep(5)
            else: return f"chr{ch}: FAIL {e}"

with ThreadPoolExecutor(max_workers=4) as ex:
    futures = {ex.submit(dl, c): c for c in CHRS}
    for f in as_completed(futures):
        print(f.result(), flush=True)

print("\n=== Verify ===")
for ch in CHRS:
    out = os.path.join(DEST, FMT.format(ch))
    ok = "OK" if verify(out) else "CORRUPT"
    print(f"chr{ch}: {ok}")
