"""S3 parallel download for 16 pQTL files (4 SomaScan + 12 Olink).
Uses chunked streaming to avoid memory issues.
"""
import boto3, os, sys, time
from concurrent.futures import ThreadPoolExecutor, as_completed
from botocore.config import Config

ENDPOINT = "https://s3-ext.decode.is:10443"
ACCESS   = "7WOFN1QZ6J23BE0I6SV5"
SECRET   = "o0N6nczkjDtxZNKeMeAVkBNg3UH6PEhLj4cHX7GV"
BUCKET   = "largescaleplasma-2023"
DEST     = "D:/多组学MR数据/pqtl"
THREADS  = 3  # 3 concurrent ~900MB streams
os.makedirs(DEST, exist_ok=True)

files = [
    "final_somascan_smp/Proteomics_SMP_PC0_3434_34_FN1_FN1_3_10032022.txt.gz",
    "final_somascan_smp/Proteomics_SMP_PC0_5019_16_UCHL1_PGP9_5_10032022.txt.gz",
    "final_somascan_smp/Proteomics_SMP_PC0_2737_22_NOV_NovH_10032022.txt.gz",
    "final_somascan_smp/Proteomics_SMP_PC0_4920_10_LYZ_Lysozyme_10032022.txt.gz",
    "final_olink_ukb_bi/GBR_UKB_OLINK2_OID30787_FN1_Fibronectin_adjAgeSexPC_InvNorm_22122022.txt.gz",
    "final_olink_ukb_bi/GBR_UKB_OLINK_OID20299_CCN3_CCN_family_member_3_adjAgeSexBatPC_InvNorm_22122022.txt.gz",
    "final_olink_ukb_bi/GBR_UKB_OLINK2_OID31414_GABARAP_Gamma_aminobutyric_acid_receptor_associated_protein_adjAgeSexPC_InvNorm_22122022.txt.gz",
    "final_olink_ukb_bi/GBR_UKB_OLINK2_OID31253_GABARAPL1_Gamma_aminobutyric_acid_receptor_associated_protein_like_1_adjAgeSexPC_InvNorm_22122022.txt.gz",
    "final_olink_ukb_af/GBR_UKB_Africa_OLINK2_OID30787_FN1_Fibronectin_adjAgeSexPC_InvNorm_22122022.txt.gz",
    "final_olink_ukb_af/GBR_UKB_Africa_OLINK_OID20299_CCN3_CCN_family_member_3_adjAgeSexBatPC_InvNorm_22122022.txt.gz",
    "final_olink_ukb_af/GBR_UKB_Africa_OLINK2_OID31414_GABARAP_Gamma_aminobutyric_acid_receptor_associated_protein_adjAgeSexPC_InvNorm_22122022.txt.gz",
    "final_olink_ukb_af/GBR_UKB_Africa_OLINK2_OID31253_GABARAPL1_Gamma_aminobutyric_acid_receptor_associated_protein_like_1_adjAgeSexPC_InvNorm_22122022.txt.gz",
    "final_olink_ukb_sa/GBR_UKB_SAsia_OLINK2_OID30787_FN1_Fibronectin_adjAgeSexPC_InvNorm_22122022.txt.gz",
    "final_olink_ukb_sa/GBR_UKB_SAsia_OLINK_OID20299_CCN3_CCN_family_member_3_adjAgeSexBatPC_InvNorm_22122022.txt.gz",
    "final_olink_ukb_sa/GBR_UKB_SAsia_OLINK2_OID31414_GABARAP_Gamma_aminobutyric_acid_receptor_associated_protein_adjAgeSexPC_InvNorm_22122022.txt.gz",
    "final_olink_ukb_sa/GBR_UKB_SAsia_OLINK2_OID31253_GABARAPL1_Gamma_aminobutyric_acid_receptor_associated_protein_like_1_adjAgeSexPC_InvNorm_22122022.txt.gz",
]

todo = []
for key in files:
    local = os.path.join(DEST, key.split("/")[-1])
    if not os.path.exists(local):
        todo.append((key, local))

print(f"Total: {len(files)}, done: {len(files)-len(todo)}, remain: {len(todo)}",
      flush=True)
if not todo:
    print("All done!"); sys.exit(0)

s3 = boto3.client("s3",
    endpoint_url=ENDPOINT,
    aws_access_key_id=ACCESS,
    aws_secret_access_key=SECRET,
    config=Config(signature_version="s3v4", connect_timeout=30, read_timeout=3600),
    region_name="us-east-1")

def dl(item):
    key, out = item
    tmp = out + ".tmp"
    name = key.split("/")[-1]
    try:
        resp = s3.get_object(Bucket=BUCKET, Key=key)
        total = int(resp["ContentLength"])
        downloaded = 0
        with open(tmp, "wb") as f:
            for chunk in resp["Body"].iter_chunks(chunk_size=16*1024*1024):
                f.write(chunk)
                downloaded += len(chunk)
        os.rename(tmp, out)
        speed = total / max(time.time() - t0.get(name, time.time()), 1) / 1e6
        return f"OK  {name} ({total/1e6:.0f}MB)"
    except Exception as e:
        if os.path.exists(tmp): os.remove(tmp)
        return f"FAIL {name}: {e}"

t0 = {}
ok = fail = 0
start = time.time()
with ThreadPoolExecutor(max_workers=THREADS) as ex:
    futures = {ex.submit(dl, item): item for item in todo}
    for i, f in enumerate(as_completed(futures), 1):
        msg = f.result()
        print(f"[{i}/{len(todo)}] {msg}", flush=True)
        if msg.startswith("OK"): ok += 1
        else: fail += 1

elapsed = time.time() - start
print(f"\nDone. OK={ok} FAIL={fail}  Time={elapsed/60:.1f}min", flush=True)
