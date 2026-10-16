"""Optional public-source download. Analysis itself is offline."""
import argparse, hashlib, json, urllib.request
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--output-dir',type=Path,required=True);a=p.parse_args()
a.output_dir.mkdir(parents=True,exist_ok=True)
root=Path(__file__).resolve().parents[1]
for record in json.loads((root/'data/source_manifest.json').read_text()):
    path=a.output_dir/(record['id']+'.pdf')
    if path.exists():data=path.read_bytes()
    else:
        request=urllib.request.Request(record['url'],headers={'User-Agent':'FalloutExtractionHistory/1.0 (academic source verification)'})
        with urllib.request.urlopen(request,timeout=90) as response:data=response.read()
    digest=hashlib.sha256(data).hexdigest()
    if digest!=record['sha256']:raise ValueError(f"Source changed or wrong content: {record['id']}; expected {record['sha256']}, found {digest}")
    if not path.exists():path.write_bytes(data)
    print(record['id'],digest)
