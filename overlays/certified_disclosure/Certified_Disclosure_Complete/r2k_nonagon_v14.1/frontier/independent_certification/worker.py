"""Boundable standalone direct-tile worker."""
import json
import sys
from pathlib import Path
from direct_engine import tile

request=json.loads(Path(sys.argv[1]).read_text())
result=tile(request['function'],request['geometry'],request['settings'])
Path(sys.argv[2]).write_text(json.dumps(result,indent=2))
print(json.dumps({k:result.get(k) for k in ['status','count','reason','seconds']}),flush=True)
