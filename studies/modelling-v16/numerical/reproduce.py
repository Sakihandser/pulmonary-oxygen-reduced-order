from pathlib import Path
import hashlib,json,subprocess,sys
root=Path(__file__).resolve().parent
manifest=json.loads((root/'release_sha256.json').read_text(encoding='utf-8'))
for item in manifest['files']:
 p=root/item['path']
 if hashlib.sha256(p.read_bytes()).hexdigest()!=item['sha256']:raise RuntimeError('Hash mismatch: '+item['path'])
commands=[['11_diagnostic_reproduction/code/reproduce_diagnostic.py']]
for command in commands:subprocess.run([sys.executable,'-X','utf8',*command],cwd=root,check=True)
print('All archived numerical commands completed. Physiological validation remains absent.')
