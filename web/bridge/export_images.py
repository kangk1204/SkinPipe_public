"""Coordinate-derived thumbnails. Generated conformers never represent docking evidence."""
from __future__ import annotations
from io import BytesIO
import hashlib
import math

IMAGE_POLICY = 'xyz-projection-v1-etkdgv3-seed20261007'


def thumbnail(smiles, asset=None):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from PIL import Image, ImageDraw
    import numpy as np

    mol = None
    kind, source = 'unavailable', None
    if asset:
        source = asset.get('projection_sha256') or asset.get('source_sha256')
        if asset.get('format') == 'sdf':
            mol = Chem.MolFromMolBlock(asset['data'].split('$$$$')[0], sanitize=False, removeHs=False)
        elif asset.get('format') == 'pdb':
            # Raw XYZ projection: do not invent covalent connectivity.
            mol = Chem.MolFromPDBBlock(asset['data'], sanitize=False, removeHs=False, proximityBonding=False)
        if mol is not None and mol.GetNumConformers() and mol.GetConformer().Is3D():
            kind = 'stored_docking_pose'
        else:
            mol = None
    if mol is None and smiles:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None and 0 < mol.GetNumAtoms() <= 180:
            mol = Chem.AddHs(mol)
            params = AllChem.ETKDGv3()
            params.randomSeed = 20261007
            params.maxIterations = 200
            # Supported RDKit versions bound native embedding, too.
            if hasattr(params, 'timeout'):
                params.timeout = 2
            if AllChem.EmbedMolecule(mol, params) == 0:
                kind = 'generated_conformer_not_docked'
                source = hashlib.sha256(smiles.encode()).hexdigest()
            else:
                mol = None
        else:
            mol = None
    meta = {'kind': kind, 'source_sha256': source, 'image_policy': IMAGE_POLICY}
    if mol is None or not mol.GetNumConformers():
        return None, meta
    conf = mol.GetConformer()
    atoms = [a for a in mol.GetAtoms() if a.GetAtomicNum() != 1]
    if not atoms:
        return None, meta
    xyz = np.array([list(conf.GetAtomPosition(a.GetIdx())) for a in atoms])
    if not np.isfinite(xyz).all():
        raise ValueError('Non-finite structure coordinates')
    xyz -= xyz.mean(axis=0)
    # Camera rotation affects the illustration only; source coordinates are untouched.
    a, b = math.radians(28), math.radians(32)
    rot = np.array([[math.cos(b), 0, math.sin(b)], [math.sin(a)*math.sin(b), math.cos(a), -math.sin(a)*math.cos(b)],
                    [-math.cos(a)*math.sin(b), math.sin(a), math.cos(a)*math.cos(b)]])
    xyz = xyz @ rot.T
    scale = min(142/max(float(np.ptp(xyz[:, 0])), 1), 70/max(float(np.ptp(xyz[:, 1])), 1))
    positions = {atom.GetIdx(): (float(p[0]*scale+90), float(-p[1]*scale+49), float(p[2]))
                 for atom, p in zip(atoms, xyz)}
    image = Image.new('RGB', (180, 100), 'white')
    draw = ImageDraw.Draw(image)
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if i in positions and j in positions:
            draw.line([positions[i][:2], positions[j][:2]], fill='#7b8790', width=3)
    colors = {6:'#566571', 7:'#4266dc', 8:'#e05252', 9:'#67b882', 15:'#dfaa52', 16:'#d4bd39', 17:'#5cb970', 35:'#a8664d', 53:'#a568bc'}
    for atom in sorted(atoms, key=lambda atom: positions[atom.GetIdx()][2]):
        x, y, _ = positions[atom.GetIdx()]
        r = 4 if atom.GetAtomicNum() == 6 else 5
        draw.ellipse((x-r, y-r, x+r, y+r), fill=colors.get(atom.GetAtomicNum(), '#93aab5'), outline='white')
    output = BytesIO()
    image.save(output, format='PNG', optimize=True)
    payload = output.getvalue()
    meta['image_sha256'] = hashlib.sha256(payload).hexdigest()
    return payload, meta
