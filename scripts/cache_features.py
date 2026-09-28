"""CPU frozen-image embeddings and auditable RGB/texture summaries in patch order."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import time
import h5py
import numpy as np
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault('TORCH_HOME', str(ROOT / 'data/weights'))
import torch
from torchvision.models import resnet18, ResNet18_Weights

NAMES = ['red_mean','green_mean','blue_mean','red_std','green_std','blue_std',
         'gray_mean','gray_std','gray_gradient_x_abs','gray_gradient_y_abs','gray_entropy_32bins']

def morphology(images):
    x = images.astype(np.float32) / 255
    gray = x @ np.array([.299, .587, .114], dtype=np.float32)
    rows = []
    for rgb, g in zip(x, gray):
        hist = np.histogram(g, bins=32, range=(0,1))[0].astype(float)
        hist = hist[hist > 0] / g.size
        rows.append(np.r_[rgb.mean((0,1)), rgb.std((0,1)), g.mean(), g.std(),
                          np.abs(np.diff(g,axis=1)).mean(), np.abs(np.diff(g,axis=0)).mean(),
                          -(hist*np.log2(hist)).sum()])
    return np.asarray(rows, dtype=np.float32)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--threads', type=int, default=4)
    parser.add_argument('--batch-size', type=int, default=64)
    args = parser.parse_args()
    torch.set_num_threads(args.threads)
    torch.manual_seed(20260927)
    weights = ResNet18_Weights.IMAGENET1K_V1
    model = resnet18(weights=weights).eval()
    model.fc = torch.nn.Identity()
    for parameter in model.parameters():
        parameter.requires_grad_(False)
    transform = weights.transforms()
    checkpoint = ROOT / 'data/weights/hub/checkpoints/resnet18-f37072fd.pth'
    weight_sha = hashlib.file_digest(checkpoint.open('rb'), 'sha256').hexdigest()
    config = json.loads((ROOT/'configs/c0.json').read_text())
    for sample in sum(config['splits'].values(), []):
        target = ROOT/'data/cache'/f'{sample}_features.npz'
        if target.exists():
            print(f'{sample}: cached', flush=True)
            continue
        start = time.monotonic()
        embeddings, morph = [], []
        with h5py.File(ROOT/'data/raw/IDC/patches'/f'{sample}.h5') as f:
            barcodes = f['barcode'][:].reshape(-1).astype(str)
            for begin in range(0,len(barcodes),args.batch_size):
                images = f['img'][begin:begin+args.batch_size]
                batch = torch.stack([transform(Image.fromarray(im)) for im in images])
                with torch.inference_mode():
                    embeddings.append(model(batch).numpy())
                morph.append(morphology(images))
                if begin % (args.batch_size*20) == 0:
                    print(f'{sample}: {begin}/{len(barcodes)}', flush=True)
        np.savez_compressed(target, embeddings=np.concatenate(embeddings),
                            morphology=np.concatenate(morph), barcodes=barcodes,
                            morphology_names=np.array(NAMES))
        target.with_suffix('.json').write_text(json.dumps({
            'encoder':config['encoder'], 'weights_sha256':weight_sha,
            'torch_version':torch.__version__, 'threads':args.threads,
            'seconds':time.monotonic()-start, 'features':NAMES,
            'note':'RGB/texture statistics only; no validated cell types or nuclear segmentation.'}, indent=2)+'\n')
        print(f'{sample}: done in {time.monotonic()-start:.1f}s', flush=True)

if __name__ == '__main__':
    main()
