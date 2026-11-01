"""Fresh ignored training replay directories; never replace published fit records/states."""
from pathlib import Path
import shutil


def staging_root(root,output):
 root=Path(root).resolve();out=Path(output).resolve()
 if not out.is_relative_to(root/'outputs'):raise ValueError('Training replay must use a new repository outputs/ directory')
 if out.exists():raise FileExistsError('Use a new replay directory; published or prior replay fits are never overwritten')
 inputs=[root/'data/journal_revision/hest_training',root/'checkpoints/journal_revision/hest_learned/journal_v2_feature_scaler.npz']
 for p in inputs:
  if not p.exists():raise FileNotFoundError('Restore/prepare verified source input before training replay: '+str(p))
 out.mkdir(parents=True)
 for name in ['src','configs']:(out/name).symlink_to(root/name,target_is_directory=True)
 target=out/'data/journal_revision';target.mkdir(parents=True);(target/'hest_training').symlink_to(inputs[0],target_is_directory=True)
 target=out/'checkpoints/journal_revision/hest_learned';target.mkdir(parents=True);(target/'journal_v2_feature_scaler.npz').symlink_to(inputs[1])
 target=out/'experiments/journal_revision/transformer_egnn';target.mkdir(parents=True)
 shutil.copy2(root/'experiments/journal_revision/transformer_egnn/PREFLIGHT.json',target/'PREFLIGHT.json')
 shutil.copy2(root/'experiments/journal_revision/ENVIRONMENT.json',target.parent/'ENVIRONMENT.json')
 adjustment=root/'experiments/journal_revision/transformer_egnn/RESOURCE_ADJUSTMENT.json'
 if adjustment.exists():shutil.copy2(adjustment,target/adjustment.name)
 return out
