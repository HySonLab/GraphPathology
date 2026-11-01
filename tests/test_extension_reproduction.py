from pathlib import Path
import pytest
from graphpathology.extension_reproduction import staging_root


def test_replay_separates_states_and_never_overwrites_published_records(tmp_path):
 root=tmp_path/'repo';(root/'src').mkdir(parents=True);(root/'configs').mkdir();(root/'data/journal_revision/hest_training').mkdir(parents=True);p=root/'checkpoints/journal_revision/hest_learned';p.mkdir(parents=True);(p/'journal_v2_feature_scaler.npz').write_bytes(b'unchanged source scaler');p=root/'experiments/journal_revision/transformer_egnn';p.mkdir(parents=True);(p/'PREFLIGHT.json').write_text('{}');(p.parent/'ENVIRONMENT.json').write_text('{}');(p/'fits').mkdir();record=p/'fits/published.json';record.write_bytes(b'published evidence')
 out=staging_root(root,root/'outputs/replay');assert (out/'src').is_symlink();assert (out/'checkpoints/journal_revision/hest_learned/journal_v2_feature_scaler.npz').read_bytes()==b'unchanged source scaler';assert not (out/'experiments/journal_revision/transformer_egnn/fits').exists();assert record.read_bytes()==b'published evidence'
 with pytest.raises(FileExistsError):staging_root(root,out)
 with pytest.raises(ValueError):staging_root(root,root/'manuscript/replay')
