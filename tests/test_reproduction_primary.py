"""Replacement identities cannot silently stand in for original publication states."""
import json
from pathlib import Path

import numpy as np
import pytest
import torch

from graphpathology.assets import sha
from graphpathology.core import mlp_model, train_mlp
from graphpathology.reproduction_primary import (
    ORIGIN, checked_namespace, config, image_predictions, load_heads,
)


def replacement_fixture(tmp_path):
    folder = tmp_path / 'replacement'; folder.mkdir()
    np.savez(folder / 'feature.npz', mean=np.zeros(512), scale=np.ones(512))
    np.savez(folder / 'target.npz', mean=np.full(50, 5.), scale=np.full(50, 3.))
    manifest = {'origin': ORIGIN, 'genes': [f'g{i}' for i in range(50)],
                'config': {'model': {'hidden': 256, 'dropout': .1}}, 'heads': [],
                'scalers': {n: sha(folder / n) for n in ['feature.npz', 'target.npz']}}
    for seed in range(3):
        net = mlp_model(512, 50, manifest['config']['model'])
        with torch.no_grad():
            for p in net.parameters(): p.zero_()
            list(net.modules())[-1].bias.fill_(seed)
        name = f'seed{seed}.pt'
        torch.save({'origin': ORIGIN, 'seed': seed, 'genes': manifest['genes'],
                    'state_dict': net.state_dict()}, folder / name)
        manifest['heads'].append({'seed': seed, 'file': name, 'sha256': sha(folder / name)})
    (folder / 'MANIFEST.json').write_text(json.dumps(manifest))
    return folder


def test_restored_ensemble_applies_head_scaler_and_casts_last(tmp_path):
    folder = replacement_fixture(tmp_path)
    output, individual = image_predictions(folder, np.zeros((3, 512), np.float32))
    assert output.dtype == np.float32
    np.testing.assert_array_equal(output, np.full((3, 50), 8., np.float32))
    np.testing.assert_array_equal(individual[0], np.full((3, 50), 5.))
    np.testing.assert_array_equal(individual[2], np.full((3, 50), 11.))


def test_modified_checkpoint_is_rejected_before_loading(tmp_path):
    folder = replacement_fixture(tmp_path)
    with (folder / 'seed0.pt').open('ab') as f: f.write(b'corruption')
    with pytest.raises(ValueError, match='checksum'): load_heads(folder)


def test_original_identity_cannot_be_accepted_by_replacement_loader(tmp_path):
    folder = replacement_fixture(tmp_path)
    p = folder / 'MANIFEST.json'; r = json.loads(p.read_text()); r['origin'] = 'ORIGINAL'
    p.write_text(json.dumps(r))
    with pytest.raises(ValueError, match='identity'): load_heads(folder)


@pytest.mark.parametrize('destination', ['manuscript/new', 'checkpoints/legacy/new', 'reproduction/runs'])
def test_namespace_rejects_scientific_outputs_and_root(tmp_path, destination):
    with pytest.raises(ValueError): checked_namespace(tmp_path, tmp_path / destination, 'reproduction/runs')


def test_namespace_rejects_symlink_escape(tmp_path):
    allowed = tmp_path / 'reproduction/runs'; allowed.mkdir(parents=True)
    other = tmp_path / 'manuscript'; other.mkdir()
    (allowed / 'alias').symlink_to(other, target_is_directory=True)
    with pytest.raises(ValueError): checked_namespace(tmp_path, allowed / 'alias/new', 'reproduction/runs')


def test_external_section_cannot_enter_training_split(tmp_path):
    root = Path(__file__).resolve().parents[1]
    cfg = json.loads((root / 'configs/reproduction/primary-v1.json').read_text())
    cfg['train'].append('TENX99'); p = tmp_path / 'config.json'; p.write_text(json.dumps(cfg))
    with pytest.raises(ValueError, match='split'): config(root, p)


def test_training_save_load_retains_selected_validation_state(tmp_path):
    torch.set_num_threads(1)
    rng = np.random.default_rng(12)
    x = rng.normal(size=(12, 4)).astype(np.float32); y = x[:, :2].copy()
    training = {'learning_rate': .01, 'weight_decay': .0001, 'batch_size': 4,
                'max_epochs': 3, 'patience': 2, 'min_delta': .0001}
    net, fit = train_mlp(x, y, x[:4], y[:4], {'hidden': 5, 'dropout': .1}, training, 0)
    with torch.inference_mode(): expected = net(torch.from_numpy(x)).numpy()
    path = tmp_path / 'selected.pt'; torch.save(net.state_dict(), path)
    restored = mlp_model(4, 2, {'hidden': 5, 'dropout': .1})
    restored.load_state_dict(torch.load(path, weights_only=True)); restored.eval()
    with torch.inference_mode(): np.testing.assert_array_equal(restored(torch.from_numpy(x)).numpy(), expected)
    assert fit['best_epoch'] <= fit['epochs_run'] <= 3


def test_reused_matched_state_keeps_original_export_path(tmp_path):
    from graphpathology.reproduction_artifacts import matched_reference
    path = tmp_path / 'experiments/journal_revision/transformer_egnn/evaluation/P/frozen_old/PANEL.csv.gz'
    path.parent.mkdir(parents=True); path.write_bytes(b'reference')
    assert matched_reference(tmp_path, 'P', 'old') == path
    assert matched_reference(tmp_path, 'P', 'unselected_validation_only_trial') is None


def test_strict_table_replay_does_not_hide_changed_support_or_failed_metrics():
    import pandas as pd
    from graphpathology.reproduction_replay import compare_tables
    old = pd.DataFrame({'id':['a','b'],'value':[1.,2.]})
    new = old.copy(); new.loc[1,'value'] += .01
    assert compare_tables(new,old,['id'],['value'])['status'] == 'NOT_REPRODUCED'
    with pytest.raises(ValueError,match='support'):
        compare_tables(old.iloc[:1],old,['id'],['value'])


def test_external_alignment_preserves_float64_instead_of_source_float32_contract():
    from graphpathology.reproduction_replay import validate_external_alignment
    y = np.array([[1.123456789012345]],dtype=np.float64)
    cached = {'targets':y,'barcodes':np.array(['a']),'genes':np.array(['g']),'xy_um':np.zeros((1,2))}
    real = {'targets':y.copy(),'barcodes':cached['barcodes'],'genes':cached['genes'],'xy':cached['xy_um']}
    validate_external_alignment(cached,real)
    with pytest.raises(ValueError,match='float64'):
        validate_external_alignment({**cached,'targets':y.astype(np.float32)},real)


def test_recovered_runtime_comes_from_archived_metadata(tmp_path):
    from graphpathology.reproduction_artifacts import reference_threads
    reference = tmp_path / 'PANEL.csv.gz'
    (tmp_path / 'COMPLETE.json').write_text(json.dumps({'threads':1}))
    assert reference_threads(reference) == 1
    (tmp_path / 'COMPLETE.json').write_text('{}')
    assert reference_threads(reference) == 4
    assert reference_threads(None) == 4


@pytest.mark.parametrize('features',[np.zeros((0,512)),np.zeros((3,511)),np.full((3,512),np.nan)])
def test_reproduction_inference_rejects_invalid_feature_inputs(tmp_path,features):
    with pytest.raises(ValueError,match='feature matrix'):
        image_predictions(tmp_path/'missing-model',features)
