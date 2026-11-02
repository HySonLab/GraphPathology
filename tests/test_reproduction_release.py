import hashlib
from pathlib import Path
import zipfile

import pytest

from graphpathology.assets import sha
from graphpathology.reproduction_release import restore


def fixture_archive(tmp_path, members, expected=None):
    assets=tmp_path/'assets';assets.mkdir()
    path=assets/'weights.zip'
    with zipfile.ZipFile(path,'w') as z:
        for name,data in members.items():z.writestr(name,data)
    manifest={'assets':[{'name':path.name,'bytes':path.stat().st_size,'sha256':sha(path),
                         'files':expected or {n:hashlib.sha256(v).hexdigest() for n,v in members.items()}}]}
    root=tmp_path/'checkout';root.mkdir()
    return root,assets,manifest


def test_restore_is_verified_and_leaves_identical_existing_files(tmp_path):
    name='checkpoints/reproduced/run/seed0.pt'
    root,assets,manifest=fixture_archive(tmp_path,{name:b'head'})
    assert restore(root,assets,manifest)['restored_files']==1
    before=(root/name).stat().st_mtime_ns
    assert restore(root,assets,manifest)['restored_files']==0
    assert (root/name).stat().st_mtime_ns==before


@pytest.mark.parametrize('name',['../escape','/absolute','manuscript/main.tex',
                                 'checkpoints/reproduced/../../manuscript/main.tex'])
def test_restore_rejects_traversal_and_scientific_text_members(tmp_path,name):
    root,assets,manifest=fixture_archive(tmp_path,{name:b'bad'})
    with pytest.raises(ValueError,match='member'):restore(root,assets,manifest)
    assert not list(root.iterdir())


def test_restore_preflights_all_conflicts_before_any_write(tmp_path):
    names={'reproduction/new.txt':b'new','checkpoints/journal_revision/existing.pt':b'expected'}
    root,assets,manifest=fixture_archive(tmp_path,names)
    old=root/'checkpoints/journal_revision/existing.pt';old.parent.mkdir(parents=True);old.write_bytes(b'other')
    with pytest.raises(FileExistsError):restore(root,assets,manifest)
    assert not (root/'reproduction/new.txt').exists()
    assert old.read_bytes()==b'other'


def test_restore_rejects_modified_asset_before_extracting(tmp_path):
    root,assets,manifest=fixture_archive(tmp_path,{'reproduction/new.txt':b'new'})
    with (assets/'weights.zip').open('ab') as f:f.write(b'corrupt')
    with pytest.raises(ValueError,match='checksum'):restore(root,assets,manifest)
    assert not list(root.iterdir())


def test_restore_rejects_inner_member_hash_mismatch(tmp_path):
    root,assets,manifest=fixture_archive(tmp_path,{'reproduction/new.txt':b'new'},
                                        {'reproduction/new.txt':'0'*64})
    with pytest.raises(ValueError,match='Member checksum'):restore(root,assets,manifest)
    assert not list(root.iterdir())


def test_restore_rejects_existing_symlink_escape(tmp_path):
    root,assets,manifest=fixture_archive(tmp_path,{'reproduction/link/new.txt':b'new'})
    outside=tmp_path/'outside';outside.mkdir()
    (root/'reproduction').mkdir();(root/'reproduction/link').symlink_to(outside,target_is_directory=True)
    with pytest.raises(ValueError,match='escapes'):restore(root,assets,manifest)
    assert not list(outside.iterdir())


def test_restore_rejects_unindexed_member(tmp_path):
    root,assets,manifest=fixture_archive(tmp_path,{'reproduction/a':b'a','reproduction/b':b'b'},
                                        {'reproduction/a':hashlib.sha256(b'a').hexdigest()})
    with pytest.raises(ValueError,match='member set'):restore(root,assets,manifest)
