# -*- mode: python ; coding: utf-8 -*-

from PyGFETdb import __version__

block_cipher = None

added_files = [
    ('PyGFETdb\\GuiDBView\\GuiDBView_v2.ui', '.'),
    ('PyGFETdb\\GuiDBView\\GuiDataExplorer_v2.ui', 'PyGFETdb\\GuiDBView'),
    ('PyGFETdb\\GuiDBView\\GuiNormalization.ui', 'PyGFETdb\\GuiDBView'),
    ('PyGFETdb\\GuiDBView\\key.key', '.'),
    ('PyGFETdb\\Connection.bin', 'PyGFETdb'),
]

a = Analysis(
    ['PyGFETdb\\GuiDBView\\GuiDBView_v2.py'],
    pathex=[],
    binaries=[],
    datas=added_files,
    hiddenimports=[],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)
pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.zipfiles,
    a.datas,
    [],
    name='GuiDBView_v' + __version__,
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=True,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)
