'Relative data paths, output paths and file identifiers.'

from __future__ import annotations

import hashlib

import importlib.util

import json

import os

from pathlib import Path

import sys



PROJECT_ROOT = Path(__file__).resolve().parents[1]

PROJECT = PROJECT_ROOT

os.environ.setdefault('MPLCONFIGDIR', str(Path(os.environ.get('TEMP', str(PROJECT_ROOT/'.cache'))) / 'matplotlib-cache'))

sys.dont_write_bytecode = True

if str(PROJECT) not in sys.path:

    sys.path.insert(0, str(PROJECT))



def output_path(relative: str) -> Path:

    path = (PROJECT_ROOT / relative).resolve()

    if not path.is_relative_to(PROJECT_ROOT.resolve()):

        raise ValueError(f"Output escapes project: {relative}")

    path.parent.mkdir(parents=True, exist_ok=True)

    return path



def sha256(path: Path) -> str:

    h = hashlib.sha256()

    with path.open('rb') as stream:

        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b''):

            h.update(chunk)

    return h.hexdigest()



def write_json(relative: str, payload) -> None:

    output_path(relative).write_text(json.dumps(payload, ensure_ascii=False, indent=2, allow_nan=False) + '\n', encoding='utf-8')



def import_source(name: str, relative: str):

    path = PROJECT / relative

    spec = importlib.util.spec_from_file_location(name, path)

    module = importlib.util.module_from_spec(spec)

    sys.modules[name] = module

    spec.loader.exec_module(module)

    return module
