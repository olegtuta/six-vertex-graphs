"""Verify the deposited certificates, including the two unrestricted bounds."""
import argparse
from fractions import Fraction
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parent


def verify_files():
    for line in (ROOT/'SHA256SUMS.txt').read_text(encoding='ascii').splitlines():
        digest, name = line.split('  ', 1)
        if sha256((ROOT/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Checksum mismatch: '+name)


def verify_rook_upper(work, compiler):
    spec = spec_from_file_location('flag_grammar',ROOT/'f26/unrestricted/grammar.py')
    grammar = module_from_spec(spec)
    spec.loader.exec_module(grammar)
    factors = ROOT/'f74/factors.txt'
    checked = grammar.check(factors)
    assert checked['integer_gram_denominator'] == 2**40
    cxx = compiler or os.environ.get('CXX') or next(
        (c for c in ('c++','g++','cl','clang++') if shutil.which(c)), None)
    if not cxx:
        raise SystemExit('Set --compiler to a C++20 compiler executable (see README.txt).')
    binary = work/('rook.exe' if os.name == 'nt' else 'rook')
    source = ROOT/'f74/eight.cpp'
    if Path(cxx).stem.lower() == 'cl':
        args = ['/nologo','/O2','/EHsc','/std:c++20',source,'/Fe:'+str(binary)]
    else:
        args = ['-O3','-std=c++20',source,'-o',binary]
    subprocess.run([cxx,*args],cwd=work,check=True)
    upper = Fraction(11321089900811290113624899,232113757366008801543585792)
    denominator = 40320*2**80
    assert denominator % upper.denominator == 0
    result = work/'rook.json'
    subprocess.run([binary,factors,str(int(upper*denominator)),result,work/'progress'],cwd=work,check=True)
    report = json.loads(result.read_text(encoding='utf-8'))
    assert report['host_classes'] == 12346 and report['labelled_graphs_covered'] == 2**28
    assert Fraction(int(report['numerator']),int(report['denominator'])) == upper
    assert 0 < upper-Fraction(360,7381) < Fraction(695,10**11)
    return report


def run():
    if not __debug__:
        raise SystemExit('Use ordinary Python with assertions enabled.')
    sys.dont_write_bytecode = True
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--part', choices=['all','cases','f26','f74-upper'],default='all')
    p.add_argument('--compiler',help='compiler executable, e.g. cl or g++')
    args = p.parse_args()
    verify_files()
    with tempfile.TemporaryDirectory(prefix='six-vertex-') as temp:
        work = Path(temp)
        if args.part in ('all','cases'):
            output = subprocess.check_output([sys.executable,'-I','-B',ROOT/'verify_cases.py'],cwd=work)
            result = json.loads(output)
            assert result == json.loads((ROOT/'expected-cases.json').read_text(encoding='utf-8'))
            print(json.dumps(result,indent=2),flush=True)
        if args.part in ('all','f74-upper'):
            print(json.dumps(verify_rook_upper(work,args.compiler),indent=2),flush=True)
        if args.part in ('all','f26'):
            receipt = work/'f26-results.json'
            cmd = [sys.executable,'-I','-B',ROOT/'f26/verify.py','--output',receipt]
            if args.compiler:
                cmd.extend(['--compiler',args.compiler])
            subprocess.run(cmd,cwd=work,check=True)
            assert json.loads(receipt.read_text(encoding='utf-8')) == json.loads(
                (ROOT/'f26/reference-results.json').read_text(encoding='utf-8'))
    print('Verified: '+args.part)


if __name__ == '__main__':
    run()
