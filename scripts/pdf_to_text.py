#!/usr/bin/env python3
"""pdf_to_text.py

Robust PDF -> text extractor. Tries the following methods in order:
 1. PyPDF2
 2. pdfplumber
 3. PyMuPDF (fitz)
 4. pdftotext command-line utility

Usage:
  python scripts/pdf_to_text.py /path/to/input.pdf [/path/to/output.txt]

If output path is omitted, the script writes <input>.txt next to the PDF.
"""
from __future__ import annotations
import sys
from pathlib import Path
import shutil


def extract_with_pypdf2(path: Path) -> str:
    try:
        import PyPDF2
    except Exception:
        raise
    text_parts = []
    with path.open('rb') as f:
        reader = PyPDF2.PdfReader(f)
        for p in reader.pages:
            t = p.extract_text()
            if t:
                text_parts.append(t)
    return "\n\n".join(text_parts)


def extract_with_pdfplumber(path: Path) -> str:
    try:
        import pdfplumber
    except Exception:
        raise
    parts = []
    with pdfplumber.open(path) as pdf:
        for p in pdf.pages:
            t = p.extract_text()
            if t:
                parts.append(t)
    return "\n\n".join(parts)


def extract_with_fitz(path: Path) -> str:
    try:
        import fitz
    except Exception:
        raise
    doc = fitz.open(path)
    parts = []
    for p in doc:
        parts.append(p.get_text())
    doc.close()
    return "\n\n".join(parts)


def extract_with_pdftotext_cli(path: Path) -> str:
    bin_path = shutil.which('pdftotext')
    if not bin_path:
        raise FileNotFoundError('pdftotext not found')
    import subprocess
    out = subprocess.check_output([bin_path, '-layout', str(path), '-']).decode('utf-8', errors='ignore')
    return out


def normalize_text(s: str) -> str:
    # Simple whitespace normalization while keeping paragraphs
    lines = [ln.rstrip() for ln in s.splitlines()]
    # Collapse sequences of blank lines to a single blank line
    out_lines = []
    blank = False
    for ln in lines:
        if not ln.strip():
            if not blank:
                out_lines.append('')
            blank = True
        else:
            out_lines.append(ln)
            blank = False
    return '\n'.join(out_lines).strip() + '\n'


def main(argv=None):
    import argparse
    parser = argparse.ArgumentParser(description='Extract text from a PDF to a .txt file')
    parser.add_argument('input', type=Path, help='Input PDF path')
    parser.add_argument('output', type=Path, nargs='?', help='Output text path (optional)')
    args = parser.parse_args(argv)

    inp = Path(args.input)
    if not inp.exists():
        print('Input file not found:', inp, file=sys.stderr)
        return 2
    out = Path(args.output) if args.output else inp.with_suffix('.txt')

    extractors = [
        ('PyPDF2', extract_with_pypdf2),
        ('pdfplumber', extract_with_pdfplumber),
        ('PyMuPDF (fitz)', extract_with_fitz),
        ('pdftotext CLI', extract_with_pdftotext_cli),
    ]

    last_exc = None
    text = None
    for name, fn in extractors:
        try:
            print(f'Trying extractor: {name}...', file=sys.stderr)
            text = fn(inp)
            if text and text.strip():
                print(f'Extractor succeeded: {name}', file=sys.stderr)
                break
        except Exception as e:
            last_exc = e
            print(f'Extractor {name} failed: {e}', file=sys.stderr)
            continue

    if not text:
        print('All extractors failed. Last error:', last_exc, file=sys.stderr)
        return 3

    text = normalize_text(text)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(text, encoding='utf-8')
    print('Wrote output to', out)
    return 0

if __name__ == '__main__':
    raise SystemExit(main())
