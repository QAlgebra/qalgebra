"""Script for building the documentation in PDF format.

Usage:

    python docs/build_pdf.py TEXFILE

where TEXFILE is produced by Sphinx:

    tox -e docs -- -b latex _build/latex
"""

import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).parent


def get_version(filename):
    """Extract the package version."""
    with open(filename, encoding='utf8') as in_fh:
        for line in in_fh:
            if line.startswith('__version__'):
                return line.split('=')[1].strip()[1:-1]
    raise ValueError("Cannot extract version from %s" % filename)


def _patch_line(line):
    if line.startswith(r'\sphinxhref{') and line.endswith(".svg}}\n"):
        return None
    if line.startswith(r'\sphinxhref{') and line.endswith(
        ".svg?logo=github}}\n"
    ):
        return None
    if line == r'\chapter{References}' + "\n":
        return None
    line = line.replace(r'\sphinxhyphen{}', '-')
    if line.startswith(r'\(\newcommand'):
        return None
    if line.startswith(r'\section{'):
        # don't put section numbers in the HISTORY, in front of version numbers
        line = line.replace(r'\sphinxhyphen{}', '-')
        match = re.match(
            r'\\section\{(\d+\.\d+\.\d+(-rc\d+)?\s+'
            r'\([\d-]+\)|[\s(]*next version[)\s]*)\}',
            line.strip(),
        )
        if match:
            line = r'\section*{' + match.group(1) + "}\n"
    return line


def patch_tex_lines(texfile):
    """Fix errors line-by-line in the given texfile."""
    with tempfile.TemporaryDirectory() as tmpdir:
        orig = Path(tmpdir) / texfile.name
        shutil.copyfile(texfile, orig)
        with orig.open() as in_fh, texfile.open("w") as out_fh:
            for line in in_fh:
                line = _patch_line(line)
                if line is None:
                    continue
                out_fh.write(line)


def _multiline_str(*lines):
    return "\n".join(lines)


def patch_tex(texfile):
    """Fix errors in the given texfile, acting on the whole text."""
    tex = texfile.read_text(encoding='utf8')
    tex = tex.replace(
        r'\begin{equation*}' + "\n" + r'\begin{split}\begin{align}',
        r'\begin{align*}',
    )
    tex = tex.replace(
        r'\end{align}\end{split}' + "\n" + r'\end{equation*}', r'\end{align*}'
    )
    tex = tex.replace(
        r'{\Op{a}_{\rm hs}^\dagger}^2',
        r'\left(\Op{a}_{\rm hs}^\dagger\right)^2',
    )
    tex = tex.replace(
        _multiline_str(
            r'\chapter{Indices and tables}',
            r'\label{\detokenize{index:indices-and-tables}}\begin{itemize}',
            r'\item {} ',
            r'\DUrole{xref,std,std-ref}{genindex}',
            r'',
            r'\item {} ',
            r'\DUrole{xref,std,std-ref}{modindex}',
            r'',
            r'\end{itemize}',
        ),
        '',
    )
    texfile.write_text(tex, encoding='utf8')


def lualatex(texfile):
    """Run lualatex to compile the given texfile."""
    subprocess.run(
        [
            'lualatex',
            '--interaction=nonstopmode',
            '--halt-on-error',
            texfile.name,
        ],
        cwd=texfile.parent,
        check=True,
    )


def main(argv=None):
    """Main function."""
    if argv is None:
        argv = sys.argv
    texfile = argv[-1]
    if not texfile.endswith(".tex"):
        print(__doc__)
        sys.exit(1)
    texfile = Path(texfile)
    if not texfile.is_file():
        print("%s does not exist" % texfile)
        sys.exit(1)
    print("Patching %s..." % texfile)
    patch_tex_lines(texfile)
    patch_tex(texfile)
    if '--patch-only' in argv:
        print("Skipping compilation")
    else:
        print("Compiling %s..." % texfile)
        lualatex(texfile)
        lualatex(texfile)
        print("Done compiling %s" % texfile)
    sys.exit(0)


if __name__ == "__main__":
    main()
