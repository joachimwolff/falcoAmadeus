#!/usr/bin/env python3
"""Write the Read the Docs page (documentation/docs/index.md) from README.md.

The README is the single source of the documentation. Links to files of the
repository are rewritten to GitHub URLs, and the images the README shows are
copied next to the page, so that `mkdocs build --strict` passes.

    python3 documentation/sync_from_readme.py          # write
    python3 documentation/sync_from_readme.py --check  # exit 1 if out of date
"""
import pathlib
import re
import shutil
import sys

REPO = "https://github.com/joachimwolff/falcoAmadeus"
ROOT = pathlib.Path(__file__).resolve().parent.parent
DOCS = ROOT / "documentation" / "docs"


def render(readme):
    """Return the page text and the image paths it needs."""
    images = sorted(set(re.findall(r"\]\((img/[^)]+)\)", readme)))

    def link(match):
        target = match.group(1)
        if target.startswith(("http://", "https://", "#", "img/")):
            return match.group(0)
        kind = "tree" if target.endswith("/") else "blob"
        return f"]({REPO}/{kind}/master/{target})"

    page = re.sub(r"\]\(([^)\s]+)\)", link, readme)
    header = ("<!-- Generated from README.md by "
              "documentation/sync_from_readme.py; do not edit. -->\n\n")
    return header + page, images


def main():
    check = "--check" in sys.argv[1:]
    page, images = render((ROOT / "README.md").read_text())
    target = DOCS / "index.md"
    stale = [] if target.exists() and target.read_text() == page \
        else [target]
    for image in images:
        src, dst = ROOT / image, DOCS / image
        if not dst.exists() or dst.read_bytes() != src.read_bytes():
            stale.append(dst)
            if not check:
                dst.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(src, dst)
    if check:
        for path in stale:
            print(f"out of date: {path.relative_to(ROOT)}", file=sys.stderr)
        if stale:
            print("run: python3 documentation/sync_from_readme.py",
                  file=sys.stderr)
        return 1 if stale else 0
    target.write_text(page)
    return 0


if __name__ == "__main__":
    sys.exit(main())
