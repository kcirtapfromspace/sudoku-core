#!/usr/bin/env python3
"""Check authored Markdown structure, local links, and STE sentence limits.

This conservative prose check does not replace the official word dictionary or
a technical review. Code, equations in code fences, tables, and source titles
are outside its sentence count.
"""
import argparse
from collections import Counter
from pathlib import Path
import re
import subprocess
from urllib.parse import unquote, urlsplit

ROOT = Path(__file__).resolve().parent.parent
FENCE = re.compile(r"^ {0,3}(`{3,}|~{3,})(.*)$")
LINK = re.compile(r"!?\[([^\]]*)\]\(([^\n]*?)\)")
INSTRUCTION = re.compile(
    r"^(?:Run|Use|Do|Set|Keep|Start|Write|Read|Check|Install|Add|Select|Compare|Record|"
    r"Preserve|Give|Assume|Let|Take|Replace|Multiply|Restrict|Require|Reject|Bind|"
    r"Reconstruct|Enumerate|Choose|Remove|Restore|Apply|Verify|Test|Limit|Define|"
    r"Include|Put|Make|Supply|Export|Repeat|Stop|Count|Solve)\b"
)
CONTRACTION = re.compile(r"\b(?:\w+n['’]t|it's|that's|we're|they're|you'll|we'll|there's)\b", re.I)


def documents():
    result = subprocess.check_output(
        ["git", "ls-files", "--cached", "--others", "--exclude-standard", "--", "*.md", "*.mdx"],
        cwd=ROOT, text=True)
    return sorted(set(ROOT / line for line in result.splitlines()))


def plain(text):
    text = re.sub(r"`[^`]+`", "IDENTIFIER", text)
    text = re.sub(r"\$[^$\n]+\$", "EQUATION", text)
    text = LINK.sub(lambda m: m[1], text)
    text = re.sub(r"https?://\S+", "URL", text)
    text = text.replace("e.g.", "for example").replace("i.e.", "that is")
    return re.sub(r"[*_#]", "", text)


def sentences(text):
    return [s for s in re.split(r"(?<=[.!?])\s+", plain(text)) if s.strip()]


def units(lines):
    opened, block, start, display_math = None, [], 0, False
    for number, line in enumerate(lines, 1):
        marker = FENCE.match(line)
        if marker:
            if block:
                yield start, " ".join(block)
                block = []
            token, info = marker.groups()
            if opened is None:
                opened = (token[0], len(token))
            elif token[0] == opened[0] and len(token) >= opened[1] and not info.strip():
                opened = None
            continue
        if opened:
            continue
        if line.strip() == "$$":
            if block:
                yield start, " ".join(block)
                block = []
            display_math = not display_math
            continue
        if display_math:
            continue
        boundary = not line.strip() or line.startswith(("#", "|", ">", "![", "---"))
        boundary |= bool(re.match(r"^\s*\[[^]]+\]:", line))
        item = bool(re.match(r"^\s*(?:[-*+] |\d+[.)] )", line))
        if boundary or item:
            if block:
                yield start, " ".join(block)
                block = []
            if boundary:
                continue
        if not block:
            start = number
        block.append(re.sub(r"^\s*(?:[-*+] |\d+[.)] )", "", line.strip()))
    if block:
        yield start, " ".join(block)
    if opened:
        yield len(lines), "UNCLOSED_FENCE"
    if display_math:
        yield len(lines), "UNCLOSED_MATH"


def check(path, descriptive=25, procedural=20):
    errors = []
    lines = path.read_text().splitlines()
    for line, text in units(lines):
        if text in ("UNCLOSED_FENCE", "UNCLOSED_MATH"):
            errors.append((line, "fence", text.replace("_", " ").capitalize()))
            continue
        # A standalone bibliographic title is source metadata.
        if re.fullmatch(r"\[[^]]+\]\([^\n]+\)[.:]?", text):
            continue
        parts = sentences(text)
        if len(parts) > 6:
            errors.append((line, "paragraph", f"{len(parts)} sentences; maximum 6"))
        for sentence in parts:
            words = len(re.findall(r"\S+", sentence))
            limit = procedural if INSTRUCTION.match(sentence) else descriptive
            if words > limit:
                errors.append((line, "sentence", f"{words} words; maximum {limit}: {sentence}"))
            if CONTRACTION.search(sentence):
                errors.append((line, "contraction", sentence))
    opened = None
    for line, text in enumerate(lines, 1):
        marker = FENCE.match(text)
        if marker:
            token, info = marker.groups()
            if opened is None:
                opened = (token[0], len(token))
            elif token[0] == opened[0] and len(token) >= opened[1] and not info.strip():
                opened = None
            continue
        if opened:
            continue
        for match in LINK.finditer(text):
            target = match[2].strip().split(' "', 1)[0].strip("<>")
            url = urlsplit(target)
            if url.scheme or url.netloc or not url.path:
                continue
            destination = ROOT / unquote(url.path[1:]) if url.path.startswith("/") else path.parent / unquote(url.path)
            if not destination.exists():
                errors.append((line, "link", f"Missing local target: {target}"))
    return errors


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sentences", type=int, default=25)
    parser.add_argument("--instructions", type=int, default=20)
    args = parser.parse_args()
    failures = []
    files = documents()
    for path in files:
        for line, kind, message in check(path, args.sentences, args.instructions):
            failures.append((kind, path, line, message))
            print(f"{path.relative_to(ROOT)}:{line}: {kind}: {message}")
    print(f"Checked {len(files)} documents; {len(failures)} findings: {dict(Counter(x[0] for x in failures))}")
    raise SystemExit(bool(failures))


if __name__ == "__main__":
    main()
