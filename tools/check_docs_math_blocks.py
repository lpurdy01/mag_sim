#!/usr/bin/env python3
"""Validate block math formatting in Markdown files."""

from __future__ import annotations

import pathlib
import sys
from typing import List


def check_file(path: pathlib.Path) -> List[str]:
    lines = path.read_text().splitlines()
    errors: List[str] = []

    in_dollar_block = False
    square_depth = 0

    for idx, line in enumerate(lines):
        stripped = line.strip()
        prev_blank = idx == 0 or lines[idx - 1].strip() == ''
        next_blank = idx == len(lines) - 1 or lines[idx + 1].strip() == ''

        if stripped == '$$':
            if not in_dollar_block:
                if not prev_blank:
                    errors.append(
                        f"{path}:{idx + 1}: block math opener '$$' must be preceded by a blank line"
                    )
                in_dollar_block = True
            else:
                if not next_blank:
                    errors.append(
                        f"{path}:{idx + 1}: block math closer '$$' must be followed by a blank line"
                    )
                in_dollar_block = False
        elif stripped == '\\[':
            if not prev_blank:
                errors.append(
                    f"{path}:{idx + 1}: block math opener '\\[' must be preceded by a blank line"
                )
            square_depth += 1
        elif stripped == '\\]':
            if square_depth == 0:
                errors.append(f"{path}:{idx + 1}: unexpected block math closer '\\]'")
            else:
                square_depth -= 1
            if not next_blank:
                errors.append(
                    f"{path}:{idx + 1}: block math closer '\\]' must be followed by a blank line"
                )

    if in_dollar_block:
        errors.append(f"{path}: unmatched '$$' block delimiter")
    if square_depth != 0:
        errors.append(f"{path}: unmatched '\\['/\\]' delimiters")

    return errors


def main() -> int:
    doc_root = pathlib.Path('docs')
    problems: List[str] = []
    for path in sorted(doc_root.rglob('*.md')):
        problems.extend(check_file(path))

    if problems:
        print('\n'.join(problems), file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
