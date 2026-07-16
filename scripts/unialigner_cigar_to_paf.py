#!/usr/bin/env python3
"""Convert a UniAligner global CIGAR to an exact-match PAF record."""

from __future__ import annotations

import argparse
import re
from pathlib import Path


CIGAR_RE = re.compile(r"(\d+)([MID=X])")


def read_single_fasta(path: Path) -> tuple[str, str]:
    name = None
    sequence = []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    raise ValueError(f"Expected one FASTA record in {path}")
                name = line[1:].split()[0]
            else:
                if name is None:
                    raise ValueError(f"Sequence precedes header in {path}")
                sequence.append(line.upper())
    if name is None:
        raise ValueError(f"No FASTA record in {path}")
    return name, "".join(sequence)


def append_operation(operations: list[tuple[int, str]], length: int, operation: str) -> None:
    if length == 0:
        return
    if operations and operations[-1][1] == operation:
        operations[-1] = (operations[-1][0] + length, operation)
    else:
        operations.append((length, operation))


def convert_cigar(cigar: str, target: str, query: str) -> tuple[str, int, int]:
    parsed = CIGAR_RE.findall(cigar)
    if not parsed or "".join(f"{length}{op}" for length, op in parsed) != cigar:
        raise ValueError(f"Unsupported CIGAR: {cigar}")

    target_pos = 0
    query_pos = 0
    matches = 0
    block_length = 0
    exact_operations: list[tuple[int, str]] = []
    for length_text, operation in parsed:
        length = int(length_text)
        block_length += length
        if operation in "=X":
            append_operation(exact_operations, length, operation)
            if operation == "=":
                matches += length
            target_pos += length
            query_pos += length
        elif operation == "M":
            run_operation = None
            run_length = 0
            for offset in range(length):
                exact_operation = "=" if target[target_pos + offset] == query[query_pos + offset] else "X"
                if exact_operation != run_operation:
                    append_operation(exact_operations, run_length, run_operation or "=")
                    run_operation = exact_operation
                    run_length = 0
                run_length += 1
                if exact_operation == "=":
                    matches += 1
            append_operation(exact_operations, run_length, run_operation or "=")
            target_pos += length
            query_pos += length
        elif operation == "I":
            append_operation(exact_operations, length, operation)
            query_pos += length
        elif operation == "D":
            append_operation(exact_operations, length, operation)
            target_pos += length

    if target_pos != len(target) or query_pos != len(query):
        raise ValueError(
            "CIGAR consumption mismatch: "
            f"target {target_pos}/{len(target)}, query {query_pos}/{len(query)}"
        )
    exact_cigar = "".join(f"{length}{operation}" for length, operation in exact_operations)
    return exact_cigar, matches, block_length


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", type=Path, required=True)
    parser.add_argument("--query", type=Path, required=True)
    parser.add_argument("--cigar", type=Path, required=True)
    args = parser.parse_args()

    target_name, target = read_single_fasta(args.target)
    query_name, query = read_single_fasta(args.query)
    cigar = args.cigar.read_text().strip()
    exact_cigar, matches, block_length = convert_cigar(cigar, target, query)
    fields = [
        query_name,
        len(query),
        0,
        len(query),
        "+",
        target_name,
        len(target),
        0,
        len(target),
        matches,
        block_length,
        255,
        f"cg:Z:{exact_cigar}",
    ]
    print("\t".join(map(str, fields)))


if __name__ == "__main__":
    main()
