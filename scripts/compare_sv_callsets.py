#!/usr/bin/env python3

import argparse
import gzip
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Event:
    chrom: str
    pos: int
    event_type: str
    size: int
    event_end: int
    record_id: str


def read_events(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    events = []
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            delta = len(fields[4]) - len(fields[3])
            event_type = "INS" if delta > 0 else "DEL"
            size = abs(delta)
            pos = int(fields[1])
            event_end = pos if event_type == "INS" else pos + size
            events.append(Event(fields[0], pos, event_type, size, event_end, fields[2]))
    return events


def candidate_edges(left, right, distance, minimum_size_ratio):
    edges = [[] for _ in left]
    for left_index, first in enumerate(left):
        for right_index, second in enumerate(right):
            if first.chrom != second.chrom or first.event_type != second.event_type:
                continue
            if abs(first.pos - second.pos) > distance:
                continue
            if abs(first.event_end - second.event_end) > distance:
                continue
            if min(first.size, second.size) / max(first.size, second.size) < minimum_size_ratio:
                continue
            edges[left_index].append(right_index)
    return edges


def maximum_matching(edges, right_count):
    left_match = [-1] * len(edges)
    right_match = [-1] * right_count

    def augment(left_index, visited):
        for right_index in edges[left_index]:
            if visited[right_index]:
                continue
            visited[right_index] = True
            if right_match[right_index] == -1 or augment(right_match[right_index], visited):
                left_match[left_index] = right_index
                right_match[right_index] = left_index
                return True
        return False

    order = sorted(range(len(edges)), key=lambda index: len(edges[index]))
    for left_index in order:
        augment(left_index, [False] * right_count)
    return left_match


def main():
    parser = argparse.ArgumentParser(description="Compare sequence-resolved SV callsets.")
    parser.add_argument("--left", required=True, type=Path)
    parser.add_argument("--right", required=True, type=Path)
    parser.add_argument("--distance", required=True, type=int)
    parser.add_argument("--minimum-size-ratio", type=float, default=0.7)
    parser.add_argument("--output-matches", required=True, type=Path)
    args = parser.parse_args()

    left = read_events(args.left)
    right = read_events(args.right)
    edges = candidate_edges(left, right, args.distance, args.minimum_size_ratio)
    matching = maximum_matching(edges, len(right))
    matched = [(index, right_index) for index, right_index in enumerate(matching) if right_index >= 0]

    args.output_matches.parent.mkdir(parents=True, exist_ok=True)
    with args.output_matches.open("w") as handle:
        handle.write(
            "left_pos\tright_pos\ttype\tleft_size\tright_size\t"
            "start_distance\tend_distance\tsize_ratio\tleft_id\tright_id\n"
        )
        for left_index, right_index in matched:
            first = left[left_index]
            second = right[right_index]
            handle.write(
                f"{first.pos}\t{second.pos}\t{first.event_type}\t{first.size}\t"
                f"{second.size}\t{abs(first.pos - second.pos)}\t"
                f"{abs(first.event_end - second.event_end)}\t"
                f"{min(first.size, second.size) / max(first.size, second.size):.6f}\t"
                f"{first.record_id}\t{second.record_id}\n"
            )

    print(f"left_total\t{len(left)}")
    print(f"right_total\t{len(right)}")
    print(f"matched\t{len(matched)}")
    print(f"left_match_rate\t{len(matched) / len(left):.6f}")
    print(f"right_match_rate\t{len(matched) / len(right):.6f}")


if __name__ == "__main__":
    main()
