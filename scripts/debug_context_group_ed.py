#!/usr/bin/env python3

import argparse
import csv
import math
from collections import Counter, defaultdict
from pathlib import Path


def read_grouped_anchors(path):
    anchors = {}
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            anchors[row["anchor_id"]] = row
    return anchors


def read_anchor_contexts(path):
    contexts = {}
    alpha_columns = []
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        alpha_columns = [
            name for name in reader.fieldnames if name.endswith("_alpha_prime")
        ]
        for row in reader:
            row["tokens"] = (
                row["forward_context_tokens"].split("|")
                if row["forward_context_tokens"]
                else []
            )
            contexts[row["anchor_id"]] = row
    return contexts, alpha_columns


def read_unmatched_assignments(path):
    if path is None:
        return []
    rows = []
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["status"] == "unmatched":
                rows.append(row)
    return rows


def edit_distance(left, right, limit=None):
    if limit is not None and abs(len(left) - len(right)) > limit:
        return limit + 1
    previous = list(range(len(right) + 1))
    current = [0] * (len(right) + 1)
    for i, left_token in enumerate(left, 1):
        current[0] = i
        row_min = current[0]
        for j, right_token in enumerate(right, 1):
            current[j] = min(
                previous[j - 1] + (0 if left_token == right_token else 1),
                previous[j] + 1,
                current[j - 1] + 1,
            )
            row_min = min(row_min, current[j])
        if limit is not None and row_min > limit:
            return limit + 1
        previous, current = current, previous
    result = previous[len(right)]
    if limit is not None and result > limit:
        return limit + 1
    return result


def traceback_alignment(left, right):
    n = len(left)
    m = len(right)
    dp = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        dp[i][0] = i
    for j in range(1, m + 1):
        dp[0][j] = j
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            dp[i][j] = min(
                dp[i - 1][j - 1] + (0 if left[i - 1] == right[j - 1] else 1),
                dp[i - 1][j] + 1,
                dp[i][j - 1] + 1,
            )

    rows = []
    i = n
    j = m
    while i or j:
        if i and j and dp[i][j] == dp[i - 1][j - 1] + (
            0 if left[i - 1] == right[j - 1] else 1
        ):
            op = "match" if left[i - 1] == right[j - 1] else "substitution"
            rows.append((op, left[i - 1], right[j - 1]))
            i -= 1
            j -= 1
        elif i and dp[i][j] == dp[i - 1][j] + 1:
            rows.append(("deletion", left[i - 1], "."))
            i -= 1
        else:
            rows.append(("insertion", ".", right[j - 1]))
            j -= 1
    rows.reverse()
    return dp[n][m], rows


def token_type(token):
    if token == ".":
        return "."
    if token.startswith("A:"):
        return "A"
    if token.startswith("D:"):
        return "D"
    return "?"


def alpha_value(context, sample_name):
    return int(context.get(f"{sample_name}_alpha_prime", "0"))


def build_indexes(grouped, contexts):
    by_sample_group = defaultdict(list)
    group_counts = Counter()
    sample_group_counts = Counter()
    for anchor_id, context in contexts.items():
        group_id = grouped[anchor_id]["anchor_group_id"]
        if group_id == ".":
            continue
        sample = context["sample"]
        by_sample_group[(sample, group_id)].append(anchor_id)
        group_counts[group_id] += 1
        sample_group_counts[(sample, group_id)] += 1
    return by_sample_group, group_counts, sample_group_counts


def score_candidate_group(group_key, unmatched_rows, by_sample_group, contexts,
                          grouped, auto_ed_limit):
    target_sample, query_sample, group_id = group_key
    target_anchors = by_sample_group.get((target_sample, group_id), [])
    if not target_anchors:
        return None

    sampled = unmatched_rows[: min(2, len(unmatched_rows))]
    best_eds = []
    for row in sampled:
        query = contexts[row["query_anchor_id"]]
        query_tokens = query["tokens"]
        reverse_tokens = list(reversed(query_tokens))
        best_ed = None
        for target_id in target_anchors:
            target = contexts[target_id]
            forward_ed = edit_distance(query_tokens, target["tokens"], auto_ed_limit)
            reverse_ed = edit_distance(reverse_tokens, target["tokens"], auto_ed_limit)
            ed = min(forward_ed, reverse_ed)
            if best_ed is None or ed < best_ed:
                best_ed = ed
        if best_ed is not None:
            best_eds.append(best_ed)
    if not best_eds:
        return None
    mean_ed = sum(best_eds) / len(best_eds)
    return {
        "target_sample": target_sample,
        "query_sample": query_sample,
        "group_id": group_id,
        "unmatched_count": len(unmatched_rows),
        "target_count": len(target_anchors),
        "mean_best_ed": mean_ed,
        "sampled_best_eds": ",".join(str(value) for value in best_eds),
    }


def choose_group(args, assignments, by_sample_group, contexts, grouped):
    grouped_unmatched = defaultdict(list)
    for row in assignments:
        if args.query_sample and row["query_sample"] != args.query_sample:
            continue
        if args.target_sample and row["target_sample"] != args.target_sample:
            continue
        group_id = row["query_anchor_group_id"]
        if group_id == ".":
            continue
        target_count = len(by_sample_group.get((row["target_sample"], group_id), []))
        if target_count < args.min_target_group_count:
            continue
        key = (row["target_sample"], row["query_sample"], group_id)
        grouped_unmatched[key].append(row)

    prescreened = []
    for key, rows in grouped_unmatched.items():
        if len(rows) < args.min_unmatched_count:
            continue
        target_count = len(by_sample_group.get((key[0], key[2]), []))
        prescreened.append((len(rows), target_count, key, rows))

    prescreened.sort(reverse=True)
    candidates = []
    for _, _, key, rows in prescreened[: args.max_auto_groups]:
        score = score_candidate_group(
            key, rows, by_sample_group, contexts, grouped, args.auto_ed_limit)
        if score is not None:
            candidates.append(score)

    if not candidates:
        raise SystemExit("no suitable unmatched high-copy group found")

    candidates.sort(
        key=lambda row: (
            row["mean_best_ed"],
            -row["unmatched_count"],
            -row["target_count"],
        )
    )
    selected = candidates[0]
    query_anchor = grouped_unmatched[
        (selected["target_sample"], selected["query_sample"], selected["group_id"])
    ][0]["query_anchor_id"]
    return selected, query_anchor, candidates[: args.report_group_candidates]


def compute_summary(query_anchor, target_sample, grouped, contexts, by_sample_group):
    query = contexts[query_anchor]
    query_sample = query["sample"]
    group_id = grouped[query_anchor]["anchor_group_id"]
    query_tokens = query["tokens"]
    reverse_tokens = list(reversed(query_tokens))
    q_alpha = alpha_value(query, target_sample)

    rows = []
    for target_anchor in by_sample_group.get((target_sample, group_id), []):
        target = contexts[target_anchor]
        forward_ed = edit_distance(query_tokens, target["tokens"])
        reverse_ed = edit_distance(reverse_tokens, target["tokens"])
        best_ed = min(forward_ed, reverse_ed)
        strand = "-" if reverse_ed < forward_ed else "+"
        r_alpha = alpha_value(target, query_sample)
        candidate_score = q_alpha + r_alpha - best_ed
        grouped_target = grouped[target_anchor]
        rows.append(
            {
                "query_anchor_id": query_anchor,
                "query_sample": query_sample,
                "target_anchor_id": target_anchor,
                "target_sample": target_sample,
                "anchor_group_id": group_id,
                "target_sequence_id": target["sequence_id"],
                "target_start": grouped_target["start"],
                "target_end": grouped_target["end"],
                "target_center_token_index": target["center_token_index"],
                "forward_ed": forward_ed,
                "reverse_ed": reverse_ed,
                "best_ed": best_ed,
                "assign_strand": strand,
                "q_alpha_prime": q_alpha,
                "r_alpha_prime": r_alpha,
                "candidate_score": candidate_score,
                "query_context_key": query["canonical_context_key"],
                "target_context_key": target["canonical_context_key"],
                "query_context_tokens": query["forward_context_tokens"],
                "target_context_tokens": target["forward_context_tokens"],
            }
        )
    rows.sort(
        key=lambda row: (
            -int(row["candidate_score"]),
            int(row["best_ed"]),
            row["target_anchor_id"],
        )
    )
    return rows


def write_summary(path, rows):
    fieldnames = [
        "rank",
        "query_anchor_id",
        "query_sample",
        "target_anchor_id",
        "target_sample",
        "anchor_group_id",
        "target_sequence_id",
        "target_start",
        "target_end",
        "target_center_token_index",
        "forward_ed",
        "reverse_ed",
        "best_ed",
        "assign_strand",
        "q_alpha_prime",
        "r_alpha_prime",
        "candidate_score",
        "query_context_key",
        "target_context_key",
        "query_context_tokens",
        "target_context_tokens",
    ]
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for rank, row in enumerate(rows):
            output = dict(row)
            output["rank"] = rank
            writer.writerow(output)


def write_alignment(path, rows, contexts, top_alignments):
    fieldnames = [
        "rank",
        "alignment_index",
        "operation",
        "query_anchor_id",
        "target_anchor_id",
        "assign_strand",
        "query_token",
        "target_token",
        "query_token_type",
        "target_token_type",
    ]
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for rank, row in enumerate(rows[:top_alignments]):
            query_tokens = contexts[row["query_anchor_id"]]["tokens"]
            if row["assign_strand"] == "-":
                query_tokens = list(reversed(query_tokens))
            target_tokens = contexts[row["target_anchor_id"]]["tokens"]
            _, alignment = traceback_alignment(query_tokens, target_tokens)
            for index, (operation, query_token, target_token) in enumerate(alignment):
                writer.writerow(
                    {
                        "rank": rank,
                        "alignment_index": index,
                        "operation": operation,
                        "query_anchor_id": row["query_anchor_id"],
                        "target_anchor_id": row["target_anchor_id"],
                        "assign_strand": row["assign_strand"],
                        "query_token": query_token,
                        "target_token": target_token,
                        "query_token_type": token_type(query_token),
                        "target_token_type": token_type(target_token),
                    }
                )


def write_group_candidates(path, rows):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "rank",
        "target_sample",
        "query_sample",
        "group_id",
        "unmatched_count",
        "target_count",
        "mean_best_ed",
        "sampled_best_eds",
    ]
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for rank, row in enumerate(rows):
            output = dict(row)
            output["rank"] = rank
            output["mean_best_ed"] = f"{row['mean_best_ed']:.3f}"
            writer.writerow(output)


def main():
    parser = argparse.ArgumentParser(
        description="Dump same-anchor-group structural contexts and token EDs."
    )
    parser.add_argument("--grouped-anchors", default="grouped_anchors.tsv")
    parser.add_argument("--anchor-contexts", default="anchor_contexts.tsv")
    parser.add_argument("--assignments", default="assignments.tsv")
    parser.add_argument("--query-anchor")
    parser.add_argument("--target-sample")
    parser.add_argument("--query-sample")
    parser.add_argument("--summary-output", required=True)
    parser.add_argument("--alignment-output", required=True)
    parser.add_argument("--group-candidates-output")
    parser.add_argument("--top-alignments", type=int, default=20)
    parser.add_argument("--min-target-group-count", type=int, default=1000)
    parser.add_argument("--min-unmatched-count", type=int, default=5)
    parser.add_argument("--max-auto-groups", type=int, default=10)
    parser.add_argument("--auto-ed-limit", type=int, default=40)
    parser.add_argument("--report-group-candidates", type=int, default=20)
    args = parser.parse_args()

    grouped = read_grouped_anchors(args.grouped_anchors)
    contexts, _ = read_anchor_contexts(args.anchor_contexts)
    by_sample_group, _, _ = build_indexes(grouped, contexts)

    selected = None
    candidates = []
    query_anchor = args.query_anchor
    target_sample = args.target_sample
    if query_anchor is None:
        assignments = read_unmatched_assignments(args.assignments)
        selected, query_anchor, candidates = choose_group(
            args, assignments, by_sample_group, contexts, grouped
        )
        target_sample = selected["target_sample"]
        if args.group_candidates_output:
            write_group_candidates(args.group_candidates_output, candidates)
    elif target_sample is None:
        raise SystemExit("--target-sample is required when --query-anchor is used")

    if query_anchor not in contexts:
        raise SystemExit(f"query anchor not found in anchor_contexts: {query_anchor}")

    rows = compute_summary(query_anchor, target_sample, grouped, contexts, by_sample_group)
    if not rows:
        raise SystemExit("no target anchors found for selected group")

    Path(args.summary_output).parent.mkdir(parents=True, exist_ok=True)
    Path(args.alignment_output).parent.mkdir(parents=True, exist_ok=True)
    write_summary(args.summary_output, rows)
    write_alignment(args.alignment_output, rows, contexts, args.top_alignments)

    print("query_anchor", query_anchor)
    print("query_sample", contexts[query_anchor]["sample"])
    print("target_sample", target_sample)
    print("anchor_group_id", grouped[query_anchor]["anchor_group_id"])
    print("target_candidate_count", len(rows))
    print("best_candidate_score", rows[0]["candidate_score"])
    print("best_ed", rows[0]["best_ed"])
    print("summary_output", args.summary_output)
    print("alignment_output", args.alignment_output)
    if selected:
        print("selected_unmatched_count", selected["unmatched_count"])
        print("selected_mean_best_ed", f"{selected['mean_best_ed']:.3f}")


if __name__ == "__main__":
    main()
