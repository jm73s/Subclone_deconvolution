#!/usr/bin/env python3

import argparse
import glob
import json
import os
import subprocess
import time
from datetime import datetime

import numpy as np
import pandas as pd


STAGES = [
    "input_generated",
    "fit_done",
    "results_written",
    "cluster_filtered",
    "recomb_filtered",
]


def pyclone_dir(pat):
    return os.path.join(pat, "mut_dyn", "pyclone")


def state_path(pat):
    return os.path.join(pyclone_dir(pat), "iterative_pyclone_state.json")


def iteration_paths(pat, iteration):
    base = pyclone_dir(pat)
    return {
        "input_generated": os.path.join(base, f"pyclone_in_{pat}_iter{iteration}.tsv"),
        "fit_done": os.path.join(base, f"pyclone_out_{pat}_iter{iteration}.h5"),
        "results_written": os.path.join(base, f"pyclone_out_{pat}_iter{iteration}.tsv"),
        "cluster_filtered": os.path.join(base, f"pyclone_out_{pat}_iter{iteration}_cluster.tsv"),
        "recomb_filtered": os.path.join(base, f"pyclone_out_{pat}_iter{iteration}_cluster_recomb.tsv"),
    }


def file_nonempty(path):
    return os.path.exists(path) and os.path.getsize(path) > 0


def validate_tsv_columns(path, required_cols):
    if not file_nonempty(path):
        return False
    try:
        header = pd.read_csv(path, sep="\t", nrows=0)
    except Exception:
        return False
    return all(col in header.columns for col in required_cols)


def validate_stage_output(pat, iteration, stage):
    paths = iteration_paths(pat, iteration)

    if stage == "input_generated":
        return validate_tsv_columns(
            paths[stage],
            [
                "mutation_id",
                "sample_id",
                "ref_counts",
                "alt_counts",
                "major_cn",
                "minor_cn",
                "normal_cn",
                "tumour_content",
                "error_rate",
            ],
        )

    if stage == "fit_done":
        return file_nonempty(paths[stage])

    if stage == "results_written":
        return validate_tsv_columns(
            paths[stage], ["mutation_id", "sample_id", "cluster_id", "cellular_prevalence"]
        )

    if stage in ("cluster_filtered", "recomb_filtered"):
        return validate_tsv_columns(paths[stage], ["mutation_id", "sample_id", "cluster_id"])

    return False


def load_state(pat):
    path = state_path(pat)
    if not os.path.exists(path):
        return None
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def save_state(
    pat,
    params,
    iteration,
    stage,
    status,
    last_error=None,
    removed_clusters_count=None,
    recomb_pos_count=None,
):
    state = {
        "patient": pat,
        "params": params,
        "iteration": iteration,
        "stage": stage,
        "status": status,
        "last_error": last_error,
        "removed_clusters_count": removed_clusters_count,
        "recomb_pos_count": recomb_pos_count,
        "updated_at": datetime.now().isoformat(timespec="seconds"),
    }

    os.makedirs(pyclone_dir(pat), exist_ok=True)
    with open(state_path(pat), "w", encoding="utf-8") as f:
        json.dump(state, f, indent=2)


def params_match(saved_params, current_params):
    return saved_params == current_params


def infer_resume_point(pat):
    iteration = 0
    while True:
        for stage in STAGES:
            if not validate_stage_output(pat, iteration, stage):
                return iteration, stage
        iteration += 1


def cleanup_previous_outputs(pat):
    base = pyclone_dir(pat)
    patterns = [
        f"pyclone_in_{pat}_iter*.tsv",
        f"pyclone_out_{pat}_iter*.h5",
        f"pyclone_out_{pat}_iter*.tsv",
        f"pyclone_out_{pat}_iter*_cluster.tsv",
        f"pyclone_out_{pat}_iter*_cluster_recomb.tsv",
        "iterative_pyclone_state.json",
    ]

    removed = 0
    for pattern in patterns:
        for path in glob.glob(os.path.join(base, pattern)):
            if os.path.isfile(path):
                os.remove(path)
                removed += 1

    print(f"Fresh start: removed {removed} previous output files for {pat}.")


def run_pyclone_fit(pat, max_cluster, restarts, iteration, num_threads=20, seed=15):
    subprocess.run(
        [
            "pyclone-vi",
            "fit",
            "--seed",
            str(seed),
            "-i",
            f"{pat}/mut_dyn/pyclone/pyclone_in_{pat}_iter{iteration}.tsv",
            "-o",
            f"{pat}/mut_dyn/pyclone/pyclone_out_{pat}_iter{iteration}.h5",
            "-c",
            str(max_cluster),
            "-d",
            "beta-binomial",
            "-r",
            str(restarts),
            "-t",
            str(num_threads),
        ],
        check=True,
    )


def write_pyclone_results(pat, iteration):
    subprocess.run(
        [
            "pyclone-vi",
            "write-results-file",
            "-i",
            f"{pat}/mut_dyn/pyclone/pyclone_out_{pat}_iter{iteration}.h5",
            "-o",
            f"{pat}/mut_dyn/pyclone/pyclone_out_{pat}_iter{iteration}.tsv",
        ],
        check=True,
    )


# remove fixed and too low clusters
def cluster_filter(pat, std_thresh, fixed_thresh, neglegible_thresh, iteration, true_mut_df):
    """
    Clean clusters by removing those that are fixed (probably differencence with pao1) and those that are mostly negligible (noise or contamination)
    after ignoring timepoints with high variability.
    pat: patient to consider
    std: threshold for variability to ignore timepoints
    fixed_thresh: threshold above which a cluster is considered fixed
    neglegible_thresh: threshold below which a cluster is considered negligible if only goes above once or less
    iteration: current iteration number
    true_mut_df: dataframe with true mutation dynamics for the patient
    """

    # load pyclone df
    mut_df = pd.read_csv(f"{pat}/mut_dyn/pyclone/pyclone_out_{pat}_iter{iteration}.tsv", sep="\t")
    mut_df_out = mut_df.copy()  # for output later
    # add sample dates
    mut_df["date"] = pd.to_datetime(
        mut_df["sample_id"].str.split("_").str[0], format="%d-%m-%Y", errors="coerce"
    )
    mut_df = mut_df.dropna(subset=["date"])
    mut_df.sort_values(by="date", inplace=True)

    # assign clusters to true mut
    true_mut_df = (
        true_mut_df.merge(
            mut_df[["mutation_id", "cluster_id"]].drop_duplicates(),
            left_on="position",
            right_on="mutation_id",
            how="left",
        )
        .drop("position", axis=1)
        .copy()
    )
    true_mut_df["cluster_id"] = true_mut_df["cluster_id"].astype("Int64")  # allow NaN
    # remove non clustered mutations
    true_mut_df = true_mut_df[true_mut_df["cluster_id"].notnull()].copy()
    # add predicted pyclone cellular prevalence
    true_mut_df = pd.merge(
        true_mut_df,
        mut_df[["mutation_id", "date", "cellular_prevalence"]],
        on=["mutation_id", "date"],
        how="inner",
    )

    # calc std for each cluster at each timepoint from true mut df
    std = true_mut_df.groupby(["cluster_id", "date"])["allele_freq"].std().reset_index()
    std.rename(columns={"allele_freq": "std"}, inplace=True)
    # add to true mut df
    true_mut_df_clean = pd.merge(true_mut_df, std, on=["cluster_id", "date"], how="inner")
    # remove entries with std bigger than threshold
    true_mut_df_clean = true_mut_df_clean[true_mut_df_clean["std"] < std_thresh]
    # find clusters where predicted cellular prevalence never bellow fixed_thresh after removing high std entries
    fixed_clusters = true_mut_df_clean.groupby("cluster_id").apply(
        lambda x: (x["cellular_prevalence"] < fixed_thresh).any() == False
    )
    # find clusters where predicted cellular prevalence only above neglegible_thresh once after removing high std entries
    neglegible_clusters = true_mut_df_clean.groupby("cluster_id").apply(
        lambda x: (x["cellular_prevalence"] > neglegible_thresh).sum() < 2
    )
    remove_clusters = (
        fixed_clusters[fixed_clusters].index.tolist()
        + neglegible_clusters[neglegible_clusters].index.tolist()
    )

    # output cleaned pyclone data
    # remove clusters
    mut_df_out = mut_df_out[~mut_df_out["cluster_id"].isin(remove_clusters)]
    # relable clusters to account for removed ones and align with pairtree numbering
    # get sorted unique cluster labels
    clusters_sorted = sorted(mut_df_out["cluster_id"].unique())
    # build mapping: old_label -> new_label (1..12)
    mapping = {old: new for new, old in enumerate(clusters_sorted, start=1)}
    # relabel
    mut_df_out["cluster_id"] = mut_df_out["cluster_id"].map(mapping)
    # write to file
    with open(
        f"{pat}/mut_dyn/pyclone/pyclone_out_{pat}_iter{iteration}_cluster.tsv", "w", encoding="utf-8"
    ) as f:
        mut_df_out.to_csv(f, sep="\t", index=False)

    return sorted(remove_clusters)


def sliding_recomb_filter(pat, window_size, threshold, iteration, skip=1):
    # load initial pyclone df
    mut_df = pd.read_csv(
        f"{pat}/mut_dyn/pyclone/pyclone_out_{pat}_iter{iteration}_cluster.tsv", sep="\t"
    )
    mut_df["date"] = mut_df["sample_id"].str.split("_").str[0]
    mut_df["date"] = pd.to_datetime(mut_df["date"], format="%d-%m-%Y", errors="coerce")
    mut_df = mut_df.dropna(subset=["date"])

    recomb_pos = set()
    for cluster_id in mut_df["cluster_id"].unique():
        cluster_df = mut_df[mut_df["cluster_id"] == cluster_id].copy()
        cluster_df["position"] = cluster_df["mutation_id"].str.split("_", expand=True)[1].astype(int)
        # get pos for cluster
        c_pos = cluster_df["position"].unique()

        if c_pos.size == 0:
            continue

        # array with length up to biggest pos, with true where mut and false where not
        numline = np.array((c_pos.max() + 1) * [False])
        numline[c_pos] = True

        # sliding window across numline and count number of mutations in each window, if above threshold, add to recomb_pos
        for i in range(1, c_pos.max() - window_size + 1, skip):
            count = numline[i : i + window_size].sum()
            if count > threshold:
                recomb_pos.update(np.where(numline[i : i + window_size])[0] + i)

    # output cleaned pyclone data
    # remove muts in recombination regions
    remove_muts = [f"pos_{pos}" for pos in list(recomb_pos)]
    mut_df = mut_df[~mut_df["mutation_id"].isin(remove_muts)]
    # write to file
    with open(
        f"{pat}/mut_dyn/pyclone/pyclone_out_{pat}_iter{iteration}_cluster_recomb.tsv",
        "w",
        encoding="utf-8",
    ) as f:
        mut_df.to_csv(f, sep="\t", index=False)

    return recomb_pos


def gen_pyclone_df(pat, output_name, iteration):
    # load mutation dynamics df
    pat_dyn = pd.read_csv(f"{pat}/mut_dyn/mut_evolution_dates.txt", sep="\t")
    pat_dyn["date_lane"] = pat_dyn["date"].str.cat(pat_dyn["lane"], sep="_")
    # remove positions where allele_freq is always over 0.9
    always_high = pat_dyn.groupby("position")["allele_freq"].transform(lambda x: (x > 0.9).all())
    # remove positions where allele_freq is above 5% 1 time or less
    only_once_over_5 = pat_dyn.groupby("position")["allele_freq"].transform(
        lambda x: (x > 0.05).sum() < 2
    )
    pat_dyn = pat_dyn.loc[~always_high & ~only_once_over_5].copy()
    df_pyclone = pat_dyn[["position", "date_lane", "DP4"]].copy()
    rename_dict = {"position": "mutation_id", "date_lane": "sample_id"}
    df_pyclone.rename(columns=rename_dict, inplace=True)
    df_pyclone["mutation_id"] = "pos_" + df_pyclone["mutation_id"].astype(str)
    df_pyclone["ref_counts"] = (
        df_pyclone["DP4"].str.split(",", expand=True)[[0, 1]].astype(int).sum(axis=1)
    )
    df_pyclone["alt_counts"] = (
        df_pyclone["DP4"].str.split(",", expand=True)[[2, 3]].astype(int).sum(axis=1)
    )
    df_pyclone["major_cn"] = 1
    df_pyclone["minor_cn"] = 0
    df_pyclone["normal_cn"] = 1
    df_pyclone["tumour_content"] = 1
    df_pyclone["error_rate"] = 0.0003  # phred 35
    df_pyclone.drop(columns="DP4", inplace=True)

    # if using existing pyclone output, extract mutation_id column to filter df_pyclone to only those mutation_ids and apply recomb filter
    if iteration > 0:
        usecols = ["mutation_id"]
        df_pyclone_old = pd.read_csv(
            f"{pat}/mut_dyn/pyclone/pyclone_out_{pat}_iter{iteration-1}_cluster_recomb.tsv",
            sep="\t",
            usecols=usecols,
        )
        mutation_ids = df_pyclone_old["mutation_id"].unique()
        df_pyclone = df_pyclone[df_pyclone.mutation_id.isin(mutation_ids)].copy()

    df_pyclone.to_csv(f"{pat}/mut_dyn/pyclone/{output_name}", sep="\t", index=False)


def run_stage(stage, pat, iteration, true_mut_df, args):
    if stage == "input_generated":
        gen_pyclone_df(pat, f"pyclone_in_{pat}_iter{iteration}.tsv", iteration)
        return None

    if stage == "fit_done":
        run_pyclone_fit(
            pat,
            args.max_cluster,
            args.restarts,
            iteration,
            num_threads=args.num_threads,
            seed=args.seed,
        )
        return None

    if stage == "results_written":
        write_pyclone_results(pat, iteration)
        return None

    if stage == "cluster_filtered":
        return cluster_filter(
            pat,
            args.std_thresh,
            args.fixed_thresh,
            args.neglegible_thresh,
            iteration,
            true_mut_df,
        )

    if stage == "recomb_filtered":
        return sliding_recomb_filter(pat, args.window_size, args.recomb_threshold, iteration)

    raise ValueError(f"Unknown stage: {stage}")


def iterative_pyclone(args):
    pat = args.pat
    params = {
        "max_cluster": args.max_cluster,
        "restarts": args.restarts,
        "num_threads": args.num_threads,
        "seed": args.seed,
        "std_thresh": args.std_thresh,
        "fixed_thresh": args.fixed_thresh,
        "neglegible_thresh": args.neglegible_thresh,
        "window_size": args.window_size,
        "recomb_threshold": args.recomb_threshold,
    }

    if args.fresh:
        cleanup_previous_outputs(pat)

    saved_state = load_state(pat)

    if args.resume and saved_state is not None:
        saved_params = saved_state.get("params", {})
        if not params_match(saved_params, params):
            raise ValueError(
                "Cannot resume: current arguments differ from saved run parameters. "
                "Use matching params or run with --fresh."
            )

        if saved_state.get("status") == "completed":
            print(f"Run already completed for {pat} at iteration {saved_state.get('iteration')}.")
            return

    if args.resume:
        iteration, start_stage = infer_resume_point(pat)
        print(f"Resuming from iteration {iteration}, stage '{start_stage}'.")
    else:
        iteration, start_stage = 0, "input_generated"

    true_mut_df = pd.read_csv(f"{pat}/mut_dyn/mut_evolution_dates.txt", sep="\t")
    true_mut_df["date"] = pd.to_datetime(true_mut_df["date"], format="%d-%m-%Y")
    true_mut_df["position"] = "pos_" + true_mut_df["position"].astype(str)

    while True:
        print(f"Iteration {iteration}")

        removed_clusters_count = None
        recomb_pos_count = None

        failure_stage = start_stage
        try:
            start_index = STAGES.index(start_stage)
            for stage in STAGES[start_index:]:
                failure_stage = stage
                print(f"Running stage: {stage}")
                t0 = time.time()
                stage_result = run_stage(stage, pat, iteration, true_mut_df, args)
                if not validate_stage_output(pat, iteration, stage):
                    raise RuntimeError(
                        f"Stage '{stage}' completed but output validation failed for iteration {iteration}."
                    )

                if stage == "cluster_filtered":
                    removed_clusters_count = len(stage_result)
                    print(f"Removed clusters: {stage_result}")
                elif stage == "recomb_filtered":
                    recomb_pos_count = len(stage_result)
                    print(f"Positions in recombination regions: {recomb_pos_count}")

                elapsed = time.time() - t0
                print(f"Execution time for {stage}: {elapsed:.2f} seconds")

                save_state(
                    pat,
                    params,
                    iteration,
                    stage,
                    status="running",
                    removed_clusters_count=removed_clusters_count,
                    recomb_pos_count=recomb_pos_count,
                )

            if removed_clusters_count is None or recomb_pos_count is None:
                latest_state = load_state(pat) or {}
                removed_clusters_count = latest_state.get("removed_clusters_count", removed_clusters_count)
                recomb_pos_count = latest_state.get("recomb_pos_count", recomb_pos_count)

            if removed_clusters_count == 0 and recomb_pos_count == 0:
                save_state(
                    pat,
                    params,
                    iteration,
                    "recomb_filtered",
                    status="completed",
                    removed_clusters_count=removed_clusters_count,
                    recomb_pos_count=recomb_pos_count,
                )
                print("No clusters or recombination positions removed; iteration complete.")
                break

            iteration += 1
            start_stage = "input_generated"

        except Exception as exc:
            save_state(
                pat,
                params,
                iteration,
                failure_stage,
                status="failed",
                last_error=str(exc),
                removed_clusters_count=removed_clusters_count,
                recomb_pos_count=recomb_pos_count,
            )
            raise


def main():
    parser = argparse.ArgumentParser(
        description="Run PyClone iteratively with resumable stage checkpoints."
    )
    parser.add_argument("pat", type=str, help="Patient name")
    parser.add_argument("--max_cluster", type=int, default=40, help="Maximum number of clusters")
    parser.add_argument("--restarts", type=int, default=50, help="Number of restarts for PyClone")
    parser.add_argument("--num_threads", type=int, default=20, help="Number of threads for PyClone fit")
    parser.add_argument("--seed", type=int, default=15, help="Random seed for PyClone fit")
    parser.add_argument(
        "--std_thresh", type=float, default=0.05, help="Standard deviation threshold for variability"
    )
    parser.add_argument(
        "--fixed_thresh", type=float, default=0.9, help="Threshold above which a cluster is considered fixed"
    )
    parser.add_argument(
        "--neglegible_thresh",
        type=float,
        default=0.05,
        help="Threshold below which a cluster is considered negligible",
    )
    parser.add_argument(
        "--window_size", type=int, default=500, help="Sliding window size for recombination filtering"
    )
    parser.add_argument(
        "--recomb_threshold",
        type=int,
        default=5,
        help="Threshold for mutations in a window to consider recombination",
    )

    mode_group = parser.add_mutually_exclusive_group()
    mode_group.add_argument("--resume", action="store_true", help="Resume from latest incomplete stage")
    mode_group.add_argument(
        "--fresh",
        action="store_true",
        help="Delete previous iterative outputs for this patient and start from iteration 0",
    )

    args = parser.parse_args()
    iterative_pyclone(args)


if __name__ == "__main__":
    main()
