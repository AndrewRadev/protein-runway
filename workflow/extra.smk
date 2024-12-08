rule report_runtime_plot:
    input:
        tsv_files=expand("04_extra/benchmarks/{rule_name}/{protein_name}.tsv",
                         rule_name=["amsm", "build_nmd", "build_pdbs",
                                    "nmd_trajectory", "geostas", "segment_intermediates",
                                    "merizo_clustering", "chainsaw_clustering"],
                         protein_name=protein_names)
    output:
        aggregated_csv=report("04_extra/benchmarks/aggregated_runtime.csv"),
        runtime_plot=report("04_extra/benchmarks/runtime_plot.png")
    run:
        import pandas as pd

        # Using Agg to avoid problems with threading:
        # https://matplotlib.org/stable/users/explain/figure/backends.html
        #
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        from pathlib import Path

        def aggregate_benchmarks(tsv_files, output_file):
            benchmarks = []
            for file in tsv_files:
                data = pd.read_csv(file, sep="\t")
                rule_name = Path(file).parts[-2]  # Extract rule name
                sample_name = Path(file).stem  # Extract protein name
                runtime = data["cpu_time"].sum() if "cpu_time" in data.columns else 0
                benchmarks.append({"rule": rule_name, "sample": sample_name, "runtime": runtime})
            df = pd.DataFrame(benchmarks)
            df.to_csv(output_file, index=False)
            print(f"Aggregated data saved to {output_file}")

        def plot_sample_runtimes(data_file, output_plot):
            df = pd.read_csv(data_file)
            plt.figure(figsize=(12, 8))
            for rule in df["rule"].unique():
                subset = df[df["rule"] == rule]
                plt.scatter(subset["runtime"], [rule] * len(subset), label=rule, alpha=0.7)
            plt.xlabel("Runtime (seconds)")
            plt.ylabel("Rule")
            plt.title("Runtime of Snakemake Rules by Sample")
            plt.tight_layout()
            plt.savefig(output_plot)
            plt.close()
            print(f"Runtime plot saved to {output_plot}")

        aggregate_benchmarks(input.tsv_files, output.aggregated_csv)
        plot_sample_runtimes(output.aggregated_csv, output.runtime_plot)
