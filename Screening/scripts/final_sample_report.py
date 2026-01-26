import pandas as pd

samples = snakemake.params.samples
all_metrics = []

for sample in samples:
    m_path = f"results/{sample}/reporting/{sample}_metrics.csv"
    c_path = f"results/{sample}/reporting/{sample}_coverage_metrics.txt"

    # Load metrics
    df_metrics = pd.read_csv(m_path)

    # Load summarized coverage values (2-row TSV)
    try:
        df_cov = pd.read_csv(c_path, sep="\t")
        mean = df_cov.iloc[0]["mean_depth"]
        even = df_cov.iloc[0]["evenness"]
    except Exception as e:
        print(f"Warning: Failed to parse coverage for {sample}: {e}")
        mean = even = 0

    df_metrics["mean_depth"] = mean
    df_metrics["evenness"] = even
    all_metrics.append(df_metrics)

# Merge and write output
pd.concat(all_metrics).to_csv(snakemake.output[0], index=False)
