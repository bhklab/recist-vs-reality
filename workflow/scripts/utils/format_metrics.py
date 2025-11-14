import pandas as pd
import numpy as np
import click

@click.command()
@click.argument('input_csv', type=click.STRING)
def format_metrics(input_csv:str) -> None:
    """Print out summary statistics from MedSAM2-RECIST evaluation metrics CSV file."""
    # Load in the metrics CSV file 
    metrics_df = pd.read_csv(input_csv, index_col=0)
    print(f"Number of samples: {metrics_df.shape[0]}")
    print(f"Model name: {metrics_df['model_name'].unique()[0]}")

    # Drop the filename column
    if 'filename' in metrics_df.columns:
        metrics_df = metrics_df.drop(columns=['filename'])
    if 'model_name' in metrics_df.columns:
        metrics_df = metrics_df.drop(columns=['model_name'])
    
    print(f"Metrics: {metrics_df.columns.to_list()}")
    for _, metric_vals in metrics_df.items():
        print(f"{np.round(metric_vals.mean(), 4)} ({np.round(metric_vals.min(), 4)}-{np.round(metric_vals.max(), 4)})")


if __name__ == "__main__":
    format_metrics()