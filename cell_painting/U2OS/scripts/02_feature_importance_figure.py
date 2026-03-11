"""
Generate feature importance figure for CCM2_Val53Ile cell painting analysis.

Creates a horizontal barplot of top N XGBoost features used to distinguish
V53I variant from CCM2 reference gene, averaged across all classifiers.

Format matches: 6_feature_importance_analyses.ipynb

Usage:
    pixi run feature-importance
"""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

# ============================================================================
# Configuration
# ============================================================================

DATA_DIR = Path(__file__).parent.parent / "data"
OUTPUT_DIR = DATA_DIR / "processed" / "manuscript_figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Feature importance settings
TOP_N_FEATURES = 10


# ============================================================================
# Figure: Feature Importance Barplot
# ============================================================================

def plot_top_n_important_feat(df: pd.DataFrame, n: int = 10, ax=None, title: str = ""):
    """
    Plot top N most important features as horizontal barplot.

    Matches the format from 6_feature_importance_analyses.ipynb.

    Parameters
    ----------
    df : pd.DataFrame
        Feature importance dataframe (rows = classifiers, cols = features)
    n : int
        Number of top features to display
    ax : matplotlib.axes.Axes, optional
        Axes to plot on. If None, creates new figure.
    title : str
        Title suffix for the plot

    Returns
    -------
    pd.DataFrame
        Sorted feature statistics (mean, std) for top N features
    """
    # Select only GFP feature columns (exclude metadata)
    gfp_cols = [col for col in df.columns if "GFP" in col]

    # Calculate mean and std for each feature across classifiers
    column_stats = df[gfp_cols].agg(['mean', 'std']).T

    # Sort by mean importance (descending) and take top N
    sorted_columns = column_stats.sort_values(by="mean", ascending=False).head(n)

    # Create figure if not provided
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 6))

    # Plot horizontal barplot
    ax.barh(
        y=sorted_columns.index,
        width=sorted_columns["mean"],
        color='skyblue'
    )

    # Format axes
    ax.set_yticks(range(len(sorted_columns.index)))
    ax.set_yticklabels(sorted_columns.index, rotation=0, fontsize=12)
    ax.set_xlabel("Mean Feature Importance", fontsize=12)
    ax.set_title(title if title else "Feature Importance", fontsize=14)

    # Set x-axis ticks
    ax.set_xticks([0, 0.05, 0.10, 0.15, 0.20])
    ax.set_xlim(0, 0.22)
    ax.tick_params(axis='x', labelsize=11)

    # Add reference line at 0.01 with label for legend
    ax.axvline(x=0.01, color='r', linestyle='--', label='Threshold (0.01)')
    ax.legend(loc='lower right', fontsize=11)

    # Invert y-axis so highest importance is at top
    ax.invert_yaxis()

    return sorted_columns


def create_feature_importance_figure():
    """Create barplot of top 10 features averaged across all classifiers."""
    print("Creating feature importance figure...")

    # Load feature importance (rows = classifiers, cols = features)
    feat_imp_path = DATA_DIR / "interim" / "classification_results" / "feature_importance.csv"
    feat_imp = pd.read_csv(feat_imp_path)

    print(f"  Loaded: {feat_imp.shape[0]} classifiers, {feat_imp.shape[1]} columns")

    # Create figure
    fig, ax = plt.subplots(figsize=(10, 6))

    # Plot top N features
    top_features = plot_top_n_important_feat(
        feat_imp,
        n=TOP_N_FEATURES,
        ax=ax,
        title="Mean XGBoost Feature Importance\n(CCM2_Val53Ile vs CCM2)"
    )

    print(f"  Top feature: {top_features.index[0]} ({top_features['mean'].iloc[0]:.4f})")

    plt.tight_layout()

    # Save figure
    output_path = OUTPUT_DIR / "FigA_feat_importance.png"
    fig.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)

    print(f"  Saved: {output_path}")

    return top_features


# ============================================================================
# Main
# ============================================================================

def main():
    print("=" * 60)
    print("Generating feature importance figure")
    print("=" * 60)

    top_features = create_feature_importance_figure()

    print("\n" + "=" * 60)
    print("Done!")
    print("=" * 60)


if __name__ == "__main__":
    main()
