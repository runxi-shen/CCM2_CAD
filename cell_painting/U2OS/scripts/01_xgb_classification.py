#!/usr/bin/env python3
"""
XGBoost Classification for CCM2_Val53Ile Variant Analysis

Runs binary classification between CCM2 (reference) and CCM2_Val53Ile (variant).
Each batch is processed separately with leave-one-plate-out CV.

Expected: 4 well pairs × 4 plates = 16 classifiers per batch
Total: 32 classifiers for Batch 18 + 19

Usage:
    pixi run classify
"""
import sys
sys.stdout.reconfigure(line_buffering=True)

import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import polars as pl
import xgboost as xgb
from sklearn.metrics import roc_auc_score

warnings.filterwarnings("ignore")

# ============================================================================
# Configuration
# ============================================================================
CONFIG = {
    "input_dir": Path(__file__).parent.parent / "data/raw/single_cell_profiles",
    "output_dir": Path(__file__).parent.parent / "data/interim/classification_results",
    "final_output_dir": Path(__file__).parent.parent / "data/processed",
    "reference_allele": "CCM2",
    "variant_allele": "CCM2_Val53Ile",
    "cc_threshold": 20,
}

GFP_INTENSITY_COLUMN = "Cells_Intensity_IntegratedIntensity_GFP"


# ============================================================================
# Utility Functions
# ============================================================================
def find_feat_cols(df):
    return [c for c in df.columns if not c.startswith("Metadata_")]

def find_meta_cols(df):
    return [c for c in df.columns if c.startswith("Metadata_")]

def get_gfp_features(df):
    """Extract GFP channel features only."""
    meta_cols = find_meta_cols(df)
    feat_cols = [c for c in find_feat_cols(df)
                 if "gfp" in c.lower() and "Brightfield" not in c and "TxControl" not in c]
    return df[meta_cols + feat_cols].copy()


# ============================================================================
# XGBoost Classifier (CPU, parallelized)
# ============================================================================
def run_classifier(df_train, df_test, feat_cols):
    """Train XGBoost and return AUROC, feature importance, and predictions."""
    x_train = df_train[feat_cols].values
    x_test = df_test[feat_cols].values
    y_train = df_train["Label"].values
    y_test = df_test["Label"].values

    num_pos = (y_train == 1).sum()
    num_neg = (y_train == 0).sum()

    if num_pos == 0 or num_neg == 0:
        return None, None, None, None

    scale_pos_weight = num_neg / num_pos

    model = xgb.XGBClassifier(
        objective="binary:logistic",
        n_estimators=150,
        tree_method="hist",
        learning_rate=0.05,
        scale_pos_weight=scale_pos_weight,
        n_jobs=4,  # Fixed number of cores
        verbosity=0,
    )
    model.fit(x_train, y_train)

    pred_score = model.predict_proba(x_test)[:, 1]
    auroc = roc_auc_score(y_test, pred_score)
    feat_importance = pd.Series(model.feature_importances_, index=feat_cols)

    # Classifier info
    info_0 = df_test[df_test["Label"] == 0].iloc[0]
    info_1 = df_test[df_test["Label"] == 1].iloc[0]
    class_id = f"{info_0['Metadata_Plate']}_{info_0['Metadata_well_position']}_{info_1['Metadata_well_position']}"

    classifier_info = {
        "Classifier_ID": class_id,
        "Plate": info_0["Metadata_Plate"],
        "trainsize_0": int((y_train == 0).sum()),
        "testsize_0": int((y_test == 0).sum()),
        "well_0": info_0["Metadata_well_position"],
        "allele_0": info_0["Metadata_gene_allele"],
        "trainsize_1": int((y_train == 1).sum()),
        "testsize_1": int((y_test == 1).sum()),
        "well_1": info_1["Metadata_well_position"],
        "allele_1": info_1["Metadata_gene_allele"],
        "AUROC": auroc,
    }

    # Predictions
    predictions = pd.DataFrame({
        "Classifier_ID": class_id,
        "Label": y_test,
        "Prediction": pred_score,
    })

    return auroc, feat_importance, classifier_info, predictions


# ============================================================================
# Standard Classification (Leave-one-plate-out CV)
# ============================================================================
def run_batch_classification(df, batch_name):
    """Run classification for a single batch."""
    print(f"\n  Processing {batch_name}...")

    # Extract GFP features
    df = get_gfp_features(df)
    feat_cols = [c for c in find_feat_cols(df) if c != "Label"]
    print(f"    GFP features: {len(feat_cols)}")

    # Separate ref and var
    df_ref = df[df["Metadata_node_type"] != "allele"].copy()
    df_ref["Label"] = 1
    df_var = df[df["Metadata_node_type"] == "allele"].copy()
    df_var["Label"] = 0

    ref_wells = df_ref["Metadata_well_position"].unique().tolist()
    var_wells = df_var["Metadata_well_position"].unique().tolist()
    plates = sorted(df["Metadata_Plate"].unique().tolist())

    print(f"    Ref wells: {ref_wells}, Var wells: {var_wells}")
    print(f"    Plates: {len(plates)}")

    # Create all ref-var well pair combinations
    well_pairs = [(r, v) for r in ref_wells for v in var_wells]
    print(f"    Well pairs: {len(well_pairs)}, Expected classifiers: {len(well_pairs) * len(plates)}")

    results = []
    all_feat_imp = []
    all_preds = []

    classifier_count = 0
    for ref_well, var_well in well_pairs:
        df_ref_pair = df_ref[df_ref["Metadata_well_position"] == ref_well]
        df_var_pair = df_var[df_var["Metadata_well_position"] == var_well]
        df_sampled = pd.concat([df_ref_pair, df_var_pair], ignore_index=True)

        for test_plate in plates:
            classifier_count += 1
            print(f"    Running classifier {classifier_count}...", end="", flush=True)
            df_test = df_sampled[df_sampled["Metadata_Plate"] == test_plate].reset_index(drop=True)
            df_train = df_sampled[df_sampled["Metadata_Plate"] != test_plate].reset_index(drop=True)

            if len(df_test) < 20 or len(df_train) < 20:
                print(" skipped (too few cells)", flush=True)
                continue

            auroc, feat_imp, info, preds = run_classifier(df_train, df_test, feat_cols)
            print(f" AUROC={auroc:.4f}" if auroc else " skipped", flush=True)

            if auroc is not None:
                results.append(info)
                all_feat_imp.append(feat_imp)
                all_preds.append(preds)

    print(f"    Completed: {len(results)} classifiers", flush=True)

    if not results:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

    df_result = pd.DataFrame(results)
    df_feat = pd.DataFrame(all_feat_imp)
    df_feat["Metadata_Feature_Type"] = "GFP"
    df_feat["Metadata_Control"] = False
    df_preds = pd.concat(all_preds, ignore_index=True)

    return df_feat, df_result, df_preds


# ============================================================================
# GFP-Corrected Classification
# ============================================================================
def find_gfp_range(ref_gfp, var_gfp, min_cells=20):
    """Find overlapping GFP intensity range."""
    if len(ref_gfp) == 0 or len(var_gfp) == 0:
        return None, None

    for low_q, high_q in [(0.25, 0.75), (0.2, 0.8), (0.15, 0.85), (0.1, 0.9)]:
        ref_low, ref_high = np.quantile(ref_gfp, [low_q, high_q])
        var_low, var_high = np.quantile(var_gfp, [low_q, high_q])

        gfp_min = max(ref_low, var_low)
        gfp_max = min(ref_high, var_high)

        if gfp_min >= gfp_max:
            continue

        ref_count = np.sum((ref_gfp >= gfp_min) & (ref_gfp <= gfp_max))
        var_count = np.sum((var_gfp >= gfp_min) & (var_gfp <= gfp_max))

        if ref_count >= min_cells and var_count >= min_cells:
            return gfp_min, gfp_max

    return None, None


def run_batch_classification_gfp_adj(df, batch_name):
    """Run GFP-corrected classification for a single batch."""
    print(f"\n  Processing {batch_name} (GFP-adj)...")

    # Check if GFP intensity column exists
    if GFP_INTENSITY_COLUMN not in df.columns:
        print(f"    WARNING: {GFP_INTENSITY_COLUMN} not found, skipping GFP-adj")
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

    # Extract GFP features but keep intensity column
    df_gfp = get_gfp_features(df)
    if GFP_INTENSITY_COLUMN not in df_gfp.columns:
        df_gfp[GFP_INTENSITY_COLUMN] = df[GFP_INTENSITY_COLUMN]

    feat_cols = [c for c in find_feat_cols(df_gfp) if c != "Label" and c != GFP_INTENSITY_COLUMN]

    # Separate ref and var
    df_ref = df_gfp[df_gfp["Metadata_node_type"] != "allele"].copy()
    df_ref["Label"] = 1
    df_var = df_gfp[df_gfp["Metadata_node_type"] == "allele"].copy()
    df_var["Label"] = 0

    ref_wells = df_ref["Metadata_well_position"].unique().tolist()
    var_wells = df_var["Metadata_well_position"].unique().tolist()
    plates = sorted(df_gfp["Metadata_Plate"].unique().tolist())

    well_pairs = [(r, v) for r in ref_wells for v in var_wells]
    print(f"    Well pairs: {len(well_pairs)}, Plates: {len(plates)}", flush=True)

    results = []
    all_feat_imp = []
    all_preds = []
    classifier_count = 0

    for ref_well, var_well in well_pairs:
        df_ref_pair = df_ref[df_ref["Metadata_well_position"] == ref_well]
        df_var_pair = df_var[df_var["Metadata_well_position"] == var_well]

        for test_plate in plates:
            classifier_count += 1
            print(f"    Running classifier {classifier_count}...", end="", flush=True)

            train_plates = [p for p in plates if p != test_plate]

            # Filter by GFP intensity per plate
            train_ref_list, train_var_list = [], []
            for plate in train_plates:
                df_ref_p = df_ref_pair[df_ref_pair["Metadata_Plate"] == plate]
                df_var_p = df_var_pair[df_var_pair["Metadata_Plate"] == plate]

                gfp_min, gfp_max = find_gfp_range(
                    df_ref_p[GFP_INTENSITY_COLUMN].values,
                    df_var_p[GFP_INTENSITY_COLUMN].values
                )
                if gfp_min is None:
                    continue

                df_ref_f = df_ref_p[(df_ref_p[GFP_INTENSITY_COLUMN] >= gfp_min) &
                                    (df_ref_p[GFP_INTENSITY_COLUMN] <= gfp_max)]
                df_var_f = df_var_p[(df_var_p[GFP_INTENSITY_COLUMN] >= gfp_min) &
                                    (df_var_p[GFP_INTENSITY_COLUMN] <= gfp_max)]
                train_ref_list.append(df_ref_f)
                train_var_list.append(df_var_f)

            # Test plate
            df_ref_test = df_ref_pair[df_ref_pair["Metadata_Plate"] == test_plate]
            df_var_test = df_var_pair[df_var_pair["Metadata_Plate"] == test_plate]
            gfp_min, gfp_max = find_gfp_range(
                df_ref_test[GFP_INTENSITY_COLUMN].values,
                df_var_test[GFP_INTENSITY_COLUMN].values
            )

            if gfp_min is None or not train_ref_list:
                print(" skipped (no GFP range)", flush=True)
                continue

            df_ref_test_f = df_ref_test[(df_ref_test[GFP_INTENSITY_COLUMN] >= gfp_min) &
                                        (df_ref_test[GFP_INTENSITY_COLUMN] <= gfp_max)]
            df_var_test_f = df_var_test[(df_var_test[GFP_INTENSITY_COLUMN] >= gfp_min) &
                                        (df_var_test[GFP_INTENSITY_COLUMN] <= gfp_max)]

            # Combine and drop GFP intensity column
            df_train = pd.concat(train_ref_list + train_var_list, ignore_index=True)
            df_test = pd.concat([df_ref_test_f, df_var_test_f], ignore_index=True)

            df_train = df_train.drop(columns=[GFP_INTENSITY_COLUMN])
            df_test = df_test.drop(columns=[GFP_INTENSITY_COLUMN])

            if len(df_test) < 20 or len(df_train) < 20:
                print(" skipped (too few cells)", flush=True)
                continue

            auroc, feat_imp, info, preds = run_classifier(df_train, df_test, feat_cols)
            print(f" AUROC={auroc:.4f}" if auroc else " skipped", flush=True)

            if auroc is not None:
                results.append(info)
                all_feat_imp.append(feat_imp)
                all_preds.append(preds)

    print(f"    Completed: {len(results)} classifiers", flush=True)

    if not results:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

    df_result = pd.DataFrame(results)
    df_feat = pd.DataFrame(all_feat_imp)
    df_feat["Metadata_Feature_Type"] = "GFP"
    df_feat["Metadata_Control"] = False
    df_preds = pd.concat(all_preds, ignore_index=True)

    return df_feat, df_result, df_preds


# ============================================================================
# Main
# ============================================================================
def main():
    print("=" * 60)
    print("CCM2_Val53Ile XGBoost Classification")
    print("=" * 60)

    CONFIG["output_dir"].mkdir(parents=True, exist_ok=True)
    CONFIG["final_output_dir"].mkdir(parents=True, exist_ok=True)

    # Load data
    print("\nLoading profiles...")
    df18 = pl.read_parquet(CONFIG["input_dir"] / "batch18_ccm2_profiles.parquet").to_pandas()
    df19 = pl.read_parquet(CONFIG["input_dir"] / "batch19_ccm2_profiles.parquet").to_pandas()
    print(f"  Batch 18: {len(df18)} cells, Batch 19: {len(df19)} cells")

    # =========================================================================
    # Standard GFP Classification
    # =========================================================================
    print("\n" + "-" * 40)
    print("Standard GFP Classification")
    print("-" * 40)

    feat18, result18, preds18 = run_batch_classification(df18, "Batch 18")
    feat19, result19, preds19 = run_batch_classification(df19, "Batch 19")

    # Combine results
    df_feat = pd.concat([feat18, feat19], ignore_index=True)
    df_result = pd.concat([result18, result19], ignore_index=True)
    df_preds = pd.concat([preds18, preds19], ignore_index=True)

    if not df_result.empty:
        df_feat.to_csv(CONFIG["output_dir"] / "feature_importance.csv", index=False)
        df_result.to_csv(CONFIG["output_dir"] / "classifier_info.csv", index=False)
        df_preds.to_parquet(CONFIG["output_dir"] / "predictions.parquet", index=False)

        auroc_mean = df_result["AUROC"].mean()
        auroc_std = df_result["AUROC"].std()
        print(f"\n  Overall AUROC: {auroc_mean:.4f} ± {auroc_std:.4f}")
        print(f"  Batch 18 AUROC: {result18['AUROC'].mean():.4f} ± {result18['AUROC'].std():.4f}")
        print(f"  Batch 19 AUROC: {result19['AUROC'].mean():.4f} ± {result19['AUROC'].std():.4f}")

    # =========================================================================
    # GFP-Corrected Classification (using full GFP profiles)
    # =========================================================================
    print("\n" + "-" * 40)
    print("GFP-Corrected Classification")
    print("-" * 40)

    # Load full GFP profiles with IntegratedIntensity
    full_gfp_18_path = CONFIG["input_dir"] / "batch18_ccm2_profiles_full_gfp.parquet"
    full_gfp_19_path = CONFIG["input_dir"] / "batch19_ccm2_profiles_full_gfp.parquet"

    if full_gfp_18_path.exists() and full_gfp_19_path.exists():
        print("\nLoading full GFP profiles (with IntegratedIntensity)...")
        df18_full = pl.read_parquet(full_gfp_18_path).to_pandas()
        df19_full = pl.read_parquet(full_gfp_19_path).to_pandas()
        print(f"  Batch 18: {len(df18_full)} cells, Batch 19: {len(df19_full)} cells")

        feat18_adj, result18_adj, preds18_adj = run_batch_classification_gfp_adj(df18_full, "Batch 18")
        feat19_adj, result19_adj, preds19_adj = run_batch_classification_gfp_adj(df19_full, "Batch 19")
    else:
        print("\n  Full GFP profiles not found. Run 'pixi run extract-data' first.")
        feat18_adj, result18_adj, preds18_adj = pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
        feat19_adj, result19_adj, preds19_adj = pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

    df_feat_adj = pd.concat([feat18_adj, feat19_adj], ignore_index=True)
    df_result_adj = pd.concat([result18_adj, result19_adj], ignore_index=True)
    df_preds_adj = pd.concat([preds18_adj, preds19_adj], ignore_index=True)

    if not df_result_adj.empty:
        df_feat_adj.to_csv(CONFIG["output_dir"] / "feature_importance_gfp_adj.csv", index=False)
        df_result_adj.to_csv(CONFIG["output_dir"] / "classifier_info_gfp_adj.csv", index=False)
        df_preds_adj.to_parquet(CONFIG["output_dir"] / "predictions_gfp_adj.parquet", index=False)

        auroc_adj_mean = df_result_adj["AUROC"].mean()
        auroc_adj_std = df_result_adj["AUROC"].std()
        print(f"\n  Overall AUROC: {auroc_adj_mean:.4f} ± {auroc_adj_std:.4f}")
        if not result18_adj.empty:
            print(f"  Batch 18 AUROC: {result18_adj['AUROC'].mean():.4f} ± {result18_adj['AUROC'].std():.4f}")
        if not result19_adj.empty:
            print(f"  Batch 19 AUROC: {result19_adj['AUROC'].mean():.4f} ± {result19_adj['AUROC'].std():.4f}")
    else:
        auroc_adj_mean = np.nan
        auroc_adj_std = np.nan
        print("\n  GFP-corrected classification skipped (intensity column not available)")

    # =========================================================================
    # Summary
    # =========================================================================
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)

    summary = pd.DataFrame({
        "gene_allele": [CONFIG["variant_allele"]],
        "AUROC_standard_mean": [df_result["AUROC"].mean() if not df_result.empty else np.nan],
        "AUROC_standard_std": [df_result["AUROC"].std() if not df_result.empty else np.nan],
        "AUROC_gfp_adj_mean": [df_result_adj["AUROC"].mean() if not df_result_adj.empty else np.nan],
        "AUROC_gfp_adj_std": [df_result_adj["AUROC"].std() if not df_result_adj.empty else np.nan],
        "n_classifiers_standard": [len(df_result)],
        "n_classifiers_gfp_adj": [len(df_result_adj)],
    })
    summary.to_csv(CONFIG["final_output_dir"] / "ccm2_val53ile_phenotype_score.csv", index=False)

    print(f"\nVariant: {CONFIG['variant_allele']}")
    if not df_result.empty:
        print(f"Standard GFP: {df_result['AUROC'].mean():.4f} ± {df_result['AUROC'].std():.4f} ({len(df_result)} classifiers)")
    if not df_result_adj.empty:
        print(f"GFP-Adjusted: {df_result_adj['AUROC'].mean():.4f} ± {df_result_adj['AUROC'].std():.4f} ({len(df_result_adj)} classifiers)")
    else:
        print("GFP-Adjusted: N/A (intensity column not in feature-selected data)")
    print(f"\nReference targets: ~0.976 for both modes")
    print(f"\nResults saved to: {CONFIG['output_dir']}")
    print("\nDone!")


if __name__ == "__main__":
    main()
