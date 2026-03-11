"""
Generate cell crop figures for CCM2_Val53Ile manuscript.

Reproduces the exact visualization from
/home/shenrunx/igvf/varchamp/2025_Fang_CCM2_Imaging/notebooks/6_visualize_cell_crop.ipynb

Two selection modes:
1. Random selection (default): Random sampling of cells
2. Feature-based selection: Select cells based on feature percentiles
   - Lower mean allele: 25-50% percentile
   - Higher mean allele: 50-75% percentile

Usage:
    pixi run python scripts/03_cell_crop_figures.py           # Random selection
    pixi run python scripts/03_cell_crop_figures.py --mode feature  # Feature-based
    pixi run python scripts/03_cell_crop_figures.py --mode both     # Both figures
"""

import argparse
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd
import tifffile

# ============================================================================
# Configuration
# ============================================================================

DATA_DIR = Path(__file__).parent.parent / "data"
OUTPUT_DIR = DATA_DIR / "processed" / "manuscript_figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Image paths
TIFF_IMGS_DIR = "/home/shenrunx/igvf/varchamp/2025_varchamp_snakemake/1.image_preprocess_qc/inputs/cpg_imgs"

# Cell crop visualization parameters
# For feature-based selection: 50 cells per allele (10 cols x 5 rows)
N_CELLS = 50
N_COLS = 10
N_ROWS = 5

# For reproduce mode: 5 cells per allele (5 cols x 1 row each, stacked vertically)
N_CELLS_REPRODUCE = 5
N_COLS_REPRODUCE = 5

# Scale bar parameters
# Phenix 20x objective, 2×2 binning: 0.598 µm/pixel
PIXEL_SIZE_UM = 0.598
SCALE_BAR_LENGTH_UM = 20  # 20 µm scale bar
ALLELE_LIST = ["CCM2_Val53Ile", "CCM2"]  # Skip CCM2_Phe270Leu as requested
CROP_SIZE = 128
PERCENTILE_LOW = 1
PERCENTILE_HIGH = 99

# Quality filter parameters
# NOTE: The original notebook used raw intensity values with ranges [40,60] and [20,60].
# Our local data is MAD-normalized, so we use percentile-based filtering instead.
MIN_EDGE_DIST = 50
INTENSITY_FEATURE = "Cells_Intensity_IntegratedIntensity_GFP"
# Use percentile ranges to approximate the notebook's filtering (roughly middle 50% of cells)
INTENSITY_PERCENTILE_LOW = 25
INTENSITY_PERCENTILE_HIGH = 75

# Feature-based selection
SELECTION_FEATURE = "Cytoplasm_Intensity_MinIntensityEdge_GFP"

# Exact cell IDs from the reference notebook (6_visualize_cell_crop.ipynb)
# These reproduce the exact visualization from the VarChAMP paper
# Reordered for manuscript figure:
# - CCM2: last cell moved to first, then first 4
# - CCM2_Val53Ile: cells 6-10 (second half)
NOTEBOOK_CELL_IDS = {
    "CCM2_Val53Ile": [
        # Using cells 6-10 from notebook, 4th moved to first
        "2025_06_02_B18A8A10R1_P2T1_C23_636_51",   # was 4th, moved to first
        "2025_06_02_B18A8A10R1_P2T1_C23_637_151",
        "2025_06_02_B18A8A10R1_P2T1_C23_632_140",
        "2025_06_02_B18A8A10R1_P2T1_C23_635_96",
        "2025_06_02_B18A8A10R1_P2T1_C23_637_139",
    ],
    "CCM2": [
        # Reordered: cell 7 first, then first 4
        "2025_06_02_B18A8A10R1_P2T1_A19_164_46",
        "2025_06_02_B18A8A10R1_P2T1_A19_166_43",
        "2025_06_02_B18A8A10R1_P2T1_A19_166_53",
        "2025_06_02_B18A8A10R1_P2T1_A19_166_34",
        "2025_06_02_B18A8A10R1_P2T1_A19_165_33",
    ],
}

# Letter dict for well position conversion
LETTER_DICT = {chr(i): f"{i-64:02d}" for i in range(65, 81)}  # A=01, B=02, ..., P=16

# Channel mapping
CHANNEL_DICT = {"DAPI": "1", "GFP": "2", "AGP": "3", "Mito": "4"}

# Plate directory mapping for Batch 18-19
PLATE_DICT = {
    "B18A8A10R1_P2": {
        "T1": "2025_06_02_B18A8A10R1_P2T1__2025-06-02T08_19_26-Measurement1",
        "T2": "2025_06_02_B18A8A10R1_P2T2__2025-06-02T09_29_48-Measurement1",
        "T3": "2025_06_02_B18A8A10R1_P2T3__2025-06-02T11_08_44-Measurement1",
        "T4": "2025_06_02_B18A8A10R1_P2T4__2025-06-02T12_18_34-Measurement1",
    },
    "B19A8A10R1_P2": {
        "T1": "2025_06_03_B19A8A10R1_P2T1__2025-06-03T13_27_13-Measurement1",
        "T2": "2025_06_04_B19A8A10R1_P2T2__2025-06-04T08_17_36-Measurement2",
        "T3": "2025_06_04_B19A8A10R1_P2T3__2025-06-04T09_30_37-Measurement1",
        "T4": "2025_06_04_B19A8A10R1_P2T4__2025-06-04T11_49_06-Measurement1",
    },
}

# Batch directory mapping
BATCH_DICT = {
    "B18A8A10R1": "2025_06_10_Batch_18",
    "B19A8A10R1": "2025_06_10_Batch_19",
}


# ============================================================================
# Image Utilities - Exact same as VarChAMP
# ============================================================================

def channel_to_cmap(channel: str):
    """Return colormap for channel - exact same as VarChAMP img_utils."""
    if channel == "GFP":
        return mpl.colors.LinearSegmentedColormap.from_list("gfp_cmap", ["#000", "#65fe08"])
    elif channel == "DAPI":
        return mpl.colors.LinearSegmentedColormap.from_list("dapi_cmap", ["#000", "#0000FF"])
    elif channel == "Mito":
        return mpl.colors.LinearSegmentedColormap.from_list("mito_cmap", ["#000", "#FF0000"])
    elif channel == "AGP":
        return mpl.colors.LinearSegmentedColormap.from_list("agp_cmap", ["#000", "#FFFF00"])
    return "gray"


def get_image_path(metadata_plate: str, well: str, site: int, channel: str) -> Path:
    """
    Construct image file path from metadata.

    Parameters:
    -----------
    metadata_plate : str
        Plate identifier (e.g., '2025_06_02_B18A8A10R1_P2T1')
    well : str
        Well position (e.g., 'A19')
    site : int
        Site/FOV number (1-9)
    channel : str
        Channel name ('DAPI', 'GFP', 'AGP', 'Mito')
    """
    # Parse plate info: 2025_06_02_B18A8A10R1_P2T1 -> B18A8A10R1, P2, T1
    parts = metadata_plate.split("_")
    plate_prefix = parts[3]  # B18A8A10R1
    plate_num = parts[4][:2]  # P2
    timepoint = parts[4][2:]  # T1

    plate_map_name = f"{plate_prefix}_{plate_num}"  # B18A8A10R1_P2
    batch_dir = BATCH_DICT[plate_prefix]

    # Get timepoint directory
    timepoint_dir = PLATE_DICT[plate_map_name][timepoint]

    # Convert well to row/col format
    row_letter = well[0]
    col_num = well[1:].zfill(2)
    row_num = LETTER_DICT[row_letter]

    # Get channel number
    channel_num = CHANNEL_DICT[channel]

    # Construct filename
    filename = f"r{row_num}c{col_num}f{site:02d}p01-ch{channel_num}sk1fk1fl1.tiff"

    return Path(TIFF_IMGS_DIR) / batch_dir / "images" / timepoint_dir / "Images" / filename


def load_channel_image(metadata_plate: str, well: str, site: int, channel: str) -> np.ndarray:
    """Load a single channel TIFF image."""
    path = get_image_path(metadata_plate, well, site, channel)
    if not path.exists():
        raise FileNotFoundError(f"Image not found: {path}")
    return tifffile.imread(path)


def compute_distance_to_edge(x: float, y: float, img_width: int = 2160, img_height: int = 2160) -> float:
    """Compute distance to nearest image edge."""
    return min(x, img_width - x, y, img_height - y)


def extract_fixed_crop(img: np.ndarray, center_x: float, center_y: float, crop_size: int = 128) -> np.ndarray:
    """
    Extract fixed-size crop centered on given coordinates.

    This matches the 'fixed' method in VarChAMP cell_cropper.
    """
    half = crop_size // 2
    cx, cy = int(center_x), int(center_y)
    h, w = img.shape[:2]

    # Calculate crop boundaries
    x_min = cx - half
    x_max = cx + half
    y_min = cy - half
    y_max = cy + half

    # Calculate padding needed
    pad_left = max(0, -x_min)
    pad_right = max(0, x_max - w)
    pad_top = max(0, -y_min)
    pad_bottom = max(0, y_max - h)

    # Clip coordinates to image boundaries
    x_min_clip = max(0, x_min)
    x_max_clip = min(w, x_max)
    y_min_clip = max(0, y_min)
    y_max_clip = min(h, y_max)

    # Extract the valid region
    crop = img[y_min_clip:y_max_clip, x_min_clip:x_max_clip]

    # Add padding if needed
    if pad_left > 0 or pad_right > 0 or pad_top > 0 or pad_bottom > 0:
        crop = np.pad(
            crop,
            ((pad_top, pad_bottom), (pad_left, pad_right)),
            mode='constant',
            constant_values=0
        )

    return crop


def viz_cell_single_channel(
    crop: np.ndarray,
    channel: str,
    ax: plt.Axes,
    percentile_low: float = 1.0,
    percentile_high: float = 99.5,
    axis_off: bool = True,
    title: str = None,
):
    """
    Visualize a single channel cell crop.

    Exact same as VarChAMP cell_visualizer.viz_cell_single_channel.
    """
    # Calculate display bounds from raw data
    display_vmin = np.percentile(crop, percentile_low)
    display_vmax = np.percentile(crop, percentile_high)

    # Get colormap for this channel
    cmap = channel_to_cmap(channel)

    # Display RAW data with calculated bounds
    ax.imshow(crop, cmap=cmap, vmin=display_vmin, vmax=display_vmax)

    if axis_off:
        ax.axis('off')

    if title:
        ax.set_title(title, fontsize=10)


def add_scale_bar(
    ax: plt.Axes,
    pixel_size_um: float = PIXEL_SIZE_UM,
    image_width_pixels: int = CROP_SIZE,
    scale_bar_length_um: float = SCALE_BAR_LENGTH_UM,
    location: str = 'lower right',
    color: str = 'white',
    fontsize: int = 8,
    bar_height_pixels: float = 2.0,
    padding_fraction: float = 0.05,
    label_offset_pixels: float = 3.0
) -> None:
    """Add a scale bar to a microscopy image axes."""
    # Get image dimensions
    try:
        images = ax.get_images()
        if len(images) > 0:
            img_array = images[0].get_array()
            image_height_pixels, image_width_pixels = img_array.shape[:2]
    except (IndexError, AttributeError):
        image_height_pixels = image_width_pixels

    # Calculate scale bar length in pixels
    scale_bar_pixels = scale_bar_length_um / pixel_size_um

    # Calculate padding
    padding_x = image_width_pixels * padding_fraction
    padding_y = image_height_pixels * padding_fraction

    # Position calculation
    if 'right' in location.lower():
        bar_x_end = image_width_pixels - padding_x
        bar_x_start = bar_x_end - scale_bar_pixels
        text_x = (bar_x_start + bar_x_end) / 2
    else:
        bar_x_start = padding_x
        bar_x_end = bar_x_start + scale_bar_pixels
        text_x = (bar_x_start + bar_x_end) / 2

    if 'lower' in location.lower():
        bar_y = image_height_pixels - padding_y - bar_height_pixels
        text_y = bar_y - label_offset_pixels
    else:
        bar_y = padding_y
        text_y = bar_y - label_offset_pixels

    # Create scale bar rectangle
    scale_bar_rect = Rectangle(
        (bar_x_start, bar_y),
        scale_bar_pixels,
        bar_height_pixels,
        facecolor=color,
        edgecolor='none',
        transform=ax.transData,
        zorder=1000
    )
    ax.add_patch(scale_bar_rect)

    # Add text label
    label_text = f"{int(scale_bar_length_um)} µm"
    ax.text(
        text_x, text_y, label_text,
        color=color, fontsize=fontsize,
        ha='center', va='bottom',
        transform=ax.transData,
        zorder=1001, weight='bold'
    )


# ============================================================================
# Cell Selection Functions
# ============================================================================

def select_notebook_cells(var_profiles: dict) -> dict:
    """
    Select the exact cells used in the reference notebook (6_visualize_cell_crop.ipynb).

    Uses NOTEBOOK_CELL_IDS to find matching cells in local data.
    """
    selected_indices = {}

    for allele in ALLELE_LIST:
        df = var_profiles[allele]
        indices = []

        for cell_id in NOTEBOOK_CELL_IDS[allele]:
            # Parse the notebook CellID format: Plate_Well_ImageNumber_ObjectNumber
            parts = cell_id.split("_")
            plate = "_".join(parts[:5])  # 2025_06_02_B18A8A10R1_P2T1
            well = parts[5]               # C23 or A19
            img_num = int(parts[6])       # 634
            obj_num = int(parts[7])       # 94

            # Find matching cell in our data
            mask = (
                (df["Metadata_Plate"] == plate) &
                (df["Metadata_well_position"] == well) &
                (df["Metadata_ImageNumber"] == img_num) &
                (df["Metadata_ObjectNumber"] == obj_num)
            )

            matches = df[mask]
            if len(matches) > 0:
                idx = matches.index[0]
                indices.append(df.index.get_loc(idx))
                print(f"  Found: {cell_id}")
            else:
                print(f"  NOT FOUND: {cell_id}")

        selected_indices[allele] = indices
        print(f"  {allele}: {len(indices)}/{len(NOTEBOOK_CELL_IDS[allele])} cells found")

    return selected_indices


def select_feature_based_cells(var_profiles: dict, feature: str, n_cells: int, seed: int = 42) -> tuple:
    """
    Select cells based on feature percentiles.

    - Allele with lower mean: select from 10-50% percentile
    - Allele with higher mean: select from 50-90% percentile

    Returns:
        selected_indices: dict of indices for each allele
        selection_info: dict with feature means and percentile ranges
    """
    np.random.seed(seed)

    # Calculate mean feature value for each allele
    means = {allele: var_profiles[allele][feature].mean() for allele in ALLELE_LIST}

    # Determine which allele has higher/lower mean
    sorted_alleles = sorted(means.keys(), key=lambda x: means[x])
    lower_allele = sorted_alleles[0]
    higher_allele = sorted_alleles[1]

    print(f"\n  Feature: {feature}")
    print(f"  {lower_allele} mean: {means[lower_allele]:.4f} (lower)")
    print(f"  {higher_allele} mean: {means[higher_allele]:.4f} (higher)")

    selected_indices = {}
    selection_info = {
        "feature": feature,
        "means": means,
        "lower_allele": lower_allele,
        "higher_allele": higher_allele,
    }

    for allele in ALLELE_LIST:
        df = var_profiles[allele]
        feat_values = df[feature]

        if allele == lower_allele:
            # Select from 10-50% percentile
            p10 = feat_values.quantile(0.10)
            p50 = feat_values.quantile(0.50)
            mask = (feat_values >= p10) & (feat_values <= p50)
            pct_range = "10-50%"
        else:
            # Select from 50-90% percentile
            p50 = feat_values.quantile(0.50)
            p90 = feat_values.quantile(0.90)
            mask = (feat_values >= p50) & (feat_values <= p90)
            pct_range = "50-90%"

        eligible_indices = df[mask].index.tolist()
        n_available = len(eligible_indices)

        if n_available < n_cells:
            print(f"  Warning: Only {n_available} cells in {pct_range} for {allele}")

        # Sample from eligible cells (subsample to n_cells if too many)
        sampled = np.random.choice(
            eligible_indices,
            size=min(n_cells, n_available),
            replace=False
        )
        # Convert to positional indices
        selected_indices[allele] = [df.index.get_loc(idx) for idx in sampled]
        selection_info[allele] = {
            "pct_range": pct_range,
            "n_eligible": n_available,
            "n_selected": len(selected_indices[allele]),
        }
        print(f"  {allele}: {len(selected_indices[allele])} cells from {pct_range} percentile ({n_available} eligible)")

    return selected_indices, selection_info


# ============================================================================
# Figure Generation
# ============================================================================

def load_crops(var_profiles: dict, selected_indices: dict) -> dict:
    """Load cell crops for selected indices."""
    all_crops = {}
    for allele in ALLELE_LIST:
        all_crops[allele] = []
        indices = selected_indices[allele]
        for cell_idx, i in enumerate(indices[:N_CELLS]):
            cell_row = var_profiles[allele].iloc[i]
            print(f"  Loading cell {cell_idx+1}/{len(indices)} for {allele}...", end="\r")

            # Load GFP image
            img = load_channel_image(
                cell_row["Metadata_Plate"],
                cell_row["Metadata_well_position"],
                int(cell_row["site"]),
                "GFP"
            )

            # Extract crop centered on nucleus (fixed method)
            crop = extract_fixed_crop(
                img,
                cell_row["Nuclei_Location_Center_X"],
                cell_row["Nuclei_Location_Center_Y"],
                CROP_SIZE
            )
            all_crops[allele].append(crop)
        print()
    return all_crops


def create_figure_small(
    all_crops: dict,
    output_path: Path,
    suptitle: str = None,
    allele_subtitles: dict = None,
):
    """Create small cell crop figure (5 cells per allele, stacked vertically)."""
    # Create figure - stacked layout (CCM2 on top, CCM2_Val53Ile on bottom)
    # 2 rows (one per allele) x 5 cols
    fig, axes = plt.subplots(
        2, N_COLS_REPRODUCE,
        figsize=(N_COLS_REPRODUCE * 1.8, 2 * 2.0),
        gridspec_kw={'wspace': 0.02, 'hspace': 0.15}
    )

    # Allele order: CCM2 on top (row 0), CCM2_Val53Ile on bottom (row 1)
    allele_order = ["CCM2", "CCM2_Val53Ile"]

    # First, hide all axes
    for ax_row in axes:
        for ax in ax_row:
            ax.axis('off')

    # Plot crops for each allele in its row
    for allele_idx, allele in enumerate(allele_order):
        crops = all_crops[allele]

        for cell_idx, crop in enumerate(crops[:N_CELLS_REPRODUCE]):
            ax = axes[allele_idx, cell_idx]
            viz_cell_single_channel(
                crop,
                channel="GFP",
                ax=ax,
                axis_off=True,
                percentile_low=PERCENTILE_LOW,
                percentile_high=PERCENTILE_HIGH,
            )
            # Add scale bar to first cell of each row
            if cell_idx == 0:
                add_scale_bar(ax)

    # Add row labels on the left using text annotation
    for allele_idx, allele in enumerate(allele_order):
        if allele_subtitles and allele in allele_subtitles:
            label = f"{allele}\n{allele_subtitles[allele]}"
        else:
            label = allele
        axes[allele_idx, 0].text(-0.15, 0.5, label, fontsize=12, fontweight='bold',
                                  transform=axes[allele_idx, 0].transAxes,
                                  rotation=0, ha='right', va='center')

    plt.tight_layout(rect=[0.12, 0, 1, 1])  # Leave space on left for labels

    # Add suptitle if provided
    if suptitle:
        fig.suptitle(suptitle, fontsize=12, fontweight='bold', y=1.02, x=0.5, ha='center')

    # Save figure
    fig.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"\nSaved: {output_path}")


def create_figure(
    all_crops: dict,
    output_path: Path,
    suptitle: str = None,
    allele_subtitles: dict = None,
):
    """Create cell crop figure with side-by-side layout."""
    # Create figure - side by side layout (CCM2 left, CCM2_Val53Ile right)
    fig, axes = plt.subplots(
        N_ROWS, N_COLS * len(ALLELE_LIST),
        figsize=(N_COLS * len(ALLELE_LIST) * 1.2, N_ROWS * 1.5),
        gridspec_kw={'wspace': 0.02, 'hspace': 0.02}
    )

    # Reorder alleles: CCM2 on left, CCM2_Val53Ile on right
    allele_order = ["CCM2", "CCM2_Val53Ile"]

    # First, hide all axes (will show only those with data)
    for ax_row in axes:
        for ax in ax_row:
            ax.axis('off')

    # Plot crops for each allele side by side
    for allele_idx, allele in enumerate(allele_order):
        crops = all_crops[allele]
        col_offset = allele_idx * N_COLS

        for cell_idx, crop in enumerate(crops[:N_CELLS]):
            row = cell_idx // N_COLS
            col = col_offset + (cell_idx % N_COLS)

            ax = axes[row, col]
            viz_cell_single_channel(
                crop,
                channel="GFP",
                ax=ax,
                axis_off=True,
                percentile_low=PERCENTILE_LOW,
                percentile_high=PERCENTILE_HIGH,
            )
            # Add scale bar to first cell of each allele block (top-left)
            if cell_idx == 0:
                add_scale_bar(ax)

    # Add titles above each block
    for allele_idx, allele in enumerate(allele_order):
        col_offset = allele_idx * N_COLS
        middle_col = col_offset + N_COLS // 2

        if allele_subtitles and allele in allele_subtitles:
            title = f"{allele}\n{allele_subtitles[allele]}"
        else:
            title = allele
        axes[0, middle_col].set_title(title, fontsize=14, fontweight='bold', pad=10)

    plt.tight_layout()

    # Add horizontal gap between left and right blocks
    plt.subplots_adjust(wspace=0.03, top=0.88 if suptitle else 0.95)

    # Manually shift right block to create gap
    for row in range(N_ROWS):
        for col in range(N_COLS, N_COLS * 2):
            pos = axes[row, col].get_position()
            axes[row, col].set_position([pos.x0 + 0.02, pos.y0, pos.width, pos.height])

    # Add suptitle if provided (centered above allele titles)
    if suptitle:
        fig.suptitle(suptitle, fontsize=12, fontweight='bold', y=0.97, x=0.5, ha='center')

    # Save figure
    fig.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"\nSaved: {output_path}")


# ============================================================================
# Main
# ============================================================================

def load_and_filter_data():
    """Load cell profiles and locations, apply quality filters."""
    print("\nLoading cell profiles...")
    profiles_b18 = pd.read_parquet(
        DATA_DIR / "raw" / "single_cell_profiles" / "batch18_ccm2_profiles_full_gfp.parquet"
    )
    profiles_b19 = pd.read_parquet(
        DATA_DIR / "raw" / "single_cell_profiles" / "batch19_ccm2_profiles_full_gfp.parquet"
    )
    profiles = pd.concat([profiles_b18, profiles_b19], ignore_index=True)
    print(f"  Loaded {len(profiles)} cell profiles")

    # Load cell locations (for pixel coordinates)
    print("Loading cell locations...")
    locations = pd.read_parquet(
        DATA_DIR / "interim" / "cell_locations" / "ccm2_cell_locations.parquet"
    )
    print(f"  Loaded {len(locations)} cell locations")

    # Create CellID in profiles to match locations
    profiles["CellID"] = (
        profiles["Metadata_Plate"] + "_" +
        profiles["Metadata_well_position"] + "_" +
        profiles["Metadata_ImageNumber"].astype(str) + "_" +
        profiles["Metadata_ObjectNumber"].astype(str)
    )

    # Merge to get pixel coordinates
    merged = profiles.merge(
        locations[["CellID", "Nuclei_Location_Center_X", "Nuclei_Location_Center_Y", "site"]],
        on="CellID",
        how="inner"
    )
    print(f"  Merged {len(merged)} cells with locations")

    # Add distance to edge
    merged["dist2edge"] = merged.apply(
        lambda row: compute_distance_to_edge(
            row["Nuclei_Location_Center_X"],
            row["Nuclei_Location_Center_Y"]
        ),
        axis=1
    )

    # Filter cells by quality
    print("\nFiltering cells by quality...")
    var_profiles = {}

    for allele in ALLELE_LIST:
        allele_df = merged[merged["Metadata_gene_allele"] == allele].copy()
        print(f"  {allele}: {len(allele_df)} total cells")

        # Apply edge distance filter
        filtered = allele_df[allele_df["dist2edge"] >= MIN_EDGE_DIST]
        print(f"    After edge filter (>={MIN_EDGE_DIST}): {len(filtered)}")

        # Intensity filter using percentiles (since data is normalized)
        # This approximates the notebook's filtering of "moderate" intensity cells
        intensity_vals = filtered[INTENSITY_FEATURE]
        p_low = intensity_vals.quantile(INTENSITY_PERCENTILE_LOW / 100)
        p_high = intensity_vals.quantile(INTENSITY_PERCENTILE_HIGH / 100)

        filtered = filtered[
            (filtered[INTENSITY_FEATURE] >= p_low) &
            (filtered[INTENSITY_FEATURE] <= p_high)
        ]
        print(f"    After intensity filter ({INTENSITY_PERCENTILE_LOW}-{INTENSITY_PERCENTILE_HIGH}%): {len(filtered)}")

        var_profiles[allele] = filtered.reset_index(drop=True)
        print(f"  {allele}: {len(filtered)} cells after filtering")

    return var_profiles


def generate_reproduce_figure(var_profiles: dict):
    """Generate figure using exact cells from reference notebook (10 cells per allele, 2x5 layout)."""
    print("\n" + "=" * 60)
    print("Generating REPRODUCE cell crop figure")
    print("Exact cells from: 6_visualize_cell_crop.ipynb")
    print("=" * 60)

    selected_indices = select_notebook_cells(var_profiles)

    # Check if we have enough cells
    for allele in ALLELE_LIST:
        n_found = len(selected_indices[allele])
        n_expected = len(NOTEBOOK_CELL_IDS[allele])
        if n_found < n_expected:
            print(f"  Warning: Only {n_found}/{n_expected} notebook cells found for {allele}")

    print("\nLoading cell crops...")
    all_crops = load_crops(var_profiles, selected_indices)

    output_path = OUTPUT_DIR / "FigB_reproduce_cell_crops.png"
    create_figure_small(all_crops, output_path)


def generate_feature_figure(var_profiles: dict):
    """Generate feature-based selection cell crop figure."""
    print("\n" + "=" * 60)
    print("Generating FEATURE-BASED cell crop figure")
    print("=" * 60)

    selected_indices, selection_info = select_feature_based_cells(
        var_profiles, SELECTION_FEATURE, N_CELLS
    )

    print("\nLoading cell crops...")
    all_crops = load_crops(var_profiles, selected_indices)

    # Create subtitles with percentile ranges
    allele_subtitles = {}
    for allele in ALLELE_LIST:
        info = selection_info[allele]
        allele_subtitles[allele] = f"({info['pct_range']} percentile)"

    # Create suptitle with feature name
    suptitle = f"Selection by: {SELECTION_FEATURE}"

    output_path = OUTPUT_DIR / "FigB_feat_sel_cell_crops.png"
    create_figure(all_crops, output_path, suptitle=suptitle, allele_subtitles=allele_subtitles)


def main():
    parser = argparse.ArgumentParser(description="Generate cell crop figures")
    parser.add_argument(
        "--mode",
        choices=["reproduce", "feature", "both"],
        default="reproduce",
        help="Selection mode: reproduce (default, exact notebook cells), feature (50 cells), or both"
    )
    args = parser.parse_args()

    print("=" * 60)
    print("CCM2 Cell Crop Figure Generator")
    print("=" * 60)

    # Load and filter data (shared between modes)
    var_profiles = load_and_filter_data()

    # Generate requested figures
    if args.mode in ["reproduce", "both"]:
        generate_reproduce_figure(var_profiles)

    if args.mode in ["feature", "both"]:
        generate_feature_figure(var_profiles)

    print("\n" + "=" * 60)
    print("Done!")
    print("=" * 60)


if __name__ == "__main__":
    main()
