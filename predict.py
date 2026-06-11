import os
import sys
import argparse
import warnings
import numpy as np
import pandas as pd
import joblib

warnings.filterwarnings("ignore")

# ── RDKit import ──────────────────────────────────────────────────────────────
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
except ImportError:
    sys.exit(
        "ERROR: RDKit is not installed.\n"
        "Install with: pip install rdkit   or   conda install -c conda-forge rdkit"
    )

# ── Directory that holds all saved .h5 files ─────────────────────────────────
MODEL_DIR = ""

MORGAN_RADIUS = 2
MORGAN_NBITS  = 1024

def smiles_to_morgan(smiles_list):
    """
    Compute Morgan fingerprints for every SMILES in smiles_list.

    Returns
    -------
    feature_df    : pd.DataFrame  shape (n_valid, 1024)  columns Morgan_0..1023
    valid_indices : list of int   positions in smiles_list that succeeded
    failed_indices: list of int   positions that failed (invalid SMILES)
    """
    col_names = [f"Morgan_{i}" for i in range(MORGAN_NBITS)]
    rows          = []
    valid_indices = []
    failed_indices = []

    for idx, smi in enumerate(smiles_list):
        mol = Chem.MolFromSmiles(str(smi))
        if mol is None:
            print(f"  WARNING: Invalid SMILES at row {idx}: {smi!r}  —  skipped.")
            failed_indices.append(idx)
            continue

        fp   = AllChem.GetMorganFingerprintAsBitVect(mol, MORGAN_RADIUS, nBits=MORGAN_NBITS)
        bits = list(fp)          # list of 0/1 integers, length = MORGAN_NBITS
        rows.append(bits)
        valid_indices.append(idx)

    if not rows:
        return pd.DataFrame(), [], failed_indices

    feature_df = pd.DataFrame(rows, columns=col_names, index=valid_indices)
    return feature_df, valid_indices, failed_indices


# =============================================================================
# Step 2 — select only the training features, in training order
# =============================================================================

def select_training_features(feature_df, training_feature_names):
    """
    Keep only the Morgan bit columns that were selected during training and
    arrange them in exactly the same order the model was trained on.

    Any column missing from feature_df (should not happen with correct nBits)
    is filled with 0 and a warning is printed.
    """
    missing = [c for c in training_feature_names if c not in feature_df.columns]
    if missing:
        print(f"  WARNING: {len(missing)} expected feature column(s) not found after "
              f"fingerprint calculation — filled with 0.\n"
              f"  This usually means nBits or radius differs from training. "
              f"  Expected MORGAN_NBITS={MORGAN_NBITS}, MORGAN_RADIUS={MORGAN_RADIUS}.")
        for col in missing:
            feature_df[col] = 0

    return feature_df[training_feature_names]


# =============================================================================
# Main prediction routine
# =============================================================================

def run_prediction(input_csv, model_filename, predict_all=False):

    # ── 1. Load SMILES ────────────────────────────────────────────────────────
    if not os.path.isfile(input_csv):
        sys.exit(f"ERROR: Input file not found: {input_csv}")

    input_df = pd.read_csv(input_csv)

    if "SMILES" not in input_df.columns:
        sys.exit("ERROR: Input CSV must have a column named \'SMILES\'.")

    smiles_list = input_df["SMILES"].tolist()
    print(f"\nLoaded {len(smiles_list)} compound(s) from: {input_csv}")

    # ── 2. Compute Morgan fingerprints (identical to Step 1) ──────────────────
    print(f"\nCalculating Morgan fingerprints "
          f"(radius={MORGAN_RADIUS}, nBits={MORGAN_NBITS}) ...")

    feature_df, valid_indices, failed_indices = smiles_to_morgan(smiles_list)

    if feature_df.empty:
        sys.exit("ERROR: No valid SMILES — nothing to predict.")

    if failed_indices:
        print(f"  {len(failed_indices)} SMILES failed and were excluded.")

    # ── 3. Load training feature names (saved by the training script) ─────────
    feat_path = os.path.join(MODEL_DIR, "training_feature_names.h5")
    if not os.path.isfile(feat_path):
        sys.exit(
            f"ERROR: training_feature_names.h5 not found in \'{MODEL_DIR}\'\n"
            "  Run the Step 2 training script first so it saves this file."
        )

    training_feature_names = joblib.load(feat_path)
    print(f"  Training used {len(training_feature_names)} selected Morgan bit(s).")

    # ── 4. Align to training feature subset and order ─────────────────────────
    x_new = select_training_features(feature_df, training_feature_names)
    print(f"  Feature matrix ready: {x_new.shape}  "
          f"(compounds × selected Morgan bits)")

    # ── 5. Load model(s) ──────────────────────────────────────────────────────
    full_path = (model_filename if os.path.isabs(model_filename)
                 else os.path.join(MODEL_DIR, model_filename))

    if not os.path.isfile(full_path):
        sys.exit(f"ERROR: Model file not found: {full_path}")

    loaded = joblib.load(full_path)

    # ── 6. Predict ────────────────────────────────────────────────────────────
    output_df = pd.DataFrame({"SMILES": [smiles_list[i] for i in valid_indices]})

    if isinstance(loaded, dict):
        # models.h5  —  dict {model_name: fitted_model}
        print(f"\nRunning predictions with ALL {len(loaded)} model(s) ...")
        for name, mdl in loaded.items():
            col = "Predicted_pIC50_" + name.replace(" ", "_")
            output_df[col] = mdl.predict(x_new)
            print(f"  ✔  {name}")
    else:
        # Single .h5  (best_model.h5 or any individual model file)
        label = os.path.splitext(os.path.basename(model_filename))[0]
        print(f"\nRunning prediction with: {label}")
        output_df[f"Predicted_pIC50_{label}"] = loaded.predict(x_new)
        print(f"  ✔  {label}")

    # ── 7. Save output ────────────────────────────────────────────────────────
    out_path = "predictions_output.csv"
    output_df.to_csv(out_path, index=False)
    print(f"\n✅ Predictions saved to: {out_path}")
    print(output_df.to_string(index=False))


# =============================================================================
# Entry point
# =============================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Predict pIC50 for new compounds using a trained ML-QSAR model.\n"
            "Replicates the exact Step 1 Morgan fingerprint pipeline automatically."
        )
    )
    parser.add_argument(
        "--input", required=True,
        help="Path to input CSV file with a SMILES column."
    )
    parser.add_argument(
        "--model", default="best_model.h5",
        help=(
            "Model filename inside CS_Final/ (default: best_model.h5).\n"
            "Pass models.h5 with --all to predict with every trained model."
        )
    )
    parser.add_argument(
        "--all", action="store_true", dest="predict_all",
        help="When using models.h5, output a prediction column for every model."
    )

    args = parser.parse_args()
    run_prediction(args.input, args.model, args.predict_all)
