#!/usr/bin/env bash
set -euo pipefail

# ===============================
# PATMO Excel → reaction_network.ntw
# Ubuntu 22 兼容版本（内嵌 Python）
# ===============================

# 输入输出路径
INFILE="./tests/testfolder/Reaction_Network.xlsx"
OUTFILE="./tests/testfolder/reaction_network.ntw"

echo "==============================="
echo " PATMO Reaction Network Builder"
echo "==============================="
echo "[*] Input file : $INFILE"
echo "[*] Output file: $OUTFILE"

# -------------------------------
# ------ Dependency check -------
# -------------------------------
need_cmd() { command -v "$1" >/dev/null 2>&1; }

echo "[*] Checking python3 ..."
if ! need_cmd python3; then
  echo "python3 not found. Please install it: sudo apt install -y python3"
  exit 1
fi

echo "[*] Checking pip ..."
if ! python3 -m pip --version >/dev/null 2>&1; then
  echo "[*] Installing python3-pip ..."
  if need_cmd sudo; then
    sudo apt update && sudo apt install -y python3-pip
  else
    echo "pip not available and no sudo permission. Please install python3-pip manually." >&2
    exit 1
  fi
fi

echo "[*] Checking Python dependencies (pandas, openpyxl) ..."
if ! python3 - <<'PY' >/dev/null 2>&1
try:
    import pandas, openpyxl  # noqa
except Exception:
    raise SystemExit(1)
PY
then
  echo "[*] Installing pandas and openpyxl ..."
  if ! python3 -m pip install --user --upgrade pandas openpyxl >/dev/null 2>&1; then
    echo "[*] Retrying install with --break-system-packages ..."
    python3 -m pip install --user --upgrade --break-system-packages pandas openpyxl
  fi
fi

# -------------------------------
# Reaction Network Python Converter
# -------------------------------
python3 - "$INFILE" "$OUTFILE" <<'PY'
# -*- coding: utf-8 -*-
import sys
from pathlib import Path
import pandas as pd

EXPECTED_HEADER_LEADING = "@format:idx"
# internal column names (unique)
INTERNAL_COLUMNS = ["IDX","R1","R2","R3","P1","P2","P3","RATE"]
# exported header (duplicate names as required)
EXPORT_HEADER = ["@format:idx","R","R","R","P","P","P","rate"]
# Add a header
HEADER_COMMENT = "#@var:T=Tgas"  # <-- comment line to add

def find_columns(df):
    df.columns = [c.strip() if isinstance(c, str) else c for c in df.columns]
    if EXPECTED_HEADER_LEADING not in df.columns:
        raise ValueError(f"Required column '{EXPECTED_HEADER_LEADING}' not found.")
    idx_col = EXPECTED_HEADER_LEADING
    # handle repeated names: R, R.1, R.2 etc.
    r_cols = [c for c in df.columns if isinstance(c,str) and (c=="R" or c.startswith("R."))]
    p_cols = [c for c in df.columns if isinstance(c,str) and (c=="P" or c.startswith("P."))]
    # rate column
    if "k" in df.columns:
        rate_col = "k"
    elif "rate" in df.columns:
        rate_col = "rate"
    else:
        raise ValueError("Rate column 'k' or 'rate' not found.")
    # pad to 3
    r_cols = (r_cols + [None, None, None])[:3]
    p_cols = (p_cols + [None, None, None])[:3]
    return idx_col, r_cols, p_cols, rate_col

def s(series):
    """convert to string and keep blanks"""
    series = series.astype("string").fillna("")
    return series.apply(lambda x: x.strip() if isinstance(x,str) else x)

def build_output_df(df, idx_col, r_cols, p_cols, rate_col):
    out = pd.DataFrame()
    out["IDX"]  = s(df[idx_col])
    out["R1"]   = s(df[r_cols[0]]) if r_cols[0] else ""
    out["R2"]   = s(df[r_cols[1]]) if r_cols[1] else ""
    out["R3"]   = s(df[r_cols[2]]) if r_cols[2] else ""
    out["P1"]   = s(df[p_cols[0]]) if p_cols[0] else ""
    out["P2"]   = s(df[p_cols[1]]) if p_cols[1] else ""
    out["P3"]   = s(df[p_cols[2]]) if p_cols[2] else ""
    out["RATE"] = s(df[rate_col])
    return out[INTERNAL_COLUMNS]

def main():
    if len(sys.argv) < 3:
        raise SystemExit("Usage: python3 - <in.xlsx> <out.ntw>")
    in_path = Path(sys.argv[1]).expanduser().resolve()
    out_path = Path(sys.argv[2]).expanduser().resolve()
    if not in_path.exists():
        raise SystemExit(f"Input file not found: {in_path}")
    
    df = pd.read_excel(in_path, dtype="string", engine="openpyxl")
    if df.empty:
        raise SystemExit("Excel file is empty.")
    
    idx_col, r_cols, p_cols, rate_col = find_columns(df)
    out_df = build_output_df(df, idx_col, r_cols, p_cols, rate_col)
    
    tmp_path = out_path.with_suffix(".tmp")

    # write to temporary file first (without the header line)
    try:
        out_df.to_csv(tmp_path, index=False, header=EXPORT_HEADER, encoding="utf-8", line_terminator="\n")
    except TypeError:
        out_df.to_csv(tmp_path, index=False, header=EXPORT_HEADER, encoding="utf-8", lineterminator="\n")

    # prepend the comment line
    with open(out_path, "w", encoding="utf-8") as f_out, open(tmp_path, "r", encoding="utf-8") as f_in:
        f_out.write(f"{HEADER_COMMENT}\n")
        f_out.writelines(f_in.readlines())

    tmp_path.unlink(missing_ok=True)
    print(f"File generated successfully: {out_path}")


if __name__ == "__main__":
    main()
PY

echo "[✓] Conversion completed successfully."
echo "-----------------------------------"
echo " PATMO Reaction Network Compile Done"
echo "-----------------------------------"
