#!/usr/bin/env bash
set -euo pipefail

# --- Prompt user for subfolder name under ./tests ---
while true; do
  read -r -p "[?] Enter the folder name you prepared under ./tests folder: " USER_SUBDIR
  if [[ -n "${USER_SUBDIR// }" ]]; then
    break
  fi
  echo "[-] Folder name cannot be empty. Please try again."
done

BASE_DIR="./tests/${USER_SUBDIR}"
if [[ ! -d "$BASE_DIR" ]]; then
  echo "Folder not found: $BASE_DIR"
  echo "    Please create it and rerun."
  exit 1
fi

INFILE="${BASE_DIR}/reaction_network.xlsx"
OUTFILE="${BASE_DIR}/reaction_network.ntw"

echo "[*] Using base folder : $BASE_DIR"
echo "[*] Input Excel       : $INFILE"
echo "[*] Output NTW        : $OUTFILE"

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
INTERNAL_COLUMNS = ["IDX","R1","R2","R3","P1","P2","P3","RATE"]
EXPORT_HEADER = ["@format:idx","R","R","R","P","P","P","rate"]
HEADER_COMMENT = "#@var:T=Tgas"

def find_columns(df):
    df.columns = [c.strip() if isinstance(c, str) else c for c in df.columns]
    if EXPECTED_HEADER_LEADING not in df.columns:
        raise ValueError(f"Required column '{EXPECTED_HEADER_LEADING}' not found.")
    idx_col = EXPECTED_HEADER_LEADING
    r_cols = [c for c in df.columns if isinstance(c,str) and (c=="R" or c.startswith("R."))]
    p_cols = [c for c in df.columns if isinstance(c,str) and (c=="P" or c.startswith("P."))]
    if "k" in df.columns:
        rate_col = "k"
    elif "rate" in df.columns:
        rate_col = "rate"
    else:
        raise ValueError("Rate column 'k' or 'rate' not found.")
    r_cols = (r_cols + [None, None, None])[:3]
    p_cols = (p_cols + [None, None, None])[:3]
    return idx_col, r_cols, p_cols, rate_col

def s(series):
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
    try:
        out_df.to_csv(tmp_path, index=False, header=EXPORT_HEADER, encoding="utf-8", line_terminator="\n")
    except TypeError:
        out_df.to_csv(tmp_path, index=False, header=EXPORT_HEADER, encoding="utf-8", lineterminator="\n")
    with open(out_path, "w", encoding="utf-8") as f_out, open(tmp_path, "r", encoding="utf-8") as f_in:
        f_out.write(f"{HEADER_COMMENT}\n")
        f_out.writelines(f_in.readlines())
    tmp_path.unlink(missing_ok=True)
    print(f"File generated successfully: {out_path}")

if __name__ == "__main__":
    main()
PY

echo "Reaction network conversion completed successfully."

# -------------------------------
# Write copylist.pcp file
# -------------------------------
COPYLIST_FILE="${BASE_DIR}/copylist.pcp"

echo "[*] Creating copylist.pcp ..."
cat > "$COPYLIST_FILE" <<EOF
test.f90
profile.dat
solar_flux.txt
EOF
echo "copylist.pcp created at: $COPYLIST_FILE"

# -------------------------------
# Write test.f90 file
# -------------------------------
#TEST_F90_FILE="${BASE_DIR}/test.f90"
#
#echo "[*] Creating test.f90 ..."
#cat > "$TEST_F90_FILE" <<'EOF'
#program test
#  use patmo
#  use patmo_commons
#  use patmo_constants
#  use patmo_parameters
#  implicit none
#  real*8::dt,x(speciesNumber),t,tend,imass,one_year
#  real*8::heff(chemSpeciesNumber)
#  real*8::dep(chemSpeciesNumber)
#  integer::icell,i,j
#  real*8::convergence = 100.0
#
#  call patmo_init()
#  call patmo_loadInitialProfile("profile.dat",unitH="km",unitX="1/cm3")
#  call patmo_setFluxBB()
#  call patmo_setGravity(9.8d2)
#  wetdep(:,:) = 0d0 
#
#  imass = patmo_getTotalMass()
#  print*,"mass:",imass
# 
#  dt = secondsPerDay 
#  tend = secondsPerDay*365*20d0 
#  one_year = secondsPerDay*365
#  t = 0d0
#
#  do
#     call patmo_run(dt,convergence)
#     t = t + dt        
#     if (t==tend .or. abs(convergence) < 1e-10) then
#        call patmo_run(dt,convergence)
#        call patmo_dumpDensityToFile(29,t,patmo_idx_O2)
#     endif  
#     print '(F11.2,a2)',t/tend*1d2," %"
#     if(t>=tend) exit
#  end do
#  
#  imass = patmo_getTotalMass()
#  print *,"mass:",imass
#  call patmo_dumpHydrostaticProfile("hydrostatEnd.out")
#  call patmo_dumpJValue("jvalue.dat")
#end program test
#EOF
#echo "test.f90 created at: $TEST_F90_FILE"

# -------------------------------
# Convert settings.xlsx → options.opt
# -------------------------------
SETTINGS_XLSX="${BASE_DIR}/settings.xlsx"
OPTIONS_OPT="${BASE_DIR}/options.opt"

echo "[*] Converting settings.xlsx to options.opt ..."
python3 - "$SETTINGS_XLSX" "$OPTIONS_OPT" "$BASE_DIR" <<'PY'
import sys
from pathlib import Path
import pandas as pd
import os

if len(sys.argv) < 4:
    raise SystemExit("Usage: python3 - <settings.xlsx> <options.opt> <BASE_DIR>")

in_path = Path(sys.argv[1]).expanduser().resolve()
out_path = Path(sys.argv[2]).expanduser().resolve()
base_dir = Path(sys.argv[3]).expanduser().resolve()

if not in_path.exists():
    raise SystemExit(f"settings.xlsx not found: {in_path}")

df = pd.read_excel(in_path, dtype="string", engine="openpyxl")
expected = ["Parameter", "unit", "input"]
if not all(col in df.columns for col in expected):
    raise SystemExit("settings.xlsx must have columns: Parameter, unit, input")
df = df.fillna("")
relative_path = os.path.relpath(base_dir / "reaction_network.ntw", Path.cwd()).replace("./", "")
with open(out_path, "w", encoding="utf-8") as f:
    f.write(f"network = {relative_path}\n")
    for _, row in df.iterrows():
        param = str(row["Parameter"]).strip()
        value = str(row["input"]).strip() if isinstance(row["input"], str) else ""
        f.write(f"{param} = {value}\n")
print(f"options.opt created at: {out_path}")
PY
echo "options.opt created successfully."

# -------------------------------
# Convert profile.xlsx → profile.dat
# -------------------------------
PROFILE_XLSX="${BASE_DIR}/profile.xlsx"
PROFILE_DAT="${BASE_DIR}/profile.dat"

echo "[*] Converting profile.xlsx to profile.dat ..."
python3 - "$PROFILE_XLSX" "$PROFILE_DAT" <<'PY'
import sys
from pathlib import Path
import pandas as pd
import numpy as np

if len(sys.argv) < 3:
    raise SystemExit("Usage: python3 - <profile.xlsx> <profile.dat>")

in_path = Path(sys.argv[1]).expanduser().resolve()
out_path = Path(sys.argv[2]).expanduser().resolve()

if not in_path.exists():
    raise SystemExit(f"profile.xlsx not found: {in_path}")

df = pd.read_excel(in_path, engine="openpyxl")
df = df.where(pd.notnull(df), "")

cols_keep = [0, 1]
for i, col in enumerate(df.columns):
    if i not in cols_keep:
        df[col] = df[col].apply(lambda x: f"{float(x):.4E}" if isinstance(x,(int,float,np.integer,np.floating)) else x)

total_columns = len(df.columns)
env_params = 5
species_count = total_columns - env_params

tmp_path = out_path.with_suffix(".tmp")
try:
    df.to_csv(tmp_path, sep="\t", index=False, header=True, encoding="utf-8", line_terminator="\n")
except TypeError:
    df.to_csv(tmp_path, sep="\t", index=False, header=True, encoding="utf-8", lineterminator="\n")

with open(out_path, "w", encoding="utf-8") as fout:
    fout.write(f"{env_params}\t{species_count}\n")
    with open(tmp_path, "r", encoding="utf-8") as fin:
        fout.writelines(fin.readlines())
tmp_path.unlink(missing_ok=True)
print(f"profile.dat created with header line: {env_params} {species_count}")
PY
echo "profile.dat created successfully."

# -------------------------------
# Build solar_flux.txt from solar_flux.xlsx + options.opt
# -------------------------------
SOLAR_XLSX="${BASE_DIR}/solar_flux.xlsx"
SOLAR_TXT="${BASE_DIR}/solar_flux.txt"
OPTIONS_OPT="${BASE_DIR}/options.opt"

echo "[*] Building solar_flux.txt by interpolating to wavelength grid (nm) ..."

python3 - "$SOLAR_XLSX" "$SOLAR_TXT" "$OPTIONS_OPT" <<'PY'
import sys
from pathlib import Path
import pandas as pd
import numpy as np

if len(sys.argv) < 4:
    raise SystemExit("Usage: python3 - <solar_flux.xlsx> <solar_flux.txt> <options.opt>")

in_xlsx = Path(sys.argv[1]).expanduser().resolve()
out_txt = Path(sys.argv[2]).expanduser().resolve()
opt_path = Path(sys.argv[3]).expanduser().resolve()

if not in_xlsx.exists():
    raise SystemExit(f"solar_flux.xlsx not found: {in_xlsx}")
if not opt_path.exists():
    raise SystemExit(f"options.opt not found: {opt_path}")

params = {}
with open(opt_path,"r",encoding="utf-8") as f:
    for line in f:
        line=line.strip()
        if not line or line.startswith("#"):
            continue
        if "=" in line:
            k,v=line.split("=",1)
            params[k.strip()] = v.strip()

def get_num(key,alt_keys=(),cast=float,required_name=None):
    keys_to_try = (key,)+tuple(alt_keys)
    for kk in keys_to_try:
        if kk in params and params[kk] != "":
            try:
                return cast(params[kk])
            except Exception:
                raise SystemExit(f"Invalid value for '{kk}': {params[kk]}")
    readable = required_name or key
    alts = ", ".join(alt_keys) if alt_keys else ""
    hint = f" (or {alts})" if alts else ""
    raise SystemExit(f"Missing '{readable}'{hint} in options.opt")

lam_min = get_num("wavelengMin",alt_keys=("wavelengMin",),cast=float,required_name="wavelengMin")
lam_max = get_num("wavelengMax",cast=float,required_name="wavelengMax")
NBIN = get_num("photoBinsNumber",cast=int)

if NBIN <= 1:
    raise SystemExit("photoBinsNumber must be > 1")
if lam_min > lam_max:
    lam_min, lam_max = lam_max, lam_min

df = pd.read_excel(in_xlsx,engine="openpyxl")
cols = {c.strip():c for c in df.columns if isinstance(c,str)}
wcol = next((cols[c] for c in cols if c.lower().startswith("wavelength")),None)
icol = next((cols[c] for c in cols if c.lower().startswith("irradiance")),None)
if wcol is None or icol is None:
    raise SystemExit("Expected two columns: 'Wavelength (nm)' and 'Irradiance photon/(cm2 s nm)'")

df = df[[wcol,icol]].rename(columns={wcol:"lambda_nm",icol:"Iph"})
df = df.dropna(subset=["lambda_nm","Iph"])
df["lambda_nm"] = pd.to_numeric(df["lambda_nm"],errors="coerce")
df["Iph"] = pd.to_numeric(df["Iph"],errors="coerce")
df = df.dropna(subset=["lambda_nm","Iph"])
df = df.sort_values("lambda_nm")

lambda_grid = np.linspace(lam_min,lam_max,NBIN)
lam_src = df["lambda_nm"].to_numpy()
I_src = df["Iph"].to_numpy()
I_out = np.interp(lambda_grid,lam_src,I_src,left=0.0,right=0.0)

np.savetxt(out_txt,I_out,fmt="%.8e")
print(f"solar_flux.txt created at: {out_txt}")
PY

echo "solar_flux.txt created successfully."

# -------------------------------
# Run PATMO test with user-defined folder, then cd ./build
# -------------------------------
echo "[*] Running PATMO test with: -test=\"${USER_SUBDIR}\" ..."
python3 patmo -test="${USER_SUBDIR}"

echo "[*] Entering ./build ..."
cd ./build
