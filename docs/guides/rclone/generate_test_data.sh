#!/usr/bin/env bash
# Generates a tree of dummy data (random, incompressible) that mimics
# raw sequencing data, to test an rclone transfer.
# Usage: ./generate_test_data.sh TARGET_DIR [VOLUME_IN_GB]
#   VOLUME_IN_GB < 10  -> "quick test" set (~1 GB, a few minutes)
#   VOLUME_IN_GB >= 10 -> ~VOLUME_IN_GB GB (70 % as 1 GiB files, 30 % as 10 GiB files)
set -euo pipefail

DEST="${1:?Usage: $0 TARGET_DIR [VOLUME_IN_GB]}"
VOL="${2:-1}"

if [ -e "$DEST" ] && [ -n "$(ls -A "$DEST" 2>/dev/null)" ]; then
  echo "Error: $DEST exists and is not empty." >&2
  exit 1
fi

# Check free space before starting (requested volume + 10 %)
mkdir -p "$DEST"
NEEDED_GB=$(( (VOL < 10 ? 2 : VOL) * 11 / 10 ))
FREE_GB=$(( $(df -Pk "$DEST" | awk 'NR==2 {print $4}') / 1024 / 1024 ))
if [ "$FREE_GB" -lt "$NEEDED_GB" ]; then
  echo "Error: $NEEDED_GB GB of free space needed, only $FREE_GB GB available on this disk." >&2
  rmdir "$DEST" 2>/dev/null
  exit 1
fi

if [ "$VOL" -lt 10 ]; then
  N_MEDIUM=3;  S_MEDIUM=100M   # 3 x 100 MiB: below rclone's multipart threshold (200 MiB)
  N_LARGE=1;   S_LARGE=600M    # 1 x 600 MiB: triggers a multipart upload
else
  N_MEDIUM=$(( VOL * 70 / 100 )); S_MEDIUM=1G
  N_LARGE=$(( VOL * 30 / 1000 )); S_LARGE=10G
fi

mkdir -p "$DEST/run_01/small" "$DEST/run_01/medium" "$DEST/run_02/large" \
         "$DEST/folder with spaces/données accentuées"

# 2000 small files of 100 KiB (many objects, little volume)
for i in $(seq -w 1 2000); do
  head -c 100K /dev/urandom > "$DEST/run_01/small/read_$i.bin"
done

for i in $(seq -w 1 "$N_MEDIUM"); do
  head -c "$S_MEDIUM" /dev/urandom > "$DEST/run_01/medium/sample_$i.fastq.gz"
done

for i in $(seq -w 1 "$N_LARGE"); do
  head -c "$S_LARGE" /dev/urandom > "$DEST/run_02/large/alignment_$i.bam"
done

# Tricky names: spaces, accents, parentheses
echo "witness file" > "$DEST/folder with spaces/données accentuées/échantillon (copie) n°1.txt"

echo "Done:"
du -sh "$DEST"
find "$DEST" -type f | wc -l | xargs echo "files:"
