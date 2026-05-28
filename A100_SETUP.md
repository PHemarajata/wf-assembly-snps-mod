# Continue this work on the A100 (Google Drive transfer)

Step-by-step to move the downsampled *B. pseudomallei* recombination-aware workflow —
and the Claude Code session — from the original workstation to the A100
(`/home/cdcadmin/wf-assembly-snps-mod`). All data is staged in Google Drive at
`wfsnps-a100-transfer/` (uploaded with rclone from the workstation).

## Prerequisites on the A100 (as `cdcadmin`)
- **Claude Code** installed and logged in (your own account).
- **VS Code** + the **Claude Code extension** (if you Remote-SSH in, install it in the *SSH: A100* window).
- **rclone**, **git**, **docker**, and a **conda env with Nextflow + Java 17+** (equivalent of `nextflow-wf`).

## 1. Point rclone at the same Drive
```bash
# Fresh auth (headless A100: run `rclone authorize "drive"` on a browser machine, paste the token):
rclone config        # create a remote named e.g.  gdrive  (type drive, the peerah@gmail.com account)
# — or copy the workstation's ~/.config/rclone/rclone.conf to the same path on the A100 to reuse its token.
```
Below uses `gdrive:` — substitute your remote's name.

## 2. Restore the Claude Code session
```bash
SRC="gdrive:wfsnps-a100-transfer"
rclone copy "$SRC/wfsnps-claude-session.tar.gz" ~/ -P
mkdir -p ~/.claude
tar xzf ~/wfsnps-claude-session.tar.gz -C ~/.claude/
# -> ~/.claude/projects/-home-cdcadmin-wf-assembly-snps-mod/  (transcript + sub-agents + memory/)
cat ~/.claude/RESTORE.md
```

## 3. Get the code
```bash
cd /home/cdcadmin
git clone https://github.com/PHemarajata/wf-assembly-snps-mod.git
cd wf-assembly-snps-mod
```

## 4. Pull the data + reference into the repo, fix paths
```bash
SRC="gdrive:wfsnps-a100-transfer"
rclone copy "$SRC/downsample_out"         downsample_out/         -P   #  32 MB samplesheet + Mash matrix
rclone copy "$SRC/cdc_fasta"              cdc_fasta/              -P   # 2.2 GB study isolates
rclone copy "$SRC/contextual_fasta"       contextual_fasta/       -P   # 4.8 GB (the 696 used genomes)
rclone copy "$SRC/reference"              reference/              -P   #  7 MB  K96243 reference
rclone copy "$SRC/results_downsampled_v3" results_downsampled_v3/ -P   # 5.4 GB prior results (optional)

# the samplesheet has absolute /home/phemarajata paths -> repoint to the A100:
sed -i 's#/home/phemarajata/#/home/cdcadmin/#g' downsample_out/samplesheet.csv

# sanity: every referenced input now exists (no output = good)
awk -F, 'NR>1{print $2}' downsample_out/samplesheet.csv | while read f; do [ -f "$f" ] || echo "MISSING $f"; done | head
```

## 5. Open in VS Code and resume the conversation
1. VS Code → **Open Folder** → `/home/cdcadmin/wf-assembly-snps-mod`.
2. Claude Code panel → **session history** → open the `8a38d8ef…` session (full transcript + memory).
   - CLI: `cd /home/cdcadmin/wf-assembly-snps-mod && claude --resume`.

## 6. (Optional) re-run on the A100's hardware
The reference is now at `reference/k96243.fasta` (repo-relative — nothing to configure):
```bash
conda activate <your-nextflow-env>     # Nextflow + Java 17+
nextflow run main.nf \
  -profile bp,dgx_station_a100_updated,docker \
  --input downsample_out/samplesheet.csv \
  --outdir results_a100 \
  --ref reference/k96243.fasta \
  --use_global_reference true \
  --recombination_aware_mode --integrate_results \
  --mash_sketch_size 100000 --mash_threshold 0.002 --max_cluster_size 100 \
  --merge_singletons false --run_gubbins --snp_package parsnp
```
The `dgx_station_a100_updated` profile is tuned for that box (more cores/RAM → more Gubbins in
parallel). The GPU does not accelerate Gubbins/IQ-TREE. The medoid bug is fixed on `main`, so
`--use_global_reference` is optional now (kept for cross-cluster SNP comparability).

## Notes
- Claude login is per-machine (the bundle carries history/memory, not credentials).
- Old `/home/phemarajata/...` paths in the transcript are cosmetic; `CLAUDE.md` (in git) + the
  memory files restore the substantive context. See memory `bp-recombination-aware-run-config.md`.
- `contextual_fasta/`, `cdc_fasta/`, `downsample_out/`, `reference/`, and `results_downsampled*/`
  are gitignored — they come via Drive, not `git clone`.
