# Quick Start: Integrated Scalable Workflow

## 🚀 Ready to Use!

The enhanced scalable workflow is now **complete** and ready for production use. It provides:

✅ **Core SNPs integrated between clusters**  
✅ **Trees combined through grafting algorithms**  
✅ **Global phylogenetic analysis**  
✅ **Comprehensive reporting**

## 🎯 Quick Start Commands

### 1. Basic Integrated Analysis
```bash
nextflow run main.nf \
  -profile docker \
  --input /path/to/assemblies \
  --outdir results_integrated \
  --scalable_mode \
  --integrate_results
```

### 2. Full Analysis with Gubbins
```bash
nextflow run main.nf \
  -profile docker \
  --input /path/to/assemblies \
  --outdir results_integrated \
  --scalable_mode \
  --integrate_results \
  --run_gubbins \
  --mash_threshold 0.05 \
  --merge_singletons
```

### 3. Large Dataset (1000+ genomes)
```bash
nextflow run main.nf \
  -profile docker \
  --input /path/to/assemblies \
  --outdir results_integrated \
  --scalable_mode \
  --integrate_results \
  --mash_threshold 0.1 \
  --max_cluster_size 50 \
  --merge_singletons \
  --mash_sketch_size 10000
```

## 📁 What You'll Get

### Integrated Results (`Integrated_Results/`)
- `integrated_core_snps.fa` - **Combined core SNPs from all clusters**
- `grafted_supertree.tre` - **Grafted tree combining all cluster trees**
- `integrated_core_snps.treefile` - **Global phylogenetic tree**
- `sample_cluster_mapping.tsv` - Sample-to-cluster assignments

### Final Reports (`Final_Results/`)
- `Scalable_Analysis_Final_Report.html` - **Interactive dashboard**
- `Analysis_Statistics.tsv` - Summary statistics

### Per-Cluster Results (`Clusters/`)
- Individual cluster alignments, trees, and Gubbins results

### Core SNPs (`Core_SNPs/`)
- Core SNPs extracted from each cluster

## ⚙️ Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--scalable_mode` | false | Enable clustering-based scalable analysis |
| `--integrate_results` | true | Enable core SNP integration and tree grafting |
| `--mash_threshold` | 0.03 | Distance threshold for clustering (lower = more clusters) |
| `--max_cluster_size` | 100 | Maximum samples per cluster |
| `--merge_singletons` | false | Merge singleton clusters for analysis |
| `--run_gubbins` | true | Enable recombination analysis |

## 🔍 Check Your Results

After completion, look for:

1. **Integration Success**: Check `Integrated_Results/` directory
2. **Final Report**: Open `Final_Results/Scalable_Analysis_Final_Report.html`
3. **Statistics**: Review `Final_Results/Analysis_Statistics.tsv`
4. **Cluster Details**: Explore individual `Clusters/cluster_*/` directories

## 📊 Expected Timeline

- **100 genomes**: 2-4 hours
- **500 genomes**: 6-12 hours
- **1000+ genomes**: 12-24 hours

## 🆘 Troubleshooting

### No integrated results?
- Ensure `--scalable_mode` and `--integrate_results` are enabled
- Check that clusters were formed (see `Summaries/clusters.tsv`)

### Empty alignments?
- Increase `--mash_threshold` (try 0.1)
- Use `--merge_singletons` to combine small clusters

### Need help?
- Check the detailed logs in `results/pipeline_info/`
- Review cluster formation in `Summaries/cluster_summary.txt`

---

## 🎉 You're All Set!

The integrated scalable workflow provides everything you requested:
- Core SNPs integrated between clusters
- Trees grafted together using algorithms
- Overall results with comprehensive reporting

**Run the workflow and get your integrated results!** 🚀