# Bioinformatics Portfolio

A showcase of bioinformatics pipelines, ML/AI dashboards, and production ready tooling developed by Myren Sutton. This project demonstrates genomic data engineering, visualization, RAG assistants, and cloud based execution.

## Structure
- `pipelines/`: Nextflow variant calling workflows
- `api/`: Flask REST API for genomic data querying
- `dashboard/`: Streamlit dashboard to explore variant and clustering data
- `rag-assistant/`: LangChain powered assistant trained on genomics literature
- `cloud-runner/`: AWS scripts and notebooks for cloud-based processing

## Cleaning Up Nextflow Logs and Work Directories

Nextflow creates `work/` and `.nextflow/` directories to store intermediate files and pipeline metadata. These are needed for resuming or debugging runs, but can take up a lot of space.

- **Keep them** if you want to resume or debug failed runs.
- **Delete them** after a successful run if you only need the final results.

To clean up safely, run:

```bash
nextflow clean -f
rm -rf work/ .nextflow/
```

This will remove all intermediate files and logs, freeing up disk space.# Test workflow trigger
