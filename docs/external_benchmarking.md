# External benchmarking

The external benchmarking workflow maps spatial transcriptomic units to the GSE157329 CS12-CS14 human organogenesis single-cell reference.

## Main steps

1. Load the spatial transcriptomic query Seurat object.
2. Load the GSE157329 reference count matrix and annotations.
3. Keep CS12 and CS13-14 reference cells.
4. Normalize, select variable genes, run PCA and UMAP on the reference.
5. Use Seurat `FindTransferAnchors`, `TransferData` and `MapQuery` to transfer fine reference labels.
6. Collapse transferred fine labels into the 17 broad anatomical classes used in the manuscript using a curated correspondence table.
7. Export prediction scores, a row-normalized confusion matrix and spot-level predictions.

The script is provided at:

- `analysis/label_transfer_xu2023/run_label_transfer_GSE157329_CS12_CS14.R`
