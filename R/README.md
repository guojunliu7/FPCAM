# Cell Annotation Tool Usage Guide

This tool is designed to annotate single-cell clustering results based on a marker gene matrix and a dictionary file. It requires the upload of two files to perform cell type annotation analysis and supports downloading the annotated results along with downstream analysis outputs.

---

## 📁 1. Upload Files

Click **"Cell Annotation"** on the sidebar to upload the following two files:

### 🔹 Marker Dictionary File (.xlsx)

- Contains information linking clusters with marker genes 
- A standard single-cell transcriptome analysis workflow was performed using the Seurat package. First, the raw 10X-format data were loaded and used to construct a Seurat object, during which low-quality cells and lowly expressed genes were filtered out. The data were then normalized, and highly variable genes were identified. Next, the data were scaled, and dimensionality reduction was performed using principal component analysis (PCA). Marker genes for each cell cluster were identified, and finally, the average expression levels of all differentially expressed genes within each cluster were calculated to support downstream analyses.
- Please rename the file to <span style="color:blue"><b>Dictionary.xlsx</b></span> before uploading.


### 🔹 Cluster Expression Data File

- Contains gene expression profiles for each cluster (genes as row names)  
- Please rename the file to <span style="color:blue"><b>Marker_expression.xlsx</b></span> before uploading.



---

## ▶️ 2. Click the "Run" Button

After clicking the **Run** button, the system will start the annotation process. This computation may take a while.


## 📊 3. View Annotation Results

Once the process is complete, please click the **Result** button to view and download.

---

## 🔁 4. Reset and Re-upload

Click the **Reset** button to clear all uploaded files and begin a new analysis.

---

## 📌 Notes

- File names should not contain spaces or special characters.
- One set of marker + expression data can be processed at a time.
- Example files are available in the `example_data/` directory.

---

If you need further assistance, please contact the developer.

- Pro. Guojun Liu:gjliu77@gmail.com
- Mr. Yan Shi:1490261148@qq.com

