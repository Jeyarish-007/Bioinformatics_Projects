# Documentation: Flow Cytometry Data Analysis Using Non-Negative Matrix Factorization (NMF)

---

## 1. Title  
**Unsupervised Clustering of High-Dimensional Flow Cytometry Data Using Non-Negative Matrix Factorization (NMF)**  

---

## 2. Aim & Objectives  
### Aim  
To identify biologically meaningful cell subpopulations in high-dimensional flow cytometry data using **Non-Negative Matrix Factorization (NMF)**.  

### Objectives  
- **Preprocess** raw flow cytometry data (FCS files) for analysis.  
- **Reduce dimensionality** using PCA to determine optimal NMF components.  
- **Decompose** data into interpretable patterns using NMF.  
- **Visualize** cell clusters and marker contributions.  
- **Interpret** results to identify key cell populations.  

---

## 3. Introduction  
Flow cytometry generates high-dimensional data, making manual gating and analysis challenging. **Non-Negative Matrix Factorization (NMF)** is an unsupervised learning method that:  
- Decomposes data into **metagenes (W)** and **metacells (H)**.  
- Preserves non-negativity (biologically interpretable).  
- Helps identify **dominant cell populations** and **marker expression patterns**.  

This documentation details a Python-based workflow for:  
- Loading & normalizing FCS files.  
- Applying PCA to determine the optimal number of NMF components.  
- Clustering cells and visualizing results.  

---

## 4. Materials  
### Software & Libraries  
- **Python 3.10+**  
- **Libraries:**  
  - `FlowCal` (FCS file parsing)  
  - `scikit-learn` (PCA, NMF, t-SNE)  
  - `pandas`, `numpy` (data manipulation)  
  - `matplotlib`, `seaborn` (visualization)  

### Hardware  
- Minimum: **8GB RAM, 4-core CPU** (for 500k+ cells)  
- Recommended: **16GB+ RAM** for faster computation.  

### Dataset  
- **File:** `LUNG_PROJECT_11F24_24218_PRADYOT_GP1_PB_T1_ICP1.fcs`  
- **Dimensions:** **567,718 cells × 20 markers** (FSC-A, SSC-A, FITC-A, etc.)  

---

## 5. Methodology  
### Step 1: Data Loading & Preprocessing  
- **Load FCS file** using `FlowCal.io.FCSData()`.  
- **Convert to DataFrame** and **normalize** (min-max scaling).  

### Step 2: Dimensionality Reduction (PCA)  
- Apply PCA to determine the **optimal number of NMF components** (retaining 95% variance).  

### Step 3: Non-Negative Matrix Factorization (NMF)  
- Decompose data into:  
  - **W (Basis Matrix):** Cell-cluster associations (`n_cells × n_components`).  
  - **H (Coefficient Matrix):** Marker contributions (`n_components × n_markers`).  

### Step 4: Visualization  
1. **Heatmap** of marker contributions (`H` matrix).  
2. **t-SNE** projection of cell clusters.  

### Step 5: Interpretation  
- Identify **dominant markers per component**.  
- Relate findings to biological cell populations.  

---

## 6. Results  
### 1. Optimal NMF Components  
- PCA determined **3 components** explain **95% variance**.  

### 2. NMF Decomposition  
| Matrix | Shape | Description |  
|--------|-------|-------------|  
| **W** | `(567718, 3)` | Cell weights per component |  
| **H** | `(3, 20)` | Marker contributions per component |  

### 3. Marker Contributions (Heatmap)  
| Component | Dominant Markers | Interpretation |  
|-----------|------------------|-----------------|  
| **1** | SSC-A, BV786-A, BUV395-A | **Granularity & marker expression** |  
| **2** | FSC-A, FSC-H | **Cell size** |  
| **3** | Time | **Temporal effects** |  

### 4. t-SNE Clustering  
- Cells clustered into **3 distinct groups** (aligned with NMF components).  

---

## 7. Conclusion  
- NMF successfully **identified 3 major cell populations** in flow cytometry data.  
- **Component 1** (SSC-A, BV786-A) suggests **granular cells** (e.g., monocytes).  
- **Component 2** (FSC-A, FSC-H) represents **cell size variations**.  
- **Component 3** (Time) may indicate **instrument drift** or **staining artifacts**.  

### Future Work  
- **Validate clusters** with manual gating.  
- **Compare NMF with other methods** (e.g., UMAP, PhenoGraph).  
- **Apply to multi-sample datasets** for cohort analysis.  

---
