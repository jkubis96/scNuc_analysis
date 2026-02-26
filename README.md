### Organoid cell & nucleus analysis using JIMG_ncd

![Python version](https://img.shields.io/badge/python-%E2%89%A53.12%20%7C%20%3C3.13-blue?logo=python&logoColor=white.png)

</br>

<p align="right">
    <img src="https://github.com/jkubis96/Logos/blob/main/logos/jbs_current.png?raw=true" alt="drawing" width="250" />
    <img src="https://github.com/jkubis96/Logos/blob/main/logos/jbi_current.png?raw=true" alt="drawing" width="250" />
</p>

</br>

### Author: Jakub Kubiś 

<div align="left">
 Institute of Bioorganic Chemistry<br />
 Polish Academy of Sciences<br />
</div>


## Description

#### Flow Cytometry Analysis of Fused Dorsal–Ventral Organoids

This webpage presents the results of flow cytometry analysis performed on fused dorsal–ventral organoids using the **Amnis ImageStream** platform. The data were analyzed with the **JIMG_ncd Python library**.

The analysis included:

- Automated **nucleus detection**
- Assessment of **chromatin condensation status**
- Extraction and selection of **image-derived features** for further evaluation

We examined the relationships between **nuclear morphology**, **chromatin organization**, and **whole-cell structure** within organoids modeling **Huntington's disease (HD)**.

---

## CAG Repeat Variants (*HTT* / PolyQ Expansions)

The organoids represented different **CAG repeat lengths** in the *huntingtin (HTT)* gene, corresponding to **PolyQ expansions**:

- **21Q** – Control (healthy)
- **71Q** and **77Q** – Adult-onset HD
- **109Q** – Juvenile-onset HD

Feature selection enabled the identification of **distinct clusters of cells** sharing similar phenotypic profiles. This approach revealed:

- Significant differences in **cell abundance**
- Morphological differences between **HD cells** and cells derived from **healthy organoids**

</br>



<br />

## 📚 Table of Contents
- 1.[ Results](#res)
- 2.[ Prepare envireoment](#env)
- 3.[ Analysis scripts](#analysis)
 
<br />

# 1. Results <a id="res"></a>

Results analysis and description available here 👉 [Results 📄](https://jkubis96.github.io/scNuc_analysis/)





# 2. Prepare envireoment <a id="env"></a>



```bash
git clone https://github.com/jkubis96/scNuc_analysis.git

cd scNuc_analysis
```


```bash
poetry install --with dev --no-root
```


During the analysis, the following libraries are used:

* [PickShot](https://github.com/jkubis96/PickShot)  
* [JIMG_ncd](https://github.com/jkubis96/JIMG_ncd) 



# 3. Analysis scripts <a id="analysis"></a>

## Analysis Workflow

The analysis consisted of three main stages:

### 1. Initial Quality Assessment and Data Selection

The first stage involved a preliminary analysis of the results and manual selection of low-quality and high-quality images.  
These images were used to train the PickShot model to enable accurate image quality control (QC) in the main analysis pipeline.

Script available here 👉 [Results 📄](/scripts/prepare_pickshot_model.py)

### 2. Image Filtering and Feature Extraction

The second stage included the core analysis workflow.  
Images were filtered using the trained PickShot model to ensure high-quality input data.  

Subsequently, image-based features were extracted from imaging flow cytometry data acquired using the Amnis ImageStream system.  
Feature extraction was performed using the JIMG_ncd library.

Script available here 👉 [Results 📄](/scripts/data_prepare_features_selection.py)

### 3. Downstream Analysis

The final stage included:

- Clustering analysis  
- Data harmonization  
- Selection of cluster-discriminative features  
- Evaluation of the contribution of individual cell populations within each cluster  

This step enabled the assessment of cluster composition across different Huntington’s disease (HD) model organisms.

Script available here 👉 [Results 📄](/scripts/analysis_features_and_images_selection.py)


<br>
<br>


### Have fun JBS