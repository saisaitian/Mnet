# Mnet：an unsupervised multi-network driven drug discovery computational framework 💊

## These files are the key script for the GSFM Code used for the analysis of the paper:🌟
**A network medicine framework identifies Dihydroisotanshinone I as a SKP2-targeting anti-tumor agent for triple-negative breast cancer**

Saisai Tian, Jinyuan Lu, Chengyang Guo, Wenjing Gu, Pengli Huang, Xike Xu, Qun Wang, Weidong Zhang
## System Requirements 🛠

### Software Requirements
**OS:** Linux (Ubuntu 22.04 LTS or Rocky Linux 8.6+ recommended)  
**Environment Manager:** Miniconda/Mamba  

# STEP 1 Data processing
Description: Process data from TCGA and LINCS databases,and further measured the correlation between each sample and cell line, and removed samples that are not correlated to the cell lines.
## input 
     Cancer Gene Expression Matrix
     Drug-induced Gene Expression Matrix
     CCLE data
     
# STEP 2 Data transformation
Description: Calculation of the four biologically interpretable quantifiers (BIQs): GSFM_Up, GSFM_Down, GSFM_ssGSEA, and GSFM_TF.
## input
     Processing Cancer Gene Expression Matrix (step1)
     Processing Drug-induced Gene Expression Matrix (step1)
     Hallmark genesets
     Database:TF-target pairs by Garcia-Alonso et al.
## Build docker
     cd 2.data transformation/docker
     docker build --tag gsfm-script:latest --file GSFM.dockerfile .
     cd ..
     docker run -it -d  --restart=always --name GSFM-notebook   -p 12101:8888 --log-opt max-size=10m --log-opt max-file=5 -v `pwd`/project:/project   gsfm-script
     docker container exec -it GSFM-notebook bash
     jupyter-notebook --ip=0.0.0.0 --port=8888 --allow-root --no-browser
     run ./project/Notebooks BRCA_GSFM.ipynb

# STEP 3 Drug discovery
Description: Calculation of RS-GSFM, RS-Gene, RS-Viper                                                                                 
 step3.1:Calculation of RS-GSFM                                
 step3.2:Calculation of RS-Gene                                      
 step3.3:Calculation of RS-Viper                              
## input
     Processing Cancer Gene Expression Matrix (step1)
     Processing Drug-induced Gene Expression Matrix (step1) 
     Processing Cancer GSFM Active Matrix (step2)
     Processing Drug-induced GSFM Active Matrix  (step2)
#### Manual Installation Steps (Recommended):
```
1. Install Python 3.10
conda create -n unicure python=3.10
conda activate unicure

2. Install PyTorch (select appropriate CUDA version)
⚠ Check latest at: https://pytorch.org/get-started/locally/
pip3 install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118

3. Install Accelerate & DeepSpeed (recommended for reproduction)
pip install accelerate
⚠ Follow configuration (about how to configure DeepSpeed): https://github.com/huggingface/accelerate

4. Install core dependencies
pip install numpy pandas scikit-learn fastparquet tqdm anndata scanpy lora-pytorch scipy

5. Install Uni-Mol (required for testing)
You can create a new conda environment. 
https://github.com/deepmodeling/Uni-Mol
```

## Datasets Requirements 📚

### Step 1: Download Core Folders
Download and **overwrite** these folders to your local UniCure directories:

1. **[data folder](https://drive.google.com/drive/folders/1VPXl8h8iuhr8IdrAmAWQ9ldEHkWRCGWq?usp=drive_link)**  
   - Contains: LINCS, SciPlex and PTC datasets
   - Local path: `your_project_path/UniCure/data/`

2. **[requirement folder](https://drive.google.com/drive/folders/1VPXl8h8iuhr8IdrAmAWQ9ldEHkWRCGWq?usp=drive_link)**  
   - Contains: configuration files
   - Local path: `your_project_path/UniCure/requirement/`

3. **[model weights](https://drive.google.com/drive/folders/1qZ4QwEXST_FcIZDTizMu7aH-OXpdc_CB?usp=drive_link)**  
   - Contains: Pre-trained model weights & Dataset splits (Training, Validation, and Test sets)
   - Local path: `your_project_path/UniCure/result/`
   
> ⚠️ **Overwrite Notice**: Replace existing directories completely when copying <br> ⚠️ **Unzip Notice**: Unzip Unicure_best_model.rar

### Step 2: Download UCE Pretraining Files
Download these essential files to `requirement/UCE_pretraining_files/`:

| File | Size | Required Path |
|------|------|---------------|
| **[33l_8ep_1024t_1280.torch](https://figshare.com/articles/dataset/Universal_Cell_Embedding_Model_Files/24320806?file=43423236)** | 4.2 GB | `requirement/UCE_pretraining_files/` |
| **[all_tokens.torch](https://figshare.com/articles/dataset/Universal_Cell_Embedding_Model_Files/24320806?file=43423236)** | 780 MB | `requirement/UCE_pretraining_files/` |
| **[species_chrom.csv](https://figshare.com/articles/dataset/Universal_Cell_Embedding_Model_Files/24320806?file=43423236)** | 12 KB | `requirement/UCE_pretraining_files/` |

### Verification Checklist
After downloading, confirm directory structure:
```
UniCure/
├── data/
│   ├── lincs2020/
│ 	├── sciplex/
│	└── PTC/
├── result/
│   └── 11/
│       ├── lincs2020/
│       ├── sciplex3/ 
│	    └── sciplex4/
└── requirement/
    └── UCE_pretraining_files/
	    ├── protein_embeddings/
        ├── 33l_8ep_1024t_1280.torch
        ├── all_tokens.torch
        └── species_chrom.csv
```

## Quick Test :zap:
```
python quick_pred.py
```

## Training Reproduction :fire:

Our training pipeline is designed to be progressive. Later stages depend on the pre-trained weights generated by the previous stages (e.g., LINCS Step 1 -> LINCS Step 2 -> Sciplex 3 -> Sciplex 4). 

You can choose to run the entire pipeline at once or execute each stage individually.

### Option 1: Run the Entire Pipeline
To sequentially train and test all stages (LINCS Step 1 & 2, Sciplex 3, and Sciplex 4) with the default seed (`11`), simply run:

```bash
python main.py --run_all
```

### Option 2: Run Step-by-Step (Recommended)
If you want to train a specific stage or debug, you can run the stages individually using the provided flags. 

**1. LINCS Stage 1:**
```bash
python main.py --run_lincs1
```

**2. LINCS Stage 2 (Train & Test):**
*(Note: Ensure Stage 1 is completed and Cell Embedding is generated by `generate_emb.py` before running this)*
```bash
python main.py --run_lincs2
```

**3. Sciplex 3 (Train & Test):**
```bash
python main.py --run_sciplex3
```

**4. Sciplex 4 (Train & Test):**
```bash
python main.py --run_sciplex4
```

### Advanced: Customizing the Random Seed
By default, the pipeline uses `seed=11`. You can easily change the random seed for any execution by adding the `--seed` argument. For example, to run the Sciplex 3 stage with seed `42`:

```bash
python main.py --run_sciplex3 --seed 42
```

## Cell Embedding Generation (After Stage 1 Training) :fire:
```
python generate_emb.py
```

## Fine-tuning Reproduction :fire:
```
python finetune.py
```

# System requirement
R studio 4.1.1 R version 4.1.1 (2021-08-10)

Docker version:24.0.6

Platform: x86_64-pc-linux-gnu (64-bit)

Running under: Ubuntu 20.04.6 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/liblapack.so.3;  LAPACK version 3.9.0

locale:

 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8   
 [6] LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:

[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:

 [1] ggrepel_0.9.1         RColorBrewer_1.1.3    org.Hs.eg.db_3.12.0     
 [4] shinyWidgets_0.5.7    plotly_4.9.3          signatureSearch_1.4.6      
 [7] leaflet_2.1.1         DT_0.17               limma_3.46.0   
[10] shinyBS_0.61          ROCR_1.0.11           stringr_1.5.0        
[13] purrr_1.0.1           readr_1.4.0           tibble_3.1.8         
[16] tidyverse_1.3.0       tidyr_1.2.0           dplyr_1.0.8          
[19] ggplot2_3.3.5         pROC_1.17.0.1         memoise_2.0.1        
[22] cowplot_1.1.1          

loaded via a namespace (and not attached):
 [1] methods_4.3.1     rio_0.5.29        utf8_1.2.3        cellranger_1.1.0  readxl_1.4.2      magrittr_2.0.3    glue_1.6.2       
 [8] tibble_3.2.1      foreign_0.8-85    pkgconfig_2.0.3   lifecycle_1.0.3   utils_4.3.1       cli_3.6.1         zip_2.3.0        
[15] fansi_1.0.4       openxlsx_4.2.5.2  vctrs_0.6.3       graphics_4.3.1    data.table_1.14.8 grDevices_4.3.1   stats_4.3.1      
[22] compiler_4.3.1    forcats_1.0.0     haven_2.5.2       base_4.3.1        rstudioapi_0.14   tools_4.3.1       curl_5.0.2       
[29] hms_1.1.3         pillar_1.9.0      Rcpp_1.0.11       rlang_1.1.1       stringi_1.7.12    datasets_4.3.1 



