# School’s Out for the Summer: Cognition Varies Across the Calendar Year in Multiple Large-Scale Datasets
This contains the code used for all analyses and figures in the manuscript. 

## File organization 

```
/path/to/project/

    |-  Data/ Tabular data that may be shared on GitHub
    |-  Tables/        # csv outputs
    |-  Figs/        # Directory to collect figures output from R
    |- Tools/ # General functions needed  

```

## System Requirements
- **Operating System**: macOS 12+, Ubuntu 20.04+, Windows 10/11
- **R Version**: R 4.3.1 (2023-06-16) or later
- **Platform**: aarch64-apple-darwin20 (tested on macOS 12)
- **Dependencies**: The following R packages must be installed:
  - `tidyverse`, `lubridate`, `mgcv`, `gratia`, `gtsummary`, `data.table`, `purrr`, `effectsize`, `broom`, `circular`, `ggsci`, `factoextra`, `gt`
- **Tested On**:
  - R version 4.3.1 (Beagle Scouts)
  - macOS 12 (aarch64-apple-darwin20)


## Installation Guide
1. Install R from [CRAN](https://cran.r-project.org/).
2. Install RStudio (optional but recommended) from [RStudio](https://www.rstudio.com/).
3. Clone this repository:
   ```bash
   git clone https://github.com/yourusername/project-name.git
   cd project-name
   ```
4. Install required R packages
```r
install.packages(c("tidyverse", "lubridate", "mgcv", "gratia", "gtsummary", 
                   "data.table", "purrr", "effectsize", "broom", "circular", 
                   "ggsci", "factoextra", "gt"))
 ```

5. Place the repository folder in your working directory or adjust file paths as needed.

6. Typical install time: ~10 minutes on a normal desktop computer.

## Data Acquisition

The following datasets are required to run the analysis:

1. **ABCD Data (Adolescent Brain Cognitive Development Study)**
   - **Access Requirements**: The ABCD dataset can be obtained through the NIMH Data Archive (NDA).The process for accessing the data starts by creating an account at NDA (if you do not already have one), and then requesting access to the ABCD Study data through the dashboard there (https://nda.nih.gov/user/dashboard/data_permissions.html). More information can be found here (https://wiki.abcdstudy.org/faq/faq.html)
   - **Data Files**: Once you have access, download the relevant datasets for analysis. Place them in the `Data/ABCD` folder of this repository.

2. **ADHD1000 Data**
   - **Access Requirements**: The ADHD1000 dataset is hosted by the NDA
   - **Data Files**: After approval, download the data and place it in the `Data/ADHD1000` folder.

3. **GUSTO Data**
   - **Access Requirements**: The GUSTO dataset must be requested from here (https://gustodatavault.sg/about/request-for-data).
   - **Data Files**: After accessing the data, place it in the `Data/GUSTO` folder.

4. **PNC Data**
    - The PNC data needed to run the code are include in the repo in `Data/PNC`. 

---

### **Instructions for Preparing Data**

1. Once you have obtained the datasets, save them to the `Data/` directory of the project.
2. Ensure that the data is formatted correctly as specified in the project documentation or code (you may need to adjust the file paths in the code to match the location where you've saved the data).
3. Relevant data needed can be seen in the relevant code sections that load data.  

---

### **Instructions for running the code**
1. Make sure all other dependencies are installed (as mentioned in the Installation Guide).

2. Open and run the RMarkdown file (FirstSubmission_NHB.Rmd) in RStudio

3. Running all the code will take ~20 minutes on a laptop computer. 
