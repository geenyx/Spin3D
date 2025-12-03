# This repository contains the raw and processed data, as well as the analysis scripts, for the MnTe XMCD experiment described in the associated manuscript.

---

## Contents

### 1. Figure 3(b,c)

#### **Raw Data** 
file type: *.hdf5

**Energy scan (Positive / Negative helicity)**

| Energy (eV) | Positive | Negative |
|-------------|----------|----------|
| 640.55 | Sample_Image_2025-01-31_059 | Sample_Image_2025-01-31_060 |

#### **Analyzed Data**
file type: *.npy &*.png

- `cropped_norm_image_rcp_640.55eV_060`
- `cropped_aligned_image_xmcd_640.55eV_059_060`

#### **Python script**

- **`xmcd_single_overview_v2.ipynb`** – Generate overview XMCD image  

---

### 2. Figure 4(a,b)

#### **Raw Data** 
file type: *.hdf5

| Energy (eV) | Positive | Negative |
|-------------|----------|----------|
| 640.05 | Sample_Image_2025-01-31_025 | Sample_Image_2025-01-31_026 |
| 640.30 | Sample_Image_2025-01-31_027 | Sample_Image_2025-01-31_028 |
| 640.55 | Sample_Image_2025-01-31_032 | Sample_Image_2025-01-31_033 |
| 640.90 | Sample_Image_2025-01-31_035 | Sample_Image_2025-01-31_036 |
| 641.40 | Sample_Image_2025-01-31_038 | Sample_Image_2025-01-31_039 |
| 641.80 | Sample_Image_2025-01-31_043 | Sample_Image_2025-01-31_044 |
| 642.45 | Sample_Image_2025-01-31_046 | Sample_Image_2025-01-31_047 |
| 643.00 | Sample_Image_2025-01-31_050 | Sample_Image_2025-01-31_051 |
| 644.40 | Sample_Image_2025-01-31_053 | Sample_Image_2025-01-31_054 |
| 651.40 | Sample_Image_2025-01-31_066 | Sample_Image_2025-01-31_067 |
| 652.00 | Sample_Image_2025-01-31_069 | Sample_Image_2025-01-31_070 |
| 653.30 | Sample_Image_2025-01-31_071 | Sample_Image_2025-01-31_072 |

#### **Analyzed Data**
file type: *.npy &*.png

- `cropped_aligned_images_xmcd_640.05eV_025_026`
- `cropped_aligned_images_xmcd_640.30eV_027_028`
- `cropped_aligned_images_xmcd_640.55eV_032_033`
- `cropped_aligned_images_xmcd_640.90eV_035_036`
- `cropped_aligned_images_xmcd_641.40eV_038_039`
- `cropped_aligned_images_xmcd_641.80eV_043_044`
- `cropped_aligned_images_xmcd_642.45eV_046_047`
- `cropped_aligned_images_xmcd_643.00eV_050_051`
- `cropped_aligned_images_xmcd_644.40eV_053_054`
- `cropped_aligned_images_xmcd_651.40eV_066_067`
- `cropped_aligned_images_xmcd_652.00eV_069_070`
- `cropped_aligned_images_xmcd_653.30eV_071_072`

#### **Python script**

- **`xmcd_series_v2.ipynb`** – Generate XMCD images from energy series  

---

### 3. Figure 4(c)

#### **Raw Data** 
file type: *.hdf5

| Energy range (eV) | Positive | Negative |
|-------------------|----------|----------|
| 630–662 | Sample_Line_2025-01-30_093 | Sample_Line_2025-01-30_094 |

#### **Analyzed Data**
file type: *.txt

- `energy_image`
- `xmcd_int0_image`
- `xmcd_int1_image`
- `xmcd_intdiff_image`
- `xmcd_intdiff_std_image` 

- `energy_line`
- `dark_transmission_lcp_line`
- `dark_transmission_rcp_line`
- `dark_transmission_lcp_line_filtered`
- `dark_transmission_rcp_line_filtered`
- `xmcd_line.txt`
- `dark_xmcd_line.txt`

#### **Python script**

- **`transmissionSpectra_image_v1.ipynb`** – Calculate XMCD spectra from XMCD images  
- **`transmissionSpectra_v3.ipynb`** – Generate XMCD spectra from line scans  

---

### 4. Figure 5

#### **Raw Data** 
file type: *.hdf5

| Energy (eV) | Positive | Negative |
|-------------|----------|----------|
| 640.55 | Sample_Image_2025-01-31_061 | Sample_Image_2025-01-31_062 |

#### **Analyzed Data**
file type: *.npy &*.png

- `cropped_aligned_image_xmcd_640.55eV_061_062`

#### **Python script**
- **`xmcd_single_zoom_v2.ipynb`** – Generate zoomed XMCD images 

---

## Related Publication

If you use this dataset, please cite:  
> R. Yamamoto et al., arXiv:2502.18597 [cond-mat.mtrl-sci] (2025).
---