# 🐭 miceLLC
A rapid and cost-effective machine learning–aided framework for pre-clinical screening of potent compounds and drugs against the Lewis Lung Cancer (LLC) Cell Line.

---

## 🎯 Overview
miceLLC provides a fast, accessible, and computationally affordable pipeline for identifying potent anticancer candidates against the Lewis Lung Cancer (LLC) cell line — a cornerstone model in pre-clinical lung cancer research.

By leveraging a validated **Random Forest–based machine learning model**, miceLLC predicts the cytotoxic potency (pIC₅₀) of compounds directly from their molecular structure (SMILES), enabling researchers to prioritize candidates long before costly wet-lab experiments are performed.

---

## 🚀 Key Applications
- **Pre-clinical Screening:** Rapidly filter large compound libraries for LLC potency before in vitro or in vivo testing
- **Drug Repurposing:** Identify approved or investigational drugs with untapped activity against LLC cells
- **Lead Optimization:** Use predicted pIC₅₀ values to guide structural modifications
- **Cost Reduction:** Replace expensive early-stage assays with rapid computational pre-screening
- **Decision Support:** Prioritize the most promising candidates for downstream experimental validation

---

## 📊 Model Performance
The Random Forest regression model achieved the highest predictive performance among all evaluated algorithms, with strong generalization across all dataset splits:

| Dataset            | R² Score |
|--------------------|----------|
| Training Set       | 0.9591   |
| Validation Set     | 0.8222   |
| Test Set           | 0.8069   |
| Overall (Combined) | 0.9244   |

Random Forest outperformed all benchmarked algorithms on the test set, including Gradient Boosting (R² = 0.7989), Support Vector Machine (R² = 0.7928), and K-Nearest Neighbors (R² = 0.7462), confirming its superior explanatory power and generalization capability.

Additional performance metrics:

| Metric                                          | Value  |
|-------------------------------------------------|--------|
| MAE (Validation)                                | 0.4247 |
| MAE (Test)                                      | 0.4058 |
| RMSE (Test)                                     | Lowest among all models |
| Pearson Correlation (actual vs. predicted pIC₅₀, Test) | 0.9622 |

The Pearson correlation of **0.9622** on the test set confirms the real-world applicability of the model.

---

## 🔬 Why Random Forest?
Among all models evaluated — Gradient Boosting, Support Vector Machine, K-Nearest Neighbors, Decision Tree, AdaBoost, Multi-Layer Perceptron, and Deep Learning — **Random Forest consistently achieved**:

- ✅ Highest R² on the test set (0.8069)
- ✅ Lowest test RMSE and MAE across all models
- ✅ Strongest Pearson correlation (0.9622) between actual and predicted pIC₅₀
- ✅ Minimal overfitting — stable and consistent performance across training, validation, and test phases

Other models exhibited reduced test R² scores and higher error metrics, reflecting less robust fit, potential overfitting, or underfitting issues.

---

## ⚙️ Installation

### Prerequisites
- Python ≥ 3.8
- `conda` package manager

### Install

**Ubuntu/Linux Terminal / Windows Command Line**

```bash
# Clone the repository
git clone https://github.com/aidrabd/miceLLC.git
cd miceLLC

# Make prediction script executable
chmod +x predict.py
```

**Step 1: Install Miniconda (if not already installed)**

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

**Step 2: Initialize and activate conda**

```bash
conda init
```

Restart your terminal, then activate the base environment:

```bash
conda activate
```

**Step 3: Create and activate a Python environment**

```bash
conda create -n py312 python=3.12.9
conda activate py312
python --version
```

---

## 🏃‍♂️ Simple Start

### Command Line Usage

```bash
python predict.py --input new_compounds.csv --model Random_Forest_model.h5
```

| Argument  | Description |
|-----------|-------------|
| `--input` | Path to your input CSV file (required) |
| `--model` | Path to the Random Forest model file (default: `Random_Forest_model.h5`) |

---

## 🧾 Input Format

Prepare a `.csv` file with the following columns:

| Column  | Description |
|---------|-------------|
| `SMILES` | Simplified Molecular Input Line Entry System notation of the compound |
| `pIC50`  | Biological potency value — **leave empty for prediction** |

See `sample.csv` in the repository to prepare your own input file.

---

## 📖 Scientific Background

The **Lewis Lung Carcinoma (LLC)** cell line is one of the most established syngeneic mouse models in oncology, widely used to evaluate the pre-clinical efficacy of novel anticancer agents in terms of tumor growth inhibition, immune response, and drug potency. Computationally predicting compound potency against LLC cells enables faster, cheaper, and more focused prioritization of candidates before laboratory testing.

miceLLC was designed to make this pre-clinical screening step rapid and cost-effective by integrating:

- **Machine Learning (Random Forest)** for accurate and robust pIC₅₀ prediction
- **Molecular Descriptor–Based Modeling** using SMILES as structural input
- **High-Throughput Virtual Screening** for rapid prioritization of large compound sets

---

## 💡 Highlights

- ⚡ **Rapid & cost-effective** — computationally screen compounds in minutes, before any lab work
- 🏆 **Top-performing model** — Random Forest outperforms all benchmarked algorithms
- 📈 **High predictive accuracy** — test R² of 0.8069 and Pearson correlation of 0.9622
- ⚖️ **Robust generalization** — consistent performance across training, validation, and test phases
- 🔬 **LLC-specific** — purpose-built for Lewis Lung Cancer Cell Line pre-clinical screening
- 🖥️ **Ready-to-use** — single command-line call with your own compound list

---

## ⚠️ Disclaimer

miceLLC is developed for **research purposes only**.
All predictions should be validated experimentally before any clinical or commercial application.

---

## 📄 License

This project is licensed under the **MIT License** — see the [LICENSE](LICENSE) file for details.

---

## 📚 Citation

If you use miceLLC in your research, please cite:

> *(Citation to be added upon publication)*
