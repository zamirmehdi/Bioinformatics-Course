# Virus Classification Using Neural Networks

[![Course](https://img.shields.io/badge/Course-Bioinformatics-blue.svg)](#)
[![University](https://img.shields.io/badge/University-AUT-red.svg)](https://aut.ac.ir/en)

<details><summary><h2>📚 Table of Contents</h2></summary>

- [Overview](#-overview)
- [Dataset Description](#-dataset-description)
- [Feature Extraction](#-feature-extraction)
- [Model Architecture](#-model-architecture)
- [Training and Validation](#-training-and-validation)
- [Test Prediction](#-test-prediction)
- [Results and Performance](#-results-and-performance)
- [Project Structure](#-project-structure)
- [Tools and Dependencies](#-tools-and-dependencies)
- [Contact](#-contact)

</details>

---

## 📋 Overview

This project focuses on classification of six virus types using neural network models applied on viral genome nucleotide sequences. Each viral genome is represented by sequences of variable length, classified into classes labeled Class1 through Class6.

The goal is to accurately predict virus class membership from genomic data, using innovative preprocessing and deep learning techniques to handle variable-length sequence data efficiently.

---

## 📂 Dataset Description

- **Training data:** 1320 labeled samples, evenly distributed over 6 classes (220 samples per class).
- **Development data:** 180 samples (30 per class) used strictly for model evaluation, not training.
- **Test data:** 400 unlabeled genome sequences requiring prediction outputs.

All datasets are provided in CSV format with nucleotide sequences.

---

## 🔍 Feature Extraction

- Applied a **k-mer frequency extraction** method with \( k=2 \).
- Sequences transformed into normalized 2-mer frequency vectors.
- This converts variable-length sequences into fixed-length numeric inputs, preserving local nucleotide patterns essential for classification.
- 2-mer length was experimentally optimized for best model accuracy and generalization.

---

## 🧠 Model Architecture

- Implemented a **Multi-Layer Perceptron (MLP)** classifier using `scikit-learn`.
- Architecture:
  - Three hidden layers, each with 64 neurons.
  - ReLU activation functions for non-linearity.
- Hyperparameters (random state, max iterations) chosen to ensure stable and convergent training.
- MLP balances model capacity with risk of overfitting, suitable for dataset size.

---

## 🏋️‍♂️ Training and Validation

- Model trained using the training dataset (1320 samples).
- Evaluated using development dataset exclusively for unbiased accuracy assessment.
- Achieved near-perfect accuracy on development set, indicating solid generalization.
- Care taken to avoid development data leakage during training.

---

## 📊 Test Prediction

- Model applied to test dataset sequences transformed via the same 2-mer method.
- Outputs predicted classes as integers 1 through 6.
- Predictions formatted as a list for submission compliance on Quera.

---

## 📈 Results and Performance

- Achieved approximately **100% accuracy on development set**.
- Precise handling of sequence variability and balanced model architecture contributed to performance.
- Model predictions on test set show strong consistency and reliability.

---

## 🗂️ Project Structure
```
Virus Classification (Final Project)/
├── data/
│   ├── development_set.csv
│   ├── test_set.csv
│   └── training_set.csv
├── docs/
│   ├── Instruction.pdf # Project instructions (Persian)
│   └── Report.pdf # Comprehensive report (Persian)
└── src/
└── BioInformatics_FinalProject.ipynb # Jupyter notebook with preprocessing, model code, and results analysis
```
---

## 🛠️ Tools and Dependencies

- Python 3.x
- pandas, numpy for data handling
- scikit-learn for MLP implementation and evaluation
- optional: TensorFlow (if CNN or other deep models explored)
- Jupyter notebook environment used for development and analysis

---

## 📧 Contact

For questions or collaboration opportunities:  
📧 Email: amzarrinnezhad@gmail.com  
💬 GitHub Issues: [here](https://github.com/zamirmehdi/Bioinformatics-Course/issues)  
🌐 GitHub: [@zamirmehdi](https://github.com/zamirmehdi)

---

<div align="center">  
⭐ **If you found this project helpful, please consider giving it a star!** ⭐  
*Amirmehdi Zarrinnezhad*  
</div>
