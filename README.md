# Bioinformatics Course - Complete Project Collection

A comprehensive collection of **6 bioinformatics projects** covering fundamental algorithms, sequence analysis, phylogenetics, and machine learning applications in computational biology. Developed as part of the Bioinformatics course at Amirkabir University of Technology.

[![Python](https://img.shields.io/badge/Python-3.x-blue.svg)](https://www.python.org/)
[![Jupyter](https://img.shields.io/badge/Jupyter-Notebook-orange.svg)](https://jupyter.org/)
[![Course](https://img.shields.io/badge/Course-Bioinformatics-blue.svg)](#)
[![University](https://img.shields.io/badge/University-AUT-red.svg)](https://aut.ac.ir/en)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

<div align="center">

**[View Projects](#-projects) • [Technologies](#%EF%B8%8F-technologies-used) • [Getting Started](#-getting-started) • [Contact](#-contact)**

</div>

---

## 📋 Overview

This repository contains a complete journey through **computational biology and bioinformatics**, from basic biological concepts to advanced machine learning applications for viral genome classification. Each project builds upon previous knowledge, creating a comprehensive learning path in modern bioinformatics.

**Course Information:**
- **Institution**: Amirkabir University of Technology (Tehran Polytechnic) - Spring 2022
- **Author**: Amirmehdi Zarrinnezhad  

**What You'll Find:**
- ✅ **5 Theoretical & Programming Assignments** + **1 Final Machine Learning Project**
- ✅ **Complete implementations** with detailed documentation
- ✅ **Step-by-step explanations** of algorithms and methodologies
- ✅ **Real biological data** analysis and interpretation
- ✅ **Production-ready code** with comprehensive README files

---

## 🧬 Projects

### [1. Basic Biology](./1-%20Basic%20biology)
**Type:** Theoretical Analysis  
**Topics:** Stem cells, Gene expression, DNA function, Vestigial traits

A foundational exploration of molecular biology concepts essential for understanding bioinformatics algorithms.

**Key Concepts:**
- Stem cell biology and differentiation
- Gene diversity vs genetic variability
- Central Dogma: DNA → RNA → Protein
- Vestigial traits as evolutionary evidence
- Genotype-phenotype relationships

**Deliverables:**
- Comprehensive theoretical report
- Critical analysis of biological systems
- Evolutionary biology examples

---

### [2. Pairwise Sequence Alignment](./2-%20Pairwise%20Sequence%20Alignment)
**Type:** Programming + Theoretical  
**Topics:** Semi-global alignment, Dynamic programming, Scoring matrices

Implementation of sequence alignment algorithms for comparing protein sequences.

**Key Features:**
- ✅ Semi-global alignment algorithm (Needleman-Wunsch variant)
- ✅ PAM250 substitution matrix
- ✅ Multiple optimal alignment detection
- ✅ Manual calculations (Needleman-Wunsch, Smith-Waterman)

**Technologies:** Python (pure implementation, no libraries)

**Highlights:**
```python
# Finds ALL optimal alignments
aligned_seqs = semi_global_alignment(seq1, seq2, PAM250, gap_penalty=-9)

# Example output:
# Score: 20
# HEAGAWGHE-
# ---PAW-HEA
```

---

### [3. Multiple Sequence Alignment & Database Search](./3-%20Multiple%20Sequence%20Alignment%20-%20DB%20Search)
**Type:** Programming + Theoretical  
**Topics:** Star alignment, Block-based refinement, FASTA, BLAST

Advanced MSA implementation with iterative improvement and database search analysis.

**Key Features:**
- ✅ Star alignment algorithm (center-star heuristic)
- ✅ Block-based iterative refinement
- ✅ Sum-of-pairs scoring
- ✅ FASTA vs BLAST comparison
- ✅ Algorithm complexity analysis

**Technologies:** Python

**Performance:**
```
Initial Score: 51
Final Score:   60  (+17% improvement)
Method:        Iterative block realignment
Convergence:   Automatic detection
```

**Theoretical Topics:**
- FASTA k-tuple matching
- BLAST word search trees
- Heuristic vs exhaustive search

---

### [4. Profile & Hidden Markov Model](./4-%20Profile%20-%20Hidden%20Markov%20model)
**Type:** Programming + Theoretical  
**Topics:** PSSM, Profile search, PSI-BLAST, HMMs

Profile-based sequence search with pseudocount smoothing and HMM analysis.

**Key Features:**
- ✅ Profile construction from MSA
- ✅ Position-Specific Scoring Matrix (PSSM)
- ✅ Log-odds scoring with pseudocount
- ✅ Subsequence search with gap insertion
- ✅ PSI-BLAST mechanism analysis

**Technologies:** Python, math library

**Algorithm:**
```python
# Build profile from MSA
profile = build_profile(msa, pseudocount=2)

# Search query sequence
best_match = search_with_profile(query, profile)
# Output: H-L-P (with optimal gap placement)
```

**Theoretical Topics:**
- PSI-BLAST and profile drift
- Forward algorithm (HMM)
- Viterbi algorithm (optimal path)
- Sequence logos

---

### [5. Phylogenetic Trees](./5-%20Phylogenetic%20Trees)
**Type:** Theoretical Analysis  
**Topics:** UPGMA, Neighbor-Joining, Parsimony, Maximum Likelihood

Comprehensive exploration of phylogenetic tree construction methods.

**Key Topics:**
- ✅ Distance-based methods (UPGMA, NJ)
- ✅ Character-based methods (Parsimony)
- ✅ Probabilistic methods (Maximum Likelihood)
- ✅ Tree comparison and evaluation
- ✅ Algorithm complexity analysis

**Methods Compared:**

| Method | Speed | Accuracy | Molecular Clock |
|--------|-------|----------|-----------------|
| UPGMA | Fast | Moderate | Required |
| NJ | Fast | Good | Not required |
| Parsimony | Slow | Good | Not required |
| ML | Very Slow | Best | Flexible |

**Analysis Includes:**
- Manual UPGMA tree construction
- Neighbor-Joining with corrected distances
- Parsimony scoring (Fitch's algorithm)
- ML probability calculations
- Exhaustive search vs branch-and-bound

---

### [6. Virus Classification (Final Project)](./Virus%20Classification%20(Final%20Project))
**Type:** Machine Learning Project  
**Topics:** K-mer features, Neural networks, Multi-class classification

State-of-the-art virus genome classification using deep learning.

**Problem:**
- Classify DNA sequences into 6 virus classes
- Variable-length sequences (hundreds to thousands of bp)
- Small dataset (1,320 training samples)

**Solution:**
```
DNA Sequence → K-mer Extraction (k=2) → 16 Features → MLP → Class (1-6)
```

**Model Architecture:**
```
Input (16 features)
    ↓
Hidden Layer 1 (64 neurons, ReLU)
    ↓
Hidden Layer 2 (64 neurons, ReLU)
    ↓
Hidden Layer 3 (64 neurons, ReLU)
    ↓
Output (6 classes, Softmax)
```

**Technologies:** Python, scikit-learn, pandas, numpy

**Performance:**
- ✅ **100% accuracy** on development set (180 samples)
- ✅ **~100% accuracy** on test set (400 samples)
- ✅ **Top 10%** performance in class
- ✅ **CPU training** (<5 minutes)

**Key Innovation:**
- Optimal k-mer size (k=2) through systematic experimentation
- Length-normalized features for fair comparison
- Balanced architecture preventing overfitting

---

## 🛠️ Technologies Used

### Programming Languages
- **Python 3.x** - Primary language for all implementations

### Core Libraries

**Data Processing:**
- `pandas` - Data manipulation and CSV handling
- `numpy` - Numerical computations
- `itertools` - Combinatorial operations (k-mer generation)

**Machine Learning:**
- `scikit-learn` - MLP classifier, metrics, preprocessing
- `tensorflow/keras` - Alternative deep learning experiments

**Bioinformatics:**
- Custom implementations (no external bio libraries)
- Pure Python algorithms for educational purposes

**Visualization:**
- `matplotlib` - Learning curves and performance plots

### Development Tools
- **Jupyter Notebook** - Interactive development and documentation
- **Git** - Version control
- **Quera** - Submission and evaluation platform

---

## 📊 Learning Path

```
Project 1: Basic Biology
    ↓ (Understand biological foundations)
Project 2: Pairwise Alignment
    ↓ (Dynamic programming, scoring schemes)
Project 3: Multiple Alignment
    ↓ (Heuristic algorithms, database search)
Project 4: Profile & HMM
    ↓ (Position-specific scoring, probabilistic models)
Project 5: Phylogenetic Trees
    ↓ (Evolutionary relationships, tree construction)
Project 6: Machine Learning
    ↓ (Deep learning for genome classification)

Complete Bioinformatics Pipeline!
```

**Skills Progression:**
1. **Theoretical foundations** → Biological understanding
2. **Algorithm implementation** → Dynamic programming
3. **Heuristic methods** → Speed vs accuracy trade-offs
4. **Probabilistic models** → HMMs, likelihood
5. **Phylogenetic analysis** → Evolutionary inference
6. **Machine learning** → Modern AI applications

---

## 🚀 Getting Started

### Prerequisites

```bash
Python 3.7 or higher
pip (Python package manager)
```

### Installation

**Clone the repository:**
```bash
git clone https://github.com/zamirmehdi/Bioinformatics-Course.git
cd Bioinformatics-Course
```

**Install dependencies:**
```bash
pip install pandas numpy scikit-learn matplotlib jupyter
```

### Running Projects

**For Python scripts:**
```bash
cd "2- Pairwise Sequence Alignment/src"
python semi_global_alignment.py < input.txt
```

**For Jupyter notebooks:**
```bash
cd "Virus Classification (Final Project)/src"
jupyter notebook BioInformatics_FinalProject.ipynb
```

**For theoretical projects:**
- Navigate to project folder
- Review `README.md` for detailed explanations
- Check `Report.pdf` for solutions (Persian)

---

## 📁 Repository Structure

```
Bioinformatics-Course/
│
├── 1- Basic biology/
│   ├── Instruction.pdf
│   ├── Report.pdf
│   └── README.md
│
├── 2- Pairwise Sequence Alignment/
│   ├── docs/
│   │   ├── Programming Instruction.pdf
│   │   └── Theoretical/
│   │       ├── Instruction.pdf
│   │       └── Report.pdf
│   ├── src/
│   │   └── semi_global_alignment.py
│   └── README.md
│
├── 3- Multiple Sequence Alignment - DB Search/
│   ├── docs/
│   │   ├── Programming Instruction MSA.pdf
│   │   └── Theoretical/
│   │       ├── Instruction.pdf
│   │       ├── Report.pdf
│   │       └── cstar.pdf
│   ├── src/
│   │   └── main.py
│   └── README.md
│
├── 4- Profile - Hidden Markov model/
│   ├── docs/
│   │   ├── Programming Instruction - Profile.pdf
│   │   └── Theoretical/
│   │       ├── Instruction.pdf
│   │       └── Report.pdf
│   ├── src/
│   │   └── Profile.py
│   └── README.md
│
├── 5- Phylogenetic Trees/
│   ├── Instruction.pdf
│   ├── Report.pdf
│   └── README.md
│
├── Virus Classification (Final Project)/
│   ├── data/
│   │   ├── training_set.csv
│   │   ├── development_set.csv
│   │   └── test_set.csv
│   ├── docs/
│   │   ├── Instruction.pdf
│   │   └── Report.pdf
│   ├── src/
│   │   └── BioInformatics_FinalProject.ipynb
│   ├── README.md (Part 1)
│   └── README_PART2.md (Part 2)
│
└── README.md (This file)
```
---

## 📊 Quick Links

| Project | Type | Status | README |
|---------|------|--------|--------|
| 1. Basic Biology | Theory | ✅ Complete | [View](./1-%20Basic%20biology/README.md) |
| 2. Sequence Alignment | Code + Theory | ✅ Complete | [View](./2-%20Pairwise%20Sequence%20Alignment/README.md) |
| 3. MSA & DB Search | Code + Theory | ✅ Complete | [View](./3-%20Multiple%20Sequence%20Alignment%20-%20DB%20Search/README.md) |
| 4. Profile & HMM | Code + Theory | ✅ Complete | [View](./4-%20Profile%20-%20Hidden%20Markov%20model/README.md) |
| 5. Phylogenetic Trees | Theory | ✅ Complete | [View](./5-%20Phylogenetic%20Trees/README.md) |
| 6. Virus Classification | ML Project | ✅ Complete | [View](./Virus%20Classification%20(Final%20Project)/README.md) |

---

## 🎓 Key Learning Outcomes

After completing these projects, you will be able to:

### Algorithms & Data Structures
✅ Implement dynamic programming for sequence alignment  
✅ Design heuristic algorithms for NP-hard problems  
✅ Optimize time and space complexity  
✅ Handle variable-length biological data  

### Bioinformatics
✅ Perform pairwise and multiple sequence alignment  
✅ Search biological databases efficiently  
✅ Build and use profiles for sequence search  
✅ Construct phylogenetic trees  
✅ Apply machine learning to genomic data  

### Machine Learning
✅ Extract features from biological sequences  
✅ Design neural network architectures  
✅ Prevent overfitting on small datasets  
✅ Evaluate models with proper metrics  
✅ Tune hyperparameters systematically  

### Software Engineering
✅ Write clean, documented, maintainable code  
✅ Structure projects professionally  
✅ Create comprehensive documentation  
✅ Use version control (Git)  
✅ Follow best practices  

---

## 📈 Project Statistics

| Metric | Value |
|--------|-------|
| Total Projects | 6 (5 assignments + 1 final) |
| Lines of Code | ~2,000+ |
| Programming Projects | 4 |
| Theoretical Projects | 2 |
| Algorithms Implemented | 15+ |
| Documentation Pages | 100+ (combined READMEs) |
| Test Cases Passed | 100% |
| Final Grade | Excellent (Top 10%) |

---

## 🏆 Highlights & Achievements

### Technical Achievements
- ✅ **100% test case success** across all programming projects
- ✅ **Near-perfect ML model** (100% dev, ~100% test accuracy)
- ✅ **Efficient implementations** (CPU-only, fast execution)
- ✅ **Multiple optimal solutions** (Project 2 - finds ALL alignments)
- ✅ **Iterative refinement** (Project 3 - automatic improvement)

### Documentation Quality
- ✅ **Comprehensive READMEs** for every project
- ✅ **Step-by-step explanations** with examples
- ✅ **Visual diagrams** and algorithm illustrations
- ✅ **Code comments** in English
- ✅ **Bilingual support** (English docs, Persian reports)

### Academic Impact
- ✅ **Top 10% performance** in final project
- ✅ **Complete assignment portfolio** (6/6 completed)
- ✅ **High-quality reports** with detailed analysis
- ✅ **Reproducible results** with clear instructions

---

## 🔬 Real-World Applications

### Medical & Clinical
- **Disease diagnosis**: Sequence-based pathogen identification
- **Personalized medicine**: Genetic variant analysis
- **Drug discovery**: Protein target identification
- **Epidemiology**: Outbreak tracking and surveillance

### Research & Academia
- **Evolutionary biology**: Phylogenetic studies
- **Comparative genomics**: Cross-species analysis
- **Protein function**: Homology-based prediction
- **Gene discovery**: Novel sequence identification

### Biotechnology
- **Genetic engineering**: CRISPR guide design
- **Synthetic biology**: Sequence optimization
- **Bioinformatics tools**: Algorithm development
- **Data analysis**: High-throughput sequencing

---

## 📚 References & Resources

### Key Papers
- Needleman & Wunsch (1970) - Global alignment algorithm
- Smith & Waterman (1981) - Local alignment algorithm
- Henikoff & Henikoff (1992) - BLOSUM matrices
- Altschul et al. (1990) - BLAST algorithm
- Eddy (1998) - Profile HMMs
- Hemalatha Gunasekaran et al. (2021) - [Analysis of DNA Sequence Classification Using CNN and Hybrid Models (K-mer Encoding)](https://pubmed.ncbi.nlm.nih.gov/34306171/)

### Online Resources
- [NCBI BLAST](https://blast.ncbi.nlm.nih.gov/) - Database searching
- [UniProt](https://www.uniprot.org/) - Protein sequences
- [Pfam](https://pfam.xfam.org/) - Protein families
- [EMBOSS](https://www.ebi.ac.uk/Tools/emboss/) - Bioinformatics tools

### Textbooks
- *Biological Sequence Analysis* - Durbin et al.
- *Introduction to Computational Molecular Biology* - Setubal & Meidanis
- *Algorithms on Strings, Trees, and Sequences* - Gusfield

---

## 🤝 Contributing

While this is a personal academic repository, contributions are welcome!

**Ways to contribute:**
- 🐛 Report bugs or issues
- 💡 Suggest improvements
- 📖 Improve documentation
- ✨ Add new features or optimizations

**Please:**
1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Submit a pull request

---
<!--
## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

**Academic Use:**
- ✅ Free to use for educational purposes
- ✅ Citation appreciated
- ⚠️ Do not submit as your own coursework

---


## 📧 Contact

**Author:** Amirmehdi Zarrinnezhad  
**University:** Amirkabir University of Technology  
**Email:** amzarrinnezhad@gmail.com  
**GitHub:** [@zamirmehdi](https://github.com/zamirmehdi)

**Questions?**
- 💬 Open an [Issue](https://github.com/zamirmehdi/Bioinformatics-Course/issues)
- 📧 Send an email
- 🌟 Star the repository if you find it helpful!

---
-->
## ℹ️ Project Information

**Author:** Amirmehdi Zarrinnezhad  
**Course:** Bioinformatics  
**University:** Amirkabir University of Technology (Tehran Polytechnic) - Fall 2022  
**Language:** English (README), Persian (Instruction and Report PDFs)  
**GitHub Link:** [Bioinformatics Course](https://github.com/zamirmehdi/Bioinformatics-Course)

<div align="center">

**Bioinformatics Course Projects**

[1: Basic Biology](../1-%20Basic%20biology) | [2: Sequence Alignment](.) | [3: MSA & DB Search](../3-%20Multiple%20Sequence%20Alignment%20-%20DB%20Search) | [4: Profile HMM](../4-%20Profile%20-%20Hidden%20Markov%20model) | [5: Phylogenetic Trees](../5-%20Phylogenetic%20Trees) | [Final: Virus Classification](../Virus%20Classification%20(Final%20Project))

</div>

---

## 📧 Contact

Questions or collaborations? Feel free to reach out!  
📧 Email: amzarrinnezhad@gmail.com  
💬 Open an [Issue](https://github.com/zamirmehdi/Bioinformatics-Course/issues)  
🌐 GitHub: [@zamirmehdi](https://github.com/zamirmehdi)

---

<div align="center">

[⬆ Back to Main Repository](https://github.com/zamirmehdi/Bioinformatics-Course)

</div>

<p align="right">(<a href="#top">back to top</a>)</p>

<div align="center">

⭐ **If you found this project helpful, please consider giving it a star!** ⭐

*Amirmehdi Zarrinnezhad*

</div>

<!--
## 🙏 Acknowledgments

**Course Staff:**
- **Professor**: Bioinformatics Course Instructor - Amirkabir University
- **Teaching Assistants**: 
  - Mohammadreza Nafar - Programming guidance
  - Aref Mousavi - Theoretical support
  - Fateme Zahrasarafaei - General assistance

**Tools & Platforms:**
- **Quera** - Assignment submission and evaluation
- **scikit-learn** - Machine learning framework
- **Python Community** - Excellent documentation and support

**Special Thanks:**
- Classmates for discussions and collaboration (within academic integrity guidelines)
- Open-source community for tools and libraries

---
-->

<!--
<div align="center">

## ⭐ Star History

If you found this repository helpful, please consider giving it a star!

[![Star History Chart](https://api.star-history.com/svg?repos=zamirmehdi/Bioinformatics-Course&type=Date)](https://star-history.com/#zamirmehdi/Bioinformatics-Course&Date)

---

**Made with ❤️ and lots of ☕ by Amirmehdi Zarrinnezhad**


[⬆ Back to Top](#bioinformatics-course---complete-project-collection)

</div> -->
