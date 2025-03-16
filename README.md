# BiophysicsExpert

A Python package for biophysics formula lookup and exploration.

## Description

BiophysicsExpert is an educational tool that provides a searchable database of biophysics formulas. It allows users to find formulas based on known variables or desired unknown variables, making it easier to identify the right equations for specific problems.

## Features

- Search formulas by known variables
- Find formulas that solve for specific unknowns
- View all available formulas
- Interactive UI for Jupyter notebooks and Google Colab
- LaTeX rendering of mathematical formulas

## Installation

You can install the package directly from Google Colab:

```python
!pip install git+https://github.com/evanpeikon/biophysics_expert.git
```

## Usage

### Basic Usage

```python
from biophysics_expert import BiophysicsExpert

# Create an expert system
expert = BiophysicsExpert()

# Find formulas by known variables
results = expert.find_by_knowns(["δQ", "δW"], min_match_ratio=0.5)
for formula in results:
    print(f"Formula: {formula['Formula_Name']}")
    print(f"Equation: {formula['Formula']}")
    print(f"Solves for: {formula['Unknown']}")
    print("---")

# Find formulas by unknown variable
results = expert.find_by_unknown("dE")
for formula in results:
    print(f"Formula: {formula['Formula_Name']}")
    print(f"Required variables: {', '.join(formula['Variables'])}")
```

### Interactive UI in Notebooks

```python
from biophysics_expert.ui import launch_ui

# Launch the interactive UI
launch_ui()
```

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License - see the LICENSE file for details.
