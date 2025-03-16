# ☄️ BiophysicsExpert Overview
BiophysicsExpert is a comprehensive educational tool that provides a searchable database of over 150 biophysics formulas. It's designed for students and educators, enabling users to quickly identify relevant equations by searching with known variables or desired unknown variables. This streamlined approach eliminates the need to manually sift through textbooks or reference materials when solving biophysics problems.

The tool offers flexible usage options, functioning as both a command-line utility or as an interactive graphical interface within Jupyter notebooks or Google Colab environments.

## Key Features
BiophysicsExpert has the following key features:
- Variable-Based Formula Search: Locate formulas based on the variables you already know
- Solution-Oriented Search: Find formulas that can solve for specific unknown variables
- Comprehensive Database: Access over 150 biophysics formulas in a structured format
- Complete Formula Catalog: View and browse all available formulas in the database
- Interactive Notebook Interface: User-friendly GUI specifically designed for Jupyter and Google Colab environments
- Professional Presentation: Mathematical expressions rendered with LaTeX
- Flexible Usage Options: Command-line and GUI interfaces to suit different workflows

# ☄️ Installation and Useage
## Installation
Install BiophysicsExpert directly from GitHub using pip:

```python
# Install the package directly from GitHub (Remove ! when running from Bash command line. Keep ! when running from notebook environment)
!pip install git+https://github.com/evanpeikon/biophysics_expert.git
```

## Useage Method 1 (Commnad Line Interface)
BiophysicsExpert can be imported and used programmatically with simple Python code, as demonstrated below:

```python
# Import BiophysicsExpert
from biophysics_expert import BiophysicsExpert

# Initialize the expert system
expert = BiophysicsExpert()
```

The code above imports the BiophysicsExpert class and creates an instance of it called expert. This instance gives you access to the database of biophysics formulas. Next, I'll demonstrated how to search the database by known variables:

```python
# Search by Known Variables
results = expert.find_by_knowns(["variable 1", "variable 2"], min_match_ratio=0.5)

# Display Results from Known Variables Search
for formula in results:
    print(f"Formula: {formula['Formula_Name']}")
    print(f"Equation: {formula['Formula']}")
    print(f"Solves for: {formula['Unknown']}")
    print("---")
```

This code above the formula database for equations that contain the specified variables "variable 1" and "variable 2". The min_match_ratio=0.5 parameter indicates that formulas need to contain at least 50% of the provided variables to be included in the results. Then, it iterates through each formula found and prints the name of the formula, the actual equation, and what variable the formula can solve for. 

Next, I'll demonstrated how to search the database by an unknown variable:

```python
# Search by Known Variables
results = expert.find_by_unknown("unknown variable")

# Display Results from Unknown Variable Search
for formula in results:
    print(f"Formula: {formula['Formula_Name']}")
    print(f"Required variables: {', '.join(formula['Variables'])}")
```

The code above searches for formulas that can be used to calculate an uknown variable. Then, the code prints each formula that can solve for said unknown variable as well as the name of the formula and a comma-separated list of all variables needed to use this formula. 

### Example Usage
```Bash
$ pip install git+https://github.com/evanpeikon/biophysics_expert.git
```
```python
# Import BiophysicsExpert
from biophysics_expert import BiophysicsExpert

# Initialize the expert system
expert = BiophysicsExpert()

# Search by Known Variables
results = expert.find_by_knowns(["dS", "T", "p"], min_match_ratio=0.5)

for formula in results:
    print(f"Formula: {formula['Formula_Name']}")
    print(f"Equation: {formula['Formula']}")
    print(f"Solves for: {formula['Unknown']}")
    print("---")
```
```
Formula: First Law of Thermodynamics for Reversible Processes
Equation: dE = TdS - pdV
Solves for: dE
---
Formula: Entropy Differential Equation For Reversible Processes
Equation: dS = δQ/T
Solves for: dS
---
Formula: Fundamental Gibbs Equation
Equation: dE = TdS - pdV + μdN
Solves for: dE
---
Formula: Heat Transfer and Entropy Change
Equation: dS = δQ/T
Solves for: dS
---
Formula: Heat Transfer and Entropy Change For Constant Temperature Processes
Equation: ΔS = Q/T
Solves for: ΔS
---
Formula: Average Kinetic Energy
Equation: ⟨E_kin⟩ = 3k_BT/2
Solves for: ⟨E_kin⟩
```

## Useage Method 2 (Graphical User Interface)

For a more interactive experience in Jupyter notebooks or Google Colab, use the built-in graphical interface:

```python
# Install the package directly from GitHub 
!pip install git+https://github.com/evanpeikon/biophysics_expert.git

# Import BiophysicsExpert UI and launch UI
from biophysics_expert.ui import launch_ui

# Launch the interactive UI
launch_ui()
```
The code above produces the following output:

<img width="729" alt="Screenshot 2025-03-16 at 4 05 04 PM" src="https://github.com/user-attachments/assets/68733fb1-4567-4214-a582-00267a8ecf6d" />

This app-based interface provides intuitive search fields, dropdown menus, and visual formula displays with LaTeX rendering for a seamless user experience.

# ☄️ Contributors
## How to Contribute
Contributions to BiophysicsExpert are welcome! Here's how you can help:
- Add New Formulas: Expand our database with additional biophysics formulas
- Improve Documentation: Help clarify usage instructions or add examples
- Fix Bugs: Address issues in the codebase
- Enhance Features: Improve existing functionality or add new features
