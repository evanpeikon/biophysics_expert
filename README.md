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
