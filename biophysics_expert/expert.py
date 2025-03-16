"""
Core functionality of the BiophysicsExpert system.

This module provides the main BiophysicsExpert class that handles formula lookup
and management.
"""

from .formulas import get_default_formulas

class BiophysicsExpert:
    """
    Expert system for biophysics formulas.
    
    Provides methods to find formulas based on known variables or desired
    unknown variables.
    """
    
    def __init__(self, custom_formulas=None):
        """
        Initialize the expert system with formulas.
        
        Parameters
        ----------
        custom_formulas : list, optional
            List of custom formula dictionaries to add to the default set
        """
        # Initialize with default formulas
        self.formulas = get_default_formulas()
        
        # Add custom formulas if provided
        if custom_formulas:
            if isinstance(custom_formulas, list):
                self.formulas.extend(custom_formulas)
            else:
                raise TypeError("custom_formulas must be a list of formula dictionaries")
        
        # Cache all variables and unknowns for faster lookups
        self._cache_variables()
    
    def _cache_variables(self):
        """Cache all unique variables and unknowns for faster lookups."""
        self.all_variables = set()
        self.all_unknowns = set()
        
        for formula in self.formulas:
            self.all_variables.update(formula.get("Variables", []))
            self.all_unknowns.add(formula.get("Unknown", ""))
    
    def add_formula(self, formula):
        """
        Add a new formula to the expert system.
        
        Parameters
        ----------
        formula : dict
            Dictionary with formula information
            
        Returns
        -------
        bool
            True if added successfully, False otherwise
        """
        # Validate formula structure
        required_keys = ["Formula_Name", "Unknown", "Formula", "Variables"]
        if not all(key in formula for key in required_keys):
            return False
        
        # Add formula and update cache
        self.formulas.append(formula)
        self.all_variables.update(formula.get("Variables", []))
        self.all_unknowns.add(formula.get("Unknown", ""))
        
        return True
    
    def find_by_unknown(self, unknown):
        """
        Find formulas that can calculate the given unknown variable.
        
        Parameters
        ----------
        unknown : str
            The variable to solve for
            
        Returns
        -------
        list
            List of matching formula dictionaries
        """
        matches = []
        for formula in self.formulas:
            if formula["Unknown"] == unknown:
                matches.append(formula)
        return matches
    
    def find_by_knowns(self, knowns, min_match_ratio=0.0):
        """
        Find formulas where the known variables are present.
        
        Parameters
        ----------
        knowns : str or list
            Known variable(s)
        min_match_ratio : float, optional
            Minimum ratio of known variables to total variables
            
        Returns
        -------
        list
            List of matching formula dictionaries with additional match information
        """
        if isinstance(knowns, str):
            knowns = [knowns]
        known_set = set(knowns)
        
        matches = []
        for formula in self.formulas:
            formula_vars = set(formula.get("Variables", []))
            
            # Find intersection between knowns and formula variables
            common_vars = known_set.intersection(formula_vars)
            
            # Calculate match ratio
            match_ratio = len(common_vars) / len(formula_vars) if formula_vars else 0
            
            # If we have a minimum match, include this formula
            if match_ratio >= min_match_ratio:
                matches.append({
                    **formula,
                    "match_ratio": match_ratio,
                    "common_vars": common_vars,
                    "missing_vars": formula_vars - known_set
                })
        
        # Sort by match ratio (highest first)
        return sorted(matches, key=lambda x: x["match_ratio"], reverse=True)
    
    def get_all_variables(self):
        """
        Get all unique variables used in formulas.
        
        Returns
        -------
        set
            Set of all variables
        """
        return self.all_variables
    
    def get_all_unknowns(self):
        """
        Get all unique unknowns that can be solved for.
        
        Returns
        -------
        set
            Set of all unknowns
        """
        return self.all_unknowns
    
    def get_all_formulas(self):
        """
        Get all formulas in the expert system.
        
        Returns
        -------
        list
            List of all formula dictionaries
        """
        return self.formulas.copy()
