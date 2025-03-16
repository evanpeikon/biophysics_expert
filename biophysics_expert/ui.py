"""
UI components for the BiophysicsExpert package.

This module provides interactive UI components for use in Jupyter notebooks
and Google Colab.
"""

import ipywidgets as widgets
from IPython.display import display, HTML, Math
import re

from .expert import BiophysicsExpert
from .utils import get_css, get_mathjax_script, get_match_color

def launch_ui(expert=None):
    """
    Launch the interactive UI for the BiophysicsExpert system.
    
    Parameters
    ----------
    expert : BiophysicsExpert, optional
        Expert system instance to use
        
    Returns
    -------
    None
        The UI is displayed in the notebook
    """
    # Create expert system if not provided
    if expert is None:
        expert = BiophysicsExpert()
    
    # Set up the notebook environment
    display(HTML(get_mathjax_script() + get_css()))
    
    # Display title
    display(HTML("""
    <div style="background-color: #f5f5f5; padding: 20px; border-radius: 10px; margin-bottom: 20px;">
        <h1 style="color: #1a73e8;">📚 Biophysics Formula Expert System</h1>
        <p>Find thermodynamic formulas based on variables you know or want to calculate.</p>
    </div>
    """))
    
    # Create tabs
    tab_names = ['Search by Known Variables', 'Search by Unknown', 'View All Formulas']
    tabs = widgets.Tab()
    tab_contents = [widgets.Output() for _ in range(len(tab_names))]
    tabs.children = tab_contents
    for i, name in enumerate(tab_names):
        tabs.set_title(i, name)
    display(tabs)
    
    # Get all unique variables and unknowns
    all_vars = expert.get_all_variables()
    all_unknowns = expert.get_all_unknowns()
    
    # Tab 1: Search by Known Variables
    with tab_contents[0]:
        display(HTML("<h3>Select Known Variables</h3><p>Hold Ctrl/Cmd to select multiple variables</p>"))
        
        # Create multi-select for known variables with improved visibility
        known_vars_select = widgets.SelectMultiple(
            options=sorted(list(all_vars)),
            description='',
            layout=widgets.Layout(width='70%', height='200px')
        )
        display(known_vars_select)
        
        # Create slider for match percentage
        min_match_slider = widgets.IntSlider(
            value=20,
            min=0,
            max=100,
            step=5,
            description='Min Match %:',
            layout=widgets.Layout(width='70%')
        )
        display(min_match_slider)
        
        # Create search button
        search_button = widgets.Button(
            description='Find Matching Formulas',
            button_style='primary',
            layout=widgets.Layout(width='50%', margin='10px 0')
        )
        display(search_button)
        
        # Create output area for results
        output_knowns = widgets.Output()
        display(output_knowns)
        
        # Define button click handler
        def search_by_knowns(b):
            with output_knowns:
                output_knowns.clear_output()
                
                if not known_vars_select.value:
                    display(HTML("<div style='color: #d32f2f; padding: 10px;'>Please select at least one variable.</div>"))
                    return
                
                # Show selected variables for confirmation
                selected_vars = list(known_vars_select.value)
                display(HTML(f"<div style='background-color: #e8f0fe; padding: 10px; border-radius: 5px; margin-bottom: 10px;'><strong>Selected variables:</strong> {', '.join(selected_vars)}</div>"))
                
                matches = expert.find_by_knowns(
                    selected_vars,
                    min_match_ratio=min_match_slider.value/100
                )
                
                if not matches:
                    display(HTML("<div style='color: #d32f2f; padding: 10px;'>No formulas match your criteria. Try selecting different variables or lowering the match percentage.</div>"))
                else:
                    display(HTML(f"<div style='color: #1e8e3e; padding: 10px;'>Found {len(matches)} matching formula(s).</div>"))
                    for formula in matches:
                        match_percent = formula.get("match_ratio", 0) * 100
                        match_color = get_match_color(formula.get("match_ratio", 0))
                        
                        # Create individual output for each formula card
                        formula_output = widgets.Output()
                        display(formula_output)
                        
                        with formula_output:
                            # First show the formula name and container
                            display(HTML(f"""
                            <div class="result-card">
                                <h3>{formula['Formula_Name']}</h3>
                                <div class="formula-box">
                                    <strong>Solves for:</strong> {formula['Unknown']}<br>
                                    <strong>Formula:</strong>
                                </div>
                            """))
                            
                            # Render the LaTeX formula separately
                            display(Math(formula['LaTeX']))
                            
                            # Complete the card with match info and additional info
                            display(HTML(f"""
                                <p><strong>Match:</strong> <span style='color: {match_color};'>{match_percent:.1f}%</span></p>
                                <p><strong>Known:</strong> {', '.join(sorted(formula['common_vars']))}</p>
                                <p><strong>Missing:</strong> {', '.join(sorted(formula['missing_vars']))}</p>
                                {f"<p><strong>Note:</strong> {formula['Additional_Info']}</p>" if formula.get('Additional_Info') else ""}
                            </div>
                            """))
        
        search_button.on_click(search_by_knowns)
    
    # Tab 2: Search by Unknown
    with tab_contents[1]:
        # Create dropdown for unknown variables
        unknown_dropdown = widgets.Dropdown(
            options=sorted(list(all_unknowns)),
            description='Unknown:',
            layout=widgets.Layout(width='50%')
        )
        display(unknown_dropdown)
        
        # Create search button
        search_unknown_button = widgets.Button(
            description='Find Formulas',
            button_style='primary',
            layout=widgets.Layout(width='30%', margin='10px 0')
        )
        display(search_unknown_button)
        
        # Create output area for results
        output_unknown = widgets.Output()
        display(output_unknown)
        
        # Define button click handler
        def search_by_unknown(b):
            with output_unknown:
                output_unknown.clear_output()
                
                matches = expert.find_by_unknown(unknown_dropdown.value)
                
                if not matches:
                    display(HTML(f"<div style='color: #d32f2f; padding: 10px;'>No formulas found that solve for {unknown_dropdown.value}.</div>"))
                else:
                    display(HTML(f"<div style='color: #1e8e3e; padding: 10px;'>Found {len(matches)} formula(s) for calculating {unknown_dropdown.value}.</div>"))
                    
                    for formula in matches:
                        # Create individual output for each formula card
                        formula_output = widgets.Output()
                        display(formula_output)
                        
                        with formula_output:
                            # First show the formula name and container
                            display(HTML(f"""
                            <div class="result-card">
                                <h3>{formula['Formula_Name']}</h3>
                                <div class="formula-box">
                                    <strong>Solves for:</strong> {formula['Unknown']}<br>
                                    <strong>Formula:</strong>
                                </div>
                            """))
                            
                            # Render the LaTeX formula separately
                            display(Math(formula['LaTeX']))
                            
                            # Complete the card with variables and additional info
                            display(HTML(f"""
                                <p><strong>Required variables:</strong> {', '.join(formula['Variables'])}</p>
                                {f"<p><strong>Note:</strong> {formula['Additional_Info']}</p>" if formula.get('Additional_Info') else ""}
                            </div>
                            """))
        
        search_unknown_button.on_click(search_by_unknown)
    
    # Tab 3: View All Formulas
    with tab_contents[2]:
        # Group formulas by name
        formula_groups = {}
        for formula in expert.get_all_formulas():
            name = formula["Formula_Name"]
            if name not in formula_groups:
                formula_groups[name] = []
            formula_groups[name].append(formula)
        
        # Display all formulas
        for name, formulas in formula_groups.items():
            display(HTML(f"<h2>{name}</h2>"))
            
            for formula in formulas:
                # Create individual output for each formula card
                formula_output = widgets.Output()
                display(formula_output)
                
                with formula_output:
                    # First show the formula name and container
                    display(HTML(f"""
                    <div class="result-card">
                        <div class="formula-box">
                            <strong>Solves for:</strong> {formula['Unknown']}<br>
                            <strong>Formula:</strong>
                        </div>
                    """))
                    
                    # Render the LaTeX formula separately
                    display(Math(formula['LaTeX']))
                    
                    # Complete the card with variables and additional info
                    display(HTML(f"""
                        <p><strong>Required variables:</strong> {', '.join(formula['Variables'])}</p>
                        {f"<p><strong>Note:</strong> {formula['Additional_Info']}</p>" if formula.get('Additional_Info') else ""}
                    </div>
                    """))
    
    # Instructions
    display(HTML("""
    <div style="background-color: #e8f0fe; padding: 15px; border-radius: 8px; margin-top: 20px;">
        <h3 style="color: #1a73e8;">Instructions</h3>
        <h4>Search by Known Variables</h4>
        <ol>
            <li><strong>Select multiple variables</strong> by holding down Ctrl (or Cmd on Mac)
