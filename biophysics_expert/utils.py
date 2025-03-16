"""
Utility functions for the BiophysicsExpert package.

This module provides helper functions for the package.
"""

def get_css():
    """
    Get the CSS styles for the UI.
    
    Returns
    -------
    str
        CSS styles as a string
    """
    return """
    <style>
        .formula-box {
            background-color: #f8f9fa;
            border-radius: 5px;
            padding: 15px;
            margin: 10px 0;
            border-left: 3px solid #4CAF50;
        }
        .result-card {
            background-color: white;
            border-radius: 8px;
            padding: 15px;
            margin-bottom: 15px;
            box-shadow: 0 1px 3px rgba(0,0,0,0.12), 0 1px 2px rgba(0,0,0,0.24);
        }
        .match-high { color: #1e8e3e; }
        .match-medium { color: #fbbc04; }
        .match-low { color: #ea4335; }
    </style>
    """

def get_mathjax_script():
    """
    Get the MathJax script for LaTeX rendering.
    
    Returns
    -------
    str
        MathJax script as a string
    """
    return """
    <script src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.7/MathJax.js?config=TeX-MML-AM_CHTML" async></script>
    <script type="text/x-mathjax-config">
    MathJax.Hub.Config({
      tex2jax: {
        inlineMath: [['$','$'], ['\\(','\\)']],
        displayMath: [['$$','$$'], ['\\[','\\]']],
        processEscapes: true
      },
      CommonHTML: { linebreaks: { automatic: true } },
      "HTML-CSS": {
        linebreaks: { automatic: true },
        availableFonts: ["TeX"],
        scale: 110
      }
    });
    </script>
    """

def get_match_color(match_ratio):
    """
    Get the color for a match ratio.
    
    Parameters
    ----------
    match_ratio : float
        Match ratio between 0 and 1
        
    Returns
    -------
    str
        CSS color code
    """
    match_percent = match_ratio * 100
    
    if match_percent >= 70:
        return "#1e8e3e"  # green for high match
    elif match_percent >= 40:
        return "#fbbc04"  # yellow for medium match
    else:
        return "#ea4335"  # red for low match
