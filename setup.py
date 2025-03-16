from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="biophysics_expert",
    version="0.1.0",
    author="Your Name",
    author_email="your.email@example.com",
    description="A biophysics formula expert system for educational purposes",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/biophysics_expert",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
    install_requires=[
        "ipywidgets>=7.0.0",
        "ipython>=7.0.0",
    ],
    include_package_data=True,
)
