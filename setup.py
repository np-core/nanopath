from setuptools import setup, find_packages

setup(
    name="nanopath",
    url="https://github.com/np-core/nanopath",
    author="Eike J. Steinig",
    author_email="eikejoachim.steinig@my.jcu.edu.au",
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        "click",
        "pytest",
        "paramiko",
        "pymongo",
        "mongoengine",
        "flask",
        "flask-socketio",
        "flask-cors",
        "tqdm",
        "colorama",
        "pandas",
        "seaborn",
        "scipy",
        "dendropy",
        "python-dateutil",
        "scikit-learn",
        "numpy",
        "jinja2",
        "pyfastx",
        "pysankey",
        "ratelimiter",
        "biopython",
        "networkx",
        "dateparser",
        "pyyaml",
        "phylo-treetime",
        "ont_fast5_api",
        "weasyprint"
    ],
    entry_points="""
        [console_scripts]
        np=nanopath.terminal.cli:terminal_client
        nanopath=nanopath.terminal.cli:terminal_client
    """,
    version="0.1",
    license="MIT",
    description="NanoPath package for nanopore sequencing in pathology units",
)
