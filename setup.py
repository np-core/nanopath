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
        "python-dateutil",
        "numpy",
        "pyfastx"
    ],
    entry_points="""
        [console_scripts]
        np=nanopath.terminal.client:terminal_client
        nanopath=nanopath.terminal.client:terminal_client
    """,
    version="0.1",
    license="MIT",
    description="NanoPath package for nanopore sequencing in pathology units",
)
