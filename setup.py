from setuptools import setup, find_packages

with open("requeriments.txt") as f:
    requirements = f.read().splitlines()

setup(
    name="comparate",
    version="1.0.1",
    packages=find_packages(include=["comparate", "comparate.*"]),
    include_package_data=True,
    install_requires=requirements,
    entry_points={
        'console_scripts': [
            'comparate=comparate.cli:main',
        ],
    },
)