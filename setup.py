from setuptools import setup, find_packages

def readfile(filename):
    with open(filename, 'r+') as f:
        return f.read()


setup(
    name="genetargeter",
    version="2020.9.21",
    description="Creates constructs to gene edit P. falciparum; developed by the Niles Lab at MIT.",
    long_description=readfile('README.md'),
    author="Pablo Cardenas",
    author_email="pablocarderam@gmail.com",
    url="genetargeter.mit.edu",
    py_modules=['genetargeter','gRNAScores'],
    install_requires=['joblib',
                      'numpy',
                      'scikit-learn',
                      'scipy',
                      'biopython',
                      'pandas',
                      'dill',
                      'hyperopt'
                      ],
    packages=find_packages(),
    license=readfile('LICENSE'),
    include_package_data=True,
    entry_points={
        'console_scripts': [
            'genetargeter = genetargeter:runAll'
        ]
    },
)
