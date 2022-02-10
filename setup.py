from setuptools import find_packages, setup
setup(
    name='sabs_kabsch',
    packages=find_packages(include=['sabs_kabsch']),
    version='0.1.0',
    description='Implementation of the Kabsch algorithm as a Python library',
    author='Gemma Gordon',
    license='MIT',
    install_requires=['biopython==1.79', 'matplotlib==3.5.1', 'numpy==1.22.2']
)