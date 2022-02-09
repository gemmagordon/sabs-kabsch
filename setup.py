from setuptools import find_packages, setup
setup(
    name='sabs_kabsch',
    packages=find_packages(),
    version='0.1.0',
    description='Implementation of the Kabsch algorithm as a Python library',
    author='Gemma Gordon',
    license='MIT',
    install_requires=['numpy==1.22.1', 'biopython==1.79']
)