[![Documentation Status](https://readthedocs.org/projects/sabs-kabsch/badge/?version=latest&style=flat-square)](https://sabs-kabsch.readthedocs.io/en/latest/?badge=latest&style=flat-square)
# sabs-kabsch
Implementation of Kabsch algorithm as a Python library for SABS R3 Modelling &amp; Scientific Computing. 

## *Instructions for use:*
 
- **Installation**

    - Activate Python environment (*venv*)
    - Package is pip-installable:

        ```pip install sabs-kabsch```
        

- **Examples**

    - Navigate to *examples* folder to find:

        - Jupyter notebook explaining how the algorithm implementation works step-by-step.
        - Example matplotlib visualisations of results.


## *About the Kabsch algorithm:*

- based on linear algebra
- alignment algorithm 
- minimises the RMSD between two sets of points - in this case in 3D space but it can be scaled up to N dimensions
- Kabsch algorithm based on finding the optimal rotation to align two sets of points
- O(n) efficiency 


### *What can the Kabsch algorithm be used for?*

- applications in computational biology and chemistry for aligning molecules e.g. protein structures
- e.g. pymol align method? 

### *How does the Kabsch algorithm work?*

  <img src="https://github.com/gemmagordon/sabs-kabsch/blob/main/examples/A%20vs%20B%20pre-Kabsch.png" width="40%" /> <img src="https://github.com/gemmagordon/sabs-kabsch/blob/main/examples/A%20vs%20B%20post-Kabsch.png" width="40%" />

1) Translation
2) Computing a covariance matrix
3) Computing the optimal rotation matrix

### *Remaining limitations of algorithm and this implementation:*
- doesn't include scaling
- limited by necessity of matrices being same length
- rigidity in transformations? 
- other registration/alignment algorithms - e.g. ICP? 
- using SVD doesn't account for special cases e.g. matrices with no inverse? 

### *Sources*
- https://zpl.fi/aligning-point-patterns-with-kabsch-umeyama-algorithm/
- https://pymolwiki.org/index.php/Kabsch#:~:text=The%20Kabsch%20algorithm%20uses%20linear,implementation%20of%20the%20Kabsch%20algorithm.
- https://stackoverflow.com/questions/60877274/optimal-rotation-in-3d-with-kabsch-algorithm
- https://towardsdatascience.com/the-definitive-procedure-for-aligning-two-sets-of-3d-points-with-the-kabsch-algorithm-a7ec2126c87e
- http://computerandchemistry.blogspot.com/2013/04/calculate-rmsd-from-two-xyz-files.html
- https://www.hgc.jp/~tshibuya/papers/jcb10b_preprint(faster3D).pdf 