# sabs_kabsch
Implementation of Kabsch algorithm as a Python library for SABS R3 Modelling &amp; Scientific Computing. 

## *Instructions for use:*
 
- **Installation**

    - Package is pip-installable: <br/><br/>
    First, clone the repository:  

      ```git clone https://github.com/gemmagordon/sabs_kabsch.git```

      Create a new Python virtual environment, activate it & install the package requirements (*requirements.txt*).<br/>

      Then in the cloned directory, run: 
      
      ```pip install .```<br/>

      To use in a Python script, import the package with:
        
      ```import sabs_kabsch``` 

      Call functions from sabs_kabsch with, for example:

       ```sabs_kabsch.generate_matrices()```

- **Example usage** 

    - Navigate to **examples** folder to find:

        - A **Jupyter notebook** (*kabsch.ipynb*) explaining how the algorithm implementation works and how to use it, using each function available in the package step-by-step.
        - Examples of matplotlib visualisations of results, from running the example datasets in *kabsch.ipynb* or *kabsch_functions.py*


## *About the Kabsch algorithm:*

The Kabsch algorithm is an alignment algorithm which aims to find an optimal rotation (where optimal is defined by minimising the RMSD (root-mean-square deviation), a measure of difference) to align two sets of points with O(n) efficiency.

### *What can the Kabsch algorithm be used for?*

The Kabsch algorithm is applicable in computational biology, for example in aligning two protein structures and calculating the minimal RMSD between them. 


### *How does the Kabsch algorithm work?*

  <img src="https://github.com/gemmagordon/sabs-kabsch/blob/main/examples/A%20vs%20B%20pre-Kabsch.png" width="40%" /> <img src="https://github.com/gemmagordon/sabs-kabsch/blob/main/examples/A%20vs%20B%20post-Kabsch.png" width="40%" />

1) Start with two matrices, A & B, of size NxM (which in this implementation are two datasets of X,Y,Z coordinates of size 3x50). Matrix A is being aligned to matrix B, which is acting as a sort of template.

2) Translate both A and B so that their centroids lie around the origin on the 3D axes. To do this, the mean of each column of the matrix (i.e., the mean of the X, Y and Z coordinates) is found and then subtracted from each value (coordinate) in the corresponding column.

2) Compute the covariance matrix for A and B using the formula: 

      *H = B<sub>translated</sub><sup>T</sup> x A<sub>translated</sub>*

3) Computing the optimal rotation matrix using singular value decomposition (SVD) using the formula: 

      *H = U x S x V<sup>T</sup>*

      *R = V<sup>T</sup> x U<sup>T</sup>*


4) Apply this rotation to A, so that A is rotated in such a way that offers the best alignment with B (i.e. the RMSD is minimised). 

      *A<sub>rotated</sub> = A<sub>translated</sub> x R*

5) Translate the resulting matrix for A, after its transformation, back to where the centroid was originally placed. This is done by simply adding back the relevant mean value to each coordinate in the matrix. *The same is carried out for B, in that it is translated to the origin and then translated back in this implementation, but this isn't really necessary, as the result is equal to the original matrix for B (as would be expected).*


### *Remaining limitations of algorithm and/or this implementation:*
- The Kabsch algorithm is also known as 'partial' Procrustes superimposition (as opposed to full PS) as although there is translation and rotation, there is no scaling. 
- Correspondence between the sets of points is assumed - other alignment algorithms such as ICP do not make this assumption.
- This implementation is limited by the size/shapes of the two matrices having to be identical.
- SVD is used here in calculating the optimal rotation matrix but there are alternative methods - this solution might not cover all cases.

### *Sources*
- https://zpl.fi/aligning-point-patterns-with-kabsch-umeyama-algorithm/
- https://pymolwiki.org/index.php/Kabsch#:~:text=The%20Kabsch%20algorithm%20uses%20linear,implementation%20of%20the%20Kabsch%20algorithm.
- https://stackoverflow.com/questions/60877274/optimal-rotation-in-3d-with-kabsch-algorithm
- https://towardsdatascience.com/the-definitive-procedure-for-aligning-two-sets-of-3d-points-with-the-kabsch-algorithm-a7ec2126c87e
- http://computerandchemistry.blogspot.com/2013/04/calculate-rmsd-from-two-xyz-files.html
- https://www.hgc.jp/~tshibuya/papers/jcb10b_preprint(faster3D).pdf 
- https://en.wikipedia.org/wiki/Kabsch_algorithm
- https://zpl.fi/aligning-point-patterns-with-kabsch-umeyama-algorithm/ 
- https://en.wikipedia.org/wiki/Procrustes_analysis

