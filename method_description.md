Method\_description
================

Results
-------

### Reconstructing the transcriptome in 3D.

From a total of nine sections (*n* = 9) we reconstructed the embryonic heart in three-dimensions by registering individual tissue sections onto matched planes within a 3D reference atlas that we constructed based on (). We further deliminated the embrionic heart into 10 unique anatomical regions (LV: right ventricle 26022 cells 406 spots, RV: left ventricle 44416 cells 333 spots, WH: whole heart 19849 cells 201 spots, epc: epicardium 814 cells 8 spots, OT: outflow tract 13622 cells 158 spots, LA: right atria 8045 cells 122 spots, RA: left atria 11471 cells 104 spots, P: pulmonary trunk and arteries 2100 cells 13 spots, A: ascending aorta 968 cells 7 spots).

Using cell segmentation algorithms on the HNE image we detected a total of 127307 cell nuclei in all sections. Out of all the cells detected a total of 23660 (18.58%) were captured by one of the barcoded RNA-capturing spots with an average of 16 nuclei per spot (*M* = 16.13 nuclei, *SD* = 11.56 nuclei).

<!-- A cross the entire tissue there was a negative correlation between number of cell nuclei in each spot and the number of transcripts detected in the same spot _r_ = -0.35 (_t_ = -14.41, _df_ = 1465, _p_ < 0.0001). However, upon closer examination, cell density and number of transcripts seggregated into two seperate bivariate distributions: one with an average of cell nuclei per spots and another distribution with an average of less than ten cell nuclei per spot. These two clusters were located into different anatomical compartments where spots with an average number of cell nuclei .


![Showing spots in 3D.](./images/nuclei_to_count.png)
-->
Methods
-------

### Reference atlas construction.

Our atlas is based on of carnegie stage 18 (de Bakker et al., 2016), CS18-6524. The heart was cropped out from the rest of the atlas by creating a bounding box with the pixel coordinates (7, 170, 775, 1239). Systematic elongations of the heart compared to the atlas because of its cutting plane was corrected by extending the axis by bicubic interpolation with a factor of 3.25. The atlas ontology was further refined by drawing region of interests for regions which was then extended to 10 unique anatomical regions. Each anatomical region is represented in the 2D plane as vector graphics.

### Image processing

Image processing for atlas registration was done using custom code written in C++ and R using OpenCV.

### Segmentation.

H&E images were converted to grayscale and tissue outline as well as cell nuclei were segmented by running eight linearly spaced binary thresholds on the intensity and then using a connected components algorithm with subsequent filtering of components based on size, eccentricity and intensity (Fürth et al., 2018).

### Registration.

Registration of the heart atlas onto each section follows a similar scheme as previously described for brain tissue (Fürth et al., 2018). The contour of both the H&E labeled tissue section and the atlas section were reduced and matched by use of principal component analysis. Any registration errors were manually corrected by changing, adding or removing landmarks along the contour. The registration transform was obtained by thin-plate B-splines (Bookstein, 1989).

### References

Bernadette S. de Bakker, Kees H. de Jong, Jaco Hagoort, Karel de Bree, Clara T. Besselink, Froukje E. C. de Kanter, Tyas Veldhuis, Babette Bais, Reggie Schildmeijer, Jan M. Ruijter, Roelof-Jan Oostra, Vincent M. Christoffels, Antoon F. M. Moorman (2016) An interactive three-dimensional digital atlas and quantitative database of human development *Science* , Vol. 354, Issue 6315

Bookstein (1989) Principal warps: thin-plate splines and the decomposition of deformations, *IEEE Transactions on Pattern Analysis and Machine Intelligence*, vol. 11, no. 6, pp. 567-585

Daniel Fürth, Thomas Vaissière, Ourania Tzortzi, Yang Xuan, Antje Märtin, Iakovos Lazaridis, Giada Spigolon, Gilberto Fisone, Raju Tomer, Karl Deisseroth, Marie Carlén, Courtney A. Miller, Gavin Rumbaugh & Konstantinos Meletis (2018) An interactive framework for whole-brain maps at cellular resolution, *Nature Neuroscience* vol. 21, pp. 139–149
