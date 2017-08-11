# Cell-Counter
The following program detects and quantifies the fluorescence intensity in HeLa cells after immunocytochemistry staining of thymine dimers staining with Alexa 488. 

Reference: https://doi.org/10.1016/j.bpj.2011.09.031

Detailed documentation: 
Images for analysis were collected by epifluorescence imaging with a CCD camera. The program selects for bright regions above the background (damaged cell nuclei). The program then quantifies the fluorescence intensity in each cell, excludes any bright regions (due to non-specific fluorophore clumping), averages the total fluorescent intensity of all the cells, and counts the number of cells present.  
The script contains several adjustable parameters. The name of the image is set by the name string. The images contain a region of damaged cells that are selected by adjusting the variables x and y which correspond to the coordinates of the upper left corner of the bounding box shown in MatLab output Figure 1 (Figure E.5A). The length of the bounding box is set by the xsize and ysize parameters. Also in Figure 1 is a small box that is used the set the background fluorescence. Adjusting the variables av1, av2, and thresh determine the thresholding requirements to determine damaged cell nuclei.  The mask is shown in MatLab output Figure 2 (Figure E.5B). The output variable meanpix1_bkgcorr gives the average intensity of the cells in Figure 2 after subtracting the background intensity. Output Figure 4 (Figure E.5D) is the same mask as Figure 2 except any cell with an average intensity greater than cutoff has been excluded from the analysis. The output variable meanpix2_bkgcorr gives the average intensity of the cells in Figure 4 after subtracting the background intensity.  Finally, output Figure 5 (Figure E.6) counts the number of cells in the region of interest. The cell count is given by the output variable cells. The minimum area for a cell is set by cellsize and the number of cells meeting this criterion are output in cells_threshold.  
