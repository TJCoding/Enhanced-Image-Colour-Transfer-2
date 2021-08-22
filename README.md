# Enhanced Image-Colour-Transfer-2
============================================================================

##### A Further Enhanced Implementation of the Colour Transfer Method proposed by E Reinhard et al.

UNDER CONSTRUCTION

A further revised version of the standard basic processing method is proposed which utilises image representation in the l-alpha-beta  colour space and which supersedes the Enhanced Method described [here](https://github.com/TJCoding/Enhanced-Image-Colour-Transfer). As with the basic processing, account is taken of the the mean and standard deviation values of the respective l-alpha-beta components but, in addition, the correlation between the 'alpha' and 'beta' colour components is adddressed. Additional new options allow for the matching of higher moments of the colour components beyond the second moments, for the adjustment of image saturation and for the implementation of flexible image shading options.  The processing is described further in the file '[Enhanced Image-Colour-Transfer.pdf](Documents/Further%20Enhanced%20Image-Colour-Transfer.pdf)'  in the 'Documents' sub-folder. 

The program '[Main.cpp ](Main.cpp)' runs under C++ using OpenCV. Processing options may be specified by modifying statements within the 'Processing Selections' section in the main routine.

The following images illustrate the processing describe here.

*The first image in each set is the target image and the second is the colour source image.  The third image shows the effects of standard basic processing as proposed by E Reinhard et al. The final image has been subject to the further enhanced processing using the default processing options specified in the source code and in its asssociated documentation.* 

![Composite of Vase Image: Inputs and Outputs](Documents/Images/Vase_composite.jpg?raw=true)

![Composite of Autumn Image: Inputs and Outputs](Documents/Images/Autumn_composite.jpg?raw=true)

![Composite of Flowers Image: Inputs and Outputs](Documents/Images/Flowers_composite.jpg?raw=true)

An executable for the 'Further Enhanced-Image-Colour-Transfer' processing  is available [here](https://github.com/TJCoding/Image-Colour-Transfer-Processing-Executable).



