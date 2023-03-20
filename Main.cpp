//*** IMPLEMENTATION OF A FURTHER ENHANCED ADAPTATION
//*** OF THE REINHARD COLOUR TRANSFER METHOD.
//    See 'main' routine for further details.
//
// Copyright © Terry Johnson, January 2020
// Revised 27/05/2021 (Version 5)
// https://github.com/TJCoding

#include <opencv2/highgui/highgui.hpp>
#include <opencv2/photo/photo.hpp>

// Declare functions
cv::Mat CoreProcessing(cv::Mat targetf, cv::Mat sourcef,
                       float CrossCovarianceLimit,
                       int   ReshapingIterations,
                       float ShaderVal);
cv::Mat adjust_covariance(cv::Mat Lab[3], cv::Mat sLab[3],
                          float covLim);
cv::Mat ChannelCondition(cv::Mat tChan, cv::Mat sChan);
cv::Mat SaturationProcessing(cv::Mat targetf, cv::Mat savedtf,
                             float SatVal);
cv::Mat FullShading(cv::Mat targetf, cv::Mat savedtf, cv::Mat sourcef,
                    bool ExtraShading, float ShadeVal);
cv::Mat FinalAdjustment(cv::Mat targetf, cv::Mat savedtf,
                        float TintVal, float ModifiedVal);
cv::Mat convertTolab  (cv::Mat input);
cv::Mat convertFromlab(cv::Mat input);



int main(int argc, char *argv[])
{
//  Transfers the colour distribution from the source image to the
//  target image by matching the mean and standard deviation and
//  the colour cross correlation in the L-alpha-beta colour space.
//  Additional image refinement options are also provided.

//	The implementation is an enhancement of a method described in
//  "Color Transfer between Images" paper by Reinhard et al., 2001,
//  but with additional processing options as follows.

// PROCESSING OPTIONS

//  OPTION 1
//  There is an option to match the cross correlation between
//  the colour channels 'alpha' and 'beta'.  Full matching, no
//  matching or restricted matching may be specified.
//  (See the note at the end of the code).

//  OPTION 2
//  There is an option to further match the distribution of the
//  target 'alpha' and 'beta' colour channels to the corresponding
//  source channels by iteratively adjusting the higher order
//  characteristics of the channel data (skewness and kurtosis)
//  in addition to matching the mean and standard deviation values.
//  The number of iterative adjustments can be specified for this
//  'reshaping' option. A zero value specifies no reshaping processing
//  and a non-zero value specifies the total number of reshaping
//  iterations to be implemented.

//  OPTION 3
//  There is an option to pare back excess colour saturation
//  resulting from the colour transfer processing.

//  OPTION 4
//  There is an option to retain the shading of the target image
//  so that the processing is implemented as a pure colour transfer
//  or to adopt the shading of the source image or to select an
//  intermediate shading.

//  OPTION 5
//  There is an option to further adjust the shading of the output
//  image to match the target image, the source image or an
//  intermediate image (as selected under Option 4) by directly
//  matching distributions in the corresponding grey shade images.

//  OPTION 6
//  There is an option to modify the colour tint of the final image.

//  OPTION 7
//  There is an option to mix the final image with the initial image,
//  so that only part modification occurs.

// ##########################################################################
// #######################  PROCESSING SELECTIONS  ##########################
// ##########################################################################
    // Select the processing options in accordance with the
    // preceding descriptions.

    float CrossCovarianceLimit     = 0.5;    // Option 1 (Default is '0.5')
    int   ReshapingIterations      = 1;      // Option 2 (Default is '1')
    float PercentSaturationShift   = -1.0;   // Option 3 (Default is -1.0)
    float PercentShadingShift      = 50.0;   // Option 4 (Default is 50.0)
    bool  ExtraShading             = true;   // Option 5 (Default is 'true')
    float PercentTint              = 100.0;  // Option 6 (Default is 100.0)
    float PercentModified          = 100.0;  // Option 7 (Default is 100.0)

   //  Setting CrossCovarianceLimit to 0.0 inhibits cross covariance processing.
   //  Setting ReshapingIterations to 0, inhibits reshaping processing.
   //  Setting PercentShadingShift to 0, retains the target image saturation.
   //  Setting PercentShadingShift to 0, retains the target image shading.
   //  Setting ExtraShading to 'false' reverts to simple shading.
   //  Setting PercentTint to '0', gives a monochrome image.
   //  Setting PercentModified to '0', retains the target image in full.

   //  For each of the percentage parameters, defined above, a setting of '100'
   //  allows the full processing effect.  A setting of '0' suppresses the
   //  particular processing effect. Intermediate values give an intermediate outcome
   //  Percentages below '0' or above '100' are permitted although the outcome
   //  may not always be pleasing. (See note below for saturation shift.)
   //
   // Note:
   // If a negative percentage is set for saturation shift, then the actual
   // percentage is determined automatically from the image properties.


   // Specify the image files that are to be processed,
   // where 'source image' provides the colour scheme
   // that is to be applied to 'target image'.

    std::string targetname = "images/Flowers_target.jpg";
    std::string sourcename = "images/Flowers_source.jpg";

// ###########################################################################
// ###########################################################################
// ###########################################################################

    // Read in the images and convert to floating point,
    // saving a copy of the target image for later.
    cv::Mat target = cv::imread(targetname, 1);
    cv::Mat source = cv::imread(sourcename, 1);

    cv::Mat targetf(target.size(),CV_32FC3);
    cv::Mat sourcef(source.size(),CV_32FC3);
    target.convertTo(targetf, CV_32FC3, 1.0/255.f);
    source.convertTo(sourcef, CV_32FC3, 1.0/255.f);
    cv::Mat savedtf=targetf.clone();

    // Implement augmented "Reinhard Processing" in
    // L-alpha-beta colour space.
    targetf=CoreProcessing(targetf, sourcef, CrossCovarianceLimit,
                           ReshapingIterations,
                           PercentShadingShift/100.0);
    cv::cvtColor(sourcef,sourcef,CV_BGR2GRAY); // Only need mono hereafter.

    // Implement image refinements where a change is specified.
    targetf=SaturationProcessing(targetf, savedtf,
                         PercentSaturationShift/100.0);
    targetf=FullShading(targetf, savedtf, sourcef, ExtraShading,
                        PercentShadingShift/100.0);
    targetf=FinalAdjustment(targetf,savedtf,
                            PercentTint/100.0,
                            PercentModified/100.0);

    //  Convert the processed image to integer format.
    cv::Mat result;

    targetf.convertTo(result, CV_8UC3, 255.f);
    // Display and save the final image.
    cv::imshow("processed image",result);
    cv::imwrite("images/processed.jpg", result);

    // Display image until a key is pressed.
    cv::waitKey(0);
    return 0;
   }



cv::Mat CoreProcessing(cv::Mat targetf, cv::Mat sourcef,
                       float CrossCovarianceLimit,
                       int   ReshapingIterations,
                       float ShaderVal)
{
// Implements augmented "Reinhard Processing" in
// L-alpha-beta colour space.

    // First convert the images from the BGR colour
    // space to the L-alpha-beta colour space.
    // Estimate the mean and standard deviation of
    // colour channels. Split the target and source
    // images into colour channels and standardise
    // the distribution within each channel.
    //
    // The standardised data has zero mean and
    // unit standard deviation.
    cv::Mat Lab[3], sLab[3];
    cv::Scalar tmean, tdev, smean, sdev;

    targetf = convertTolab(targetf);
    sourcef = convertTolab(sourcef);

    cv::meanStdDev(targetf, tmean, tdev);
    cv::meanStdDev(sourcef, smean, sdev);
    cv::split(targetf,Lab);
    cv::split(sourcef,sLab);

    Lab[0]=(Lab[0]-tmean[0])/tdev[0];
    Lab[1]=(Lab[1]-tmean[1])/tdev[1];
    Lab[2]=(Lab[2]-tmean[2])/tdev[2];

    sLab[0]=(sLab[0]-smean[0])/sdev[0];
    sLab[1]=(sLab[1]-smean[1])/sdev[1];
    sLab[2]=(sLab[2]-smean[2])/sdev[2];


    // Implement first phase of reshaping for the colour channels
    // when one or more iteration is specified.
    int jcount=ReshapingIterations;
    while (jcount>ceil((ReshapingIterations+1)/2))
    {
         Lab[1]=ChannelCondition(Lab[1],sLab[1]);
         Lab[2]=ChannelCondition(Lab[2],sLab[2]);
         jcount--;
     }
    // Implement cross covariance processing.
    // (null if CrossCovarianceLimit=0.0)
        targetf=adjust_covariance(Lab, sLab,
                       CrossCovarianceLimit );
        cv::split(targetf,Lab);

    // Implement second phase of reshaping
    while (jcount>0)
    {
         Lab[1]=ChannelCondition(Lab[1],sLab[1]);
         Lab[2]=ChannelCondition(Lab[2],sLab[2]);
         jcount--;
     }

    // Rescale the previously standardised colour channels
    // so that the means and standard deviations now match
    // those of the source image.
    Lab[1]=Lab[1]*sdev[1]+smean[1];
    Lab[2]=Lab[2]*sdev[2]+smean[2];

    // Rescale the lightness channel (channel 0)
    // in accordance with the specified percentage
    // shading shift.
    Lab[0]=Lab[0]*(ShaderVal*sdev[0]+(1.0-ShaderVal)*tdev[0])
           +ShaderVal*smean[0]+(1.0-ShaderVal)*tmean[0];

    // Merge channels and convert back to BGR colour space.
    cv::Mat resultant;
    cv::merge(Lab,3,resultant);
    resultant = convertFromlab(resultant);
    return resultant;
}



cv::Mat adjust_covariance(cv::Mat Lab[3], cv::Mat sLab[3], float covLim)
{
// This routine adjusts colour channels 2 and 3 of
// the image within the L-alpha-beta colour space.

// The channels each have zero mean and unit
// standard deviation but their cross correlation
// value will not normally be zero.
//
// The processing reduces the cross correlation between
// the channels to zero but then reintroduces
// correlation such that the new cross correlation
// value matches that for the source image.
//
// Throughout these manipulations the mean channel values
// are maintained at zero and the standard deviations are
// maintained as unity.
//
// The manipulations are based upon the following relationship.
//
// Let z1 and z2 be two independent (zero correlation) variables
// with zero means and unit standard deviations. It can be shown
// that variables a1 and a2 have zero means, unit standard
// deviations, and mutual cross correlation 'R' when:
//
// a1=sqrt((1+R)/2)*z1 + sqrt((1-R)/2)*z2
// a2=sqrt((1+R)/2)*z1 - sqrt((1-R)/2)*z2
//
// The above relationships are applied inversely to derive
// uncorrelated standardised colour channels variables from
// the standardised but correlated input channels.
//
// The above relationships are then applied directly to obtain
// standardised correlated colour channels with correlation
// matched to that of the source image colour channels.
//
// Original processing method attributable to Dr T E Johnson Sept 2019.

    // Declare variables
    float tcrosscorr, scrosscorr;
    float W1, W2, norm;
    cv::Mat temp1;
    cv::Scalar smean, sdev, temp2;

    // No processing required if 'covLim' set to zero.
    if(covLim!=0.0)
    {
        // Compute the correlation for the target image colour channels.
        // The correlation between the standardised variables (zero mean,
        // unit standard deviation) can be computed as the mean cross
        // product for the two channels.
        cv::multiply(Lab[1],Lab[2],temp1);
        temp2=cv::mean(temp1);
        tcrosscorr =temp2[0];

        // Compute the correlation for the source
        // image colour channels.
        cv::multiply(sLab[1],sLab[2],temp1);
        temp2=cv::mean(temp1);
        scrosscorr =temp2[0];

        // Adjust the correlation between the
        // standardised input channel values.
        W1= 0.5*sqrt((1+scrosscorr)/(1+tcrosscorr))
           +0.5*sqrt((1-scrosscorr)/(1-tcrosscorr));
        W2= 0.5*sqrt((1+scrosscorr)/(1+tcrosscorr))
           -0.5*sqrt((1-scrosscorr)/(1-tcrosscorr));

        // Limit the size of W2 if required.
        // This limits the proportional amount by which
        // a given colour channel can be augmented by
        // the energy from the other colour channel.
        if(std::abs(W2)>covLim*std::abs(W1))
        {
            W2=copysign(covLim*W1,W2);
            norm=1.0/sqrt(W1*W1+W2*W2+2*W1*W2*tcrosscorr);
            W1=W1*norm;
            W2=W2*norm;
        }
        cv::Mat z1=Lab[1].clone();

        Lab[1]=W1*z1+W2*Lab[2];
        Lab[2]=W1*Lab[2]+W2*z1;
    }
    cv:merge(Lab,3,temp1);
    return temp1;
}



cv::Mat ChannelCondition(cv::Mat Chan,cv::Mat sChan)
    {
// Modifies the distribution of values in 'Chan' to more
// closely match the distribution of those in 'sChan'.
// Separate matching operations are performed for values
// above and below the mean.  The input channels have
// been standardised so the mean is equal to zero.
// Original processing method attributable to
// Dr T E Johnson Oct 2020.

    // Declare variables
    // Computations use weighted data values.
    // 'wval' is the tuning constant for the
    // weighting function.
    cv::Mat mask, ChanU, ChanL;
    cv::Mat WU, WL;
    cv::Scalar smeanU, smeanL, wmean;
    cv::Scalar tmeanU, tmeanL, tmean, tdev;
    float k, wval=0.25;

    // sChan is processed before Chan because
    // Chan data is used in later processing.

    // Processing for upper 'sChan'.

    // Determine the mask for selecting data values
    // above zero.
    cv::threshold(sChan,mask,0,1,CV_THRESH_BINARY);
    mask.convertTo(mask,CV_8U);

    // Compute the weighting function for values
    // above zero.
    // (Zero is the mean value of the input channel).
    // The weighting function is zero for sChan values
    // equal to zero and unity for large values of
    // sChan.
    smeanU=mean(sChan,mask);
    cv::exp(-sChan*wval/smeanU[0],WU);
    WU=(1-WU).mul(1-WU);
    wmean=mean(WU,mask);
    // Compute deviation from the mean
    // and raise to the power 4 so as to
    // address kurtosis.
    cv::pow(sChan,4,ChanU);
    // Find the weighted average of the
    // fourth power of the deviations.
    smeanU=mean(WU.mul(ChanU),mask)/wmean[0];

    // Processing for lower 'sChan'.

    // As for upper processing but values are
    // selected by applying the complementary
    // masking function (1-mask).
    smeanL=mean(sChan,(1-mask));
    cv::exp(-sChan*wval/smeanL[0],WL);
    WL=(1-WL).mul(1-WL);
    wmean=mean(WL,1-mask);
    cv::pow(sChan,4,ChanL);
    smeanL=mean(WL.mul(ChanL),1-mask)/wmean[0];

    // Processing for upper 'Chan'

    cv::threshold(Chan,mask,0,1,CV_THRESH_BINARY);
    mask.convertTo(mask,CV_8U);
    tmeanU=mean(Chan,mask);
    cv::exp(-Chan*wval/tmeanU[0],WU);
    WU=(1-WU).mul(1-WU);
    wmean=mean(WU,mask);
    cv::pow(Chan,4,ChanU);
    tmeanU=mean(WU.mul(ChanU),mask)/wmean[0];

    // Processing for lower 'Chan'

    tmeanL=mean(Chan,(1-mask));
    cv::exp(-Chan*wval/tmeanL[0],WL);
    WL=(1-WL).mul(1-WL);
    wmean=mean(WL,1-mask);
    cv::pow(Chan,4,ChanL);
    tmeanL=mean(WL.mul(ChanL),1-mask)/wmean[0];

    // Modify the upper 'Chan' values

    // Compute the ratio of the weighted fourth
    // power for 'sChan' relative to that for
    // 'Chan' and then take the fourth root.
    // The resultant is used to apply a shift to
    // the 'Chan' data where the shift is a
    // function of the data deviation.
    // No shift is applied to small values and full
    // shift to large values.
    k=sqrt(sqrt(smeanU[0]/tmeanU[0]));
    ChanU=(1+WU*(k-1)).mul(Chan);

    // Similarly modify the lower 'Chan' values.
    k=sqrt(sqrt(smeanL[0]/tmeanL[0]));
    ChanL=(1+WL*(k-1)).mul(Chan);

    // Combine the upper and lower 'Chan'values to form
    // a whole
    Chan=cv::Mat::zeros(Chan.rows, Chan.cols, CV_32FC1);
    cv::add(Chan,ChanU,Chan,mask);
    cv::add(Chan,ChanL,Chan,1-mask);

    // Re-standardise the modified 'Chan' data
    // before it is fed back.
    cv::meanStdDev(Chan, tmean, tdev);
    Chan=(Chan-tmean[0])/tdev[0];
    cv::meanStdDev(Chan, tmean, tdev);

    return Chan;
    }



cv::Mat SaturationProcessing(cv::Mat targetf, cv::Mat savedtf,
                             float SatVal)
{
// This routine allows a reduction of colour saturation
// to an extent specified by the parameter 'SatVal'.
// An image that is subject to full colour transfer can
// often exhibit excessive colour saturation. This function
// allows a scaling back of saturation.

// The idea of saturation processing is to adjust the saturation
// characteristics of the processed image to match the saturation
// characteristics of an artificially constructed image whose
// saturation characteristics are considered desirable.

    // Implement a saturation change unless 100% saturation
    // is specified.
    if (SatVal!=1)
    {
        cv::Mat temp, mask, Hsv[3], tmpHsv[3];
        cv::Scalar tmean, tdev, tmpmean, tmpdev;

        // Colour saturation will be computed in accordance
        // with the definition used for the HSV colour space.
        cv::cvtColor(targetf,targetf,CV_BGR2HSV);
        cv::cvtColor(savedtf,temp,CV_BGR2HSV);
        cv::split(targetf,Hsv);
        cv::split(temp,tmpHsv);

        if(SatVal<0)
        {
        //  'SatVal' is less than 0, then compute 'SatVal'
        //  as the ratio of the largest saturation value
        //  in the original image to the largest value in
        //  processed image.
            double amin;
            double amax1,amax2;
            cv::minMaxIdx(Hsv[1],&amin,&amax1);
            cv::minMaxIdx(tmpHsv[1],&amin,&amax2);
            SatVal=amax2/amax1;
        }

        // Compute a weighted mix of the processed target
        // saturation channel  and the original image
        // saturation channel to define an initial
        // reference saturation channel.
        cv::addWeighted(Hsv[1],SatVal,tmpHsv[1],
                        1-SatVal,0.0,tmpHsv[1]);

        // The initial reference saturation channel values
        // will apply only to those pixels where the
        // saturation in the processed image exceeds that
        // in the original target image.
        cv::threshold((Hsv[1]-tmpHsv[1]),mask,0,1,CV_THRESH_BINARY);
        mask.convertTo(mask,CV_8U);

        // Create a new reference saturation channel which is taken
        // from the initial reference saturation channel except
        // where the mask function indicates that the saturation
        // value of the original image should be used.
        // This gives a modified reference saturation channel.
        temp=cv::Mat::zeros(targetf.rows, targetf.cols, CV_32FC1);
        cv::add(temp,tmpHsv[1],temp,mask);
        cv::add(temp,Hsv[1],tmpHsv[1],(1-mask));

        // Now match the mean and standard deviation of the
        // saturation channel for the processed image channel
        // to the mean and standard deviation of the modified
        // reference saturation channel. This give the final
        // saturation channel which is the output from
        // the saturation processing
        cv::meanStdDev(Hsv[1], tmean, tdev);
        cv::meanStdDev(tmpHsv[1], tmpmean, tmpdev);
        Hsv[1]=(Hsv[1]-tmean[0])/tdev[0];
        Hsv[1]=Hsv[1]*tmpdev[0]+tmpmean[0];
        cv::merge(Hsv,3,targetf);
        cv::cvtColor(targetf,targetf,CV_HSV2BGR);
    }
    return targetf;
}



cv::Mat FullShading(cv::Mat targetf, cv::Mat savedtf, cv::Mat sourcef,
                    bool ExtraShading, float ShaderVal)
    {
     // Matches the grey shade distribution of the
     // modified target image to that of a notional
     // shader image which is a linear combination
     // of the original target and source image as
     // determined by the value of 'ShaderVal'.

     if(ExtraShading)
     {
         cv::Mat greyt, greys, greyp, chans[3];
         cv::Scalar smean, tmean, sdev, tdev;

         // Compute the grey shade images for the target,
         // processed and source images.
         cv::cvtColor(savedtf,greyt,CV_BGR2GRAY);
         cv::cvtColor(targetf,greyp,CV_BGR2GRAY);
         sourcef.copyTo(greys);// Already converted.

         // Standardise the greyshade images
         // for the source and target.
         cv::meanStdDev(greys, smean, sdev);
         cv::meanStdDev(greyt, tmean, tdev);
         greyt=(greyt-tmean[0])/tdev[0];

         // Rescale the previously standardised grey shade
         // target image so that the means and standard
         // deviations now match those of the notional shader image.
         greyt=greyt*(ShaderVal*sdev[0]+(1.0-ShaderVal)*tdev[0])
               +ShaderVal*smean[0]+(1.0-ShaderVal)*tmean[0];

         // Rescale each of the colour channels of the
         // processed image identically so that in grey
         // shade the processed image more closely matches
         // the grey shading of the nominated goal image.
         //
         cv::Mat min_mat = cv::Mat(greyp.size(), CV_32FC1,1/255.0);
	     cv::max(greyp, min_mat, greyp); // Guard against zero divide;
	     // Guard against negative values;
	     cv::max(greyt, cv::Mat(greyt.size(), CV_32FC1, 0.0), greyt);



         cv::split(targetf,chans);
         cv::divide(chans[0],greyp,chans[0]);
         chans[0]=chans[0].mul(greyt);
         cv::divide(chans[1],greyp,chans[1]);
         chans[1]=chans[1].mul(greyt);
         cv::divide(chans[2],greyp,chans[2]);
         chans[2]=chans[2].mul(greyt);
         cv::merge(chans,3,targetf);
     }

     return targetf;
    }



cv::Mat FinalAdjustment(cv::Mat targetf,cv::Mat savedtf,
                        float TintVal, float ModifiedVal)
{
// Implements a change to the tint of the final image and
// to its degree of modification if a change is specified.

    // If 100% tint not specified then compute a weighted average
    // of the processed image and its grey scale representation.
    if(TintVal!=1.0)
     {
         cv::Mat grey;
         cv::Mat BGR[3];
         cv::cvtColor(targetf,grey,CV_BGR2GRAY);
         cv::split(targetf,BGR);
         BGR[0]=TintVal*BGR[0]+(1.0-TintVal)*grey;
         BGR[1]=TintVal*BGR[1]+(1.0-TintVal)*grey;
         BGR[2]=TintVal*BGR[2]+(1.0-TintVal)*grey;
         cv::merge(BGR,3,targetf);
     }

    // If 100% image modification not specified then
    // compute a weighted average of the processed image
    // and the original target image.
     if(ModifiedVal!=1.0)
     {
        targetf=ModifiedVal*targetf+(1.0-ModifiedVal)*savedtf;
     }

   return targetf;
}



// ##########################################################################
// ##### IMPLEMENTATION OF L-ALPHA-BETA FORWARD AND INVERSE TRANSFORMS ######
// ##########################################################################
// Coding taken from https://github.com/ZZPot/Color-transfer
// Credit to 'ZZPot'.
// I take responsibility for any issues arising my adaptation.


// Define the transformation matrices for L-alpha-beta transformation.
cv::Mat RGB_to_LMS = (cv::Mat_<float>(3,3) <<	0.3811f, 0.5783f, 0.0402f,
										        0.1967f, 0.7244f, 0.0782f,
                                                0.0241f, 0.1288f, 0.8444f);
float i3 = 1/sqrt(3), i6 = 1/sqrt(6), i2 = 1/sqrt(2);
cv::Mat LMS_to_lab = (cv::Mat_<float>(3,3) <<	i3, i3, i3,
										        i6, i6, -2*i6,
                                                i2, -i2, 0);



cv::Mat convertTolab(cv::Mat input)
{
	cv::Mat img_RGBf(input.size(),CV_32FC3);
    cv::Mat img_lms (input.size(),CV_32FC3);
    cv::Mat img_lab (input.size(),CV_32FC3);

	// Swap channel order (so that transformation
	// matrices can be used in their familiar form).

	// Then convert to float.
	//cv::cvtColor(input, img_RGB, CV_BGR2RGB);
	//img_RGB.convertTo(img_RGBf, CV_32FC3, 1.0/255.f);

    cv::cvtColor(input, img_RGBf, CV_BGR2RGB);


	// Apply stage 1 transform.
	cv::transform(img_RGBf, img_lms, RGB_to_LMS);

	// Define smallest permitted value and implement it.
    float epsilon =1.0/255;
	cv::Scalar min_scalar(epsilon, epsilon, epsilon);
	cv::Mat min_mat = cv::Mat(input.size(), CV_32FC3, min_scalar);
	cv::max(img_lms, min_mat, img_lms); // just before log operation.

	// Compute log10(x)  as ln(x)/ln(10)
	cv::log(img_lms,img_lms);
	img_lms=img_lms/log(10);

    // Apply stage 2 transform.
	cv::transform(img_lms, img_lab, LMS_to_lab);

	return img_lab;
}



cv::Mat convertFromlab(cv::Mat input)
{
	cv::Mat img_lms (input.size(),  CV_32FC3);
	cv::Mat img_RGBf(input.size(),  CV_32FC3);
	cv::Mat img_BGRf(input.size(),  CV_32FC3);
    cv::Mat temp    (LMS_to_lab.size(),CV_32FC1);

    // Apply inverse of stage 2 transformation.
	cv::invert(LMS_to_lab,temp);
	cv::transform(input, img_lms, temp);

	// Compute 10^x as (e^x)^(ln10)
    cv::exp(img_lms,img_lms);
    cv::pow(img_lms,(double)log(10.0),img_lms);

    // Apply inverse of stage 1 transformation.
	cv::invert(RGB_to_LMS,temp);
	cv::transform(img_lms, img_RGBf, temp);

	//  Revert channel ordering to BGR.
	cv::cvtColor(img_RGBf, img_BGRf, CV_RGB2BGR);

	return img_BGRf;
}



// ##########################################################################
// ##########################################################################
// ##########################################################################


// Notes on Cross Correlation Matching.
// ====================================
// Cross correlation matching is performed by operations of the
// form.
// Channel_alpha = W1*Channel_alpha + W2*Channel_beta
// Channel_beta  = W1*Channel_beta  + W2*Channel_alpha
// as determined by the value of CrossCovarianceLimit.
//
// If CrossCovarianceLimit = 0, W2=0 and no cross correlation
// matching is performed.
// If CrossCovarianceLimit > 0, W2 is clipped if necessary so
// that it cannot lie outside the range
// -CrossCovarianceLimit*W1 to +CrossCovarianceLimit*W1.
// This mechanism may be used to guard against an overly large
// correction term.
//
// Typically CrossCovarianceLimit might be set to 0.5
//(for a maximum modification corresponding to 50%).







