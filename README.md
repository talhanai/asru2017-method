# asru2017-method

This repo contains:
1. [**Step1_processAudio.m**] feature extraction pipeline.<br>
   Input: 1234.lld.csv, 1234.rttm.txt<br>
   Output: DATA.mat, audioFeatures.csv<br>
2. [**Step2_modeling.m**] method applied to the features.<br>
   Input: audioFeatures.csv, outcomes.csv<br>
   Output: cal variable in workspace.<br>

## Reference
This work was developed to detect individuals with cognitive impairment via features found in speech and language, and to determine the predictive power of such features. A core constraint that lead to using the methods described was the low ratio of examples to features (less than 1 to 10). You can find the full details here:

   [Spoken Language Biomarkers for Detecting Cognitive Impairment](https://groups.csail.mit.edu/sls/publications/2017/ASRU17_alhanai.pdf) <br>
   T. Alhanai, R. Au, and J. Glass, IEEE ASRU, Dec 2017, Japan

```
@inproceedings{alhanai2017spoken,
  title={Spoken Language Biomarkers for Detecting Cognitive Impairment},
  author={Alhanai, Tuka and Au, Rhoda and Glass, James},
  booktitle={Automatic Speech Recognition \& Understanding, 2017. ASRU. IEEE Workshop on},
  year={2017},
  organization={IEEE}
}
```


## Feature extraction pipeline [**Step1_processAudio.m**]
The acoustic feature extraction pipeline is as follows, with the steps also denoted in the code itself:

0. It is assumed that you have timestamps for each speaker turn. This could be manually generated or through some automatic diarization pipeline. This example uses a file in the RTTM format (1234.rttm.txt).
   
   An RTTM file has the following format:
   
   ```
   SPEAKER 1234 1 0000.00 001.00 <NA> <NA> 1234-T <NA>
   SPEAKER 1234 1 0001.00 001.00 <NA> <NA> 1234-P <NA>
   SPEAKER 1234 1 0002.00 001.50 <NA> <NA> 1234-P <NA>
   ```
   which contains a file ID, channel ID, start time, duration, and speaker ID (e.g. 1234-T). 

1. Extract features (e.g. using openSMILE toolkit). These are frame-level (every 10ms) using the IS13_ComParE.conf configuration parameters ([pipeline here](https://github.com/talhanai/acousticfeatures-fhs)). The (dummy) features used in this example are 1234.lld.csv.

   The features are in openSMILE notation" mfcc_sma[1-14],F0final_sma, voicingFinalUnclipped_sma, jitterLocalDDP_sma, jitterDPP_sma, shimmerLocal_sma, logHNR_sma, pcm_RMSenergy_sma, pcm_zcr_sma as well as \_de (differential) of these features. 
   
   (if you see/have audspec\*, we will be ignoring these)

2. Zero Mean Variance Normalize (i.e. zscore) the energy-based features over the complete utterance. Specifically mfcc_sma\*, pcm_zcr\*, pcm_RMSenergy\*, logHNR*.

3. Smooth voicingFinalUnclipped_sma using an envelope of 10 points (100 ms) using a spline over the local maxima.

4. For each speaker turn (a subject), find the frames (indices) that are above a voicing threshold (0.1 standard deviations from the mean).

5. (Optional) Using these speech frames, calculate segment-level statistics for this speaker turn (mean, maximum, minimum, median, standard deviation) over all the features, including the original voicing\* (i.e. not Step 3. smoothed), and except for F0\*, shimmer\*, and jitter\*.

   For F0\* calculate the global statistics for their non-zero values within the speech segments. Do the same for shimmer\* and jitter\*.
   
   (Toss out zero valued pitch - FO\*, shimmer\*, jitter\* - frames)

6. Now calculate the the subject-level statistics (mean, maximum, minimum, median, std).

You may choose to work with the features from Step 1, Step 5, or Step 6, and vary the processing performed in the intermediary steps.

## Method [**Step2_modeling.m**]

1. An initial set of features were selected if they had a statistically significant (p < 0.01) univariate Pearson correlation with the (training set, i.e. N − 1 subjects) outcome. 

2.  An Elastic-net regularized binomial logistic regression is trained with these features, resulting in further feature selection. Features with non-zero model coefficients (β) were selected. 

   [This Elastic-net implementation](https://web.stanford.edu/~hastie/glmnet_matlab/) was used.  

3. This model was evaluated on the held out leave-one-out test subject. 

4. This training method was repeated for all leave-one-out (N) folds.

![alt tag](https://github.com/talhanai/asru2017-method/blob/master/DiagramFeature.png?raw=true)
