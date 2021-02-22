#!/usr/bin/env bash

DIR=`pwd`
mkdir $DIR/tmp
TMPDIR=$DIR/tmp
TMP=$TMPDIR
TEMP=$TMPDIR
export TMPDIR TMP TEMP

make -f repeat-clustering/makefile \
    CPUS=96 DATABASE=cottonRB23.04.fa \
    NAME_LENGTH=6 \
    LOCAL='A1_019      A1_084  A1_Nis  A2_034  A2_060  A2_079  A2_118  A2_164  A2_G35  A2_G46  A2_G57  A2_N03 A1_088  A2_001  A2_040  A2_061  A2_080  A2_119  A2_255  A2_G36  A2_G47  A2_G58  F1_lon A1_097  A2_004  A2_041  A2_062  A2_084  A2_123  A2_G01  A2_G37  A2_G48  A2_G59 A1_028      A1_108  A2_008  A2_043  A2_064  A2_085  A2_124  A2_G02  A2_G38  A2_G49  A2_G60 A1_030      A1_125  A2_016  A2_044  A2_065  A2_087  A2_131  A2_G03  A2_G39  A2_G50  A2_G61 A1_051      A1_132  A2_019  A2_045  A2_066  A2_091  A2_141  A2_G04  A2_G40  A2_G51  A2_G62 A1_054      A1_133  A2_021  A2_046  A2_067  A2_096  A2_142  A2_G30  A2_G41  A2_G52  A2_G63 A1_073      A1_155  A2_026  A2_047  A2_069  A2_100  A2_147  A2_G31  A2_G42  A2_G53  A2_G64 A1_074      A1_Afr  A2_029  A2_056  A2_073  A2_101  A2_154  A2_G32  A2_G43  A2_G54  A2_G65 A1_077      A1_G01  A2_030  A2_057  A2_074  A2_113  A2_156  A2_G33  A2_G44  A2_G55  A2_G66 A1_079      A1_G02  A2_031  A2_059  A2_076  A2_117  A2_158  A2_G34  A2_G45  A2_G56  A2_lin'  \
    A1_019=2-sample/A1_019.fa  A1_019_SIZE=185222 \
    A1_028=2-sample/A1_028.fa  A1_028_SIZE=185222 \
    A1_030=2-sample/A1_030.fa  A1_030_SIZE=185222 \
    A1_051=2-sample/A1_051.fa  A1_051_SIZE=185222 \
    A1_054=2-sample/A1_054.fa  A1_054_SIZE=185222 \
    A1_073=2-sample/A1_073.fa  A1_073_SIZE=185222 \
    A1_074=2-sample/A1_074.fa  A1_074_SIZE=185222 \
    A1_077=2-sample/A1_077.fa  A1_077_SIZE=185222 \
    A1_079=2-sample/A1_079.fa  A1_079_SIZE=185222 \
    A1_084=2-sample/A1_084.fa  A1_084_SIZE=185222 \
    A1_088=2-sample/A1_088.fa  A1_088_SIZE=185222 \
    A1_097=2-sample/A1_097.fa  A1_097_SIZE=185222 \
    A1_108=2-sample/A1_108.fa  A1_108_SIZE=185222 \
    A1_125=2-sample/A1_125.fa  A1_125_SIZE=185222 \
    A1_132=2-sample/A1_132.fa  A1_132_SIZE=185222 \
    A1_133=2-sample/A1_133.fa  A1_133_SIZE=185222 \
    A1_155=2-sample/A1_155.fa  A1_155_SIZE=185222 \
    A1_Afr=2-sample/A1_Afr.fa  A1_Afr_SIZE=185222 \
    A1_G01=2-sample/A1_G01.fa  A1_G01_SIZE=185222 \
    A1_G02=2-sample/A1_G02.fa  A1_G02_SIZE=185222 \
    A1_Nis=2-sample/A1_Nis.fa  A1_Nis_SIZE=185222 \
    A2_001=2-sample/A2_001.fa  A2_001_SIZE=190000 \
    A2_004=2-sample/A2_004.fa  A2_004_SIZE=190000 \
    A2_008=2-sample/A2_008.fa  A2_008_SIZE=190000 \
    A2_016=2-sample/A2_016.fa  A2_016_SIZE=190000 \
    A2_019=2-sample/A2_019.fa  A2_019_SIZE=190000 \
    A2_021=2-sample/A2_021.fa  A2_021_SIZE=190000 \
    A2_026=2-sample/A2_026.fa  A2_026_SIZE=190000 \
    A2_029=2-sample/A2_029.fa  A2_029_SIZE=190000 \
    A2_030=2-sample/A2_030.fa  A2_030_SIZE=190000 \
    A2_031=2-sample/A2_031.fa  A2_031_SIZE=190000 \
    A2_034=2-sample/A2_034.fa  A2_034_SIZE=190000 \
    A2_040=2-sample/A2_040.fa  A2_040_SIZE=190000 \
    A2_041=2-sample/A2_041.fa  A2_041_SIZE=190000 \
    A2_043=2-sample/A2_043.fa  A2_043_SIZE=190000 \
    A2_044=2-sample/A2_044.fa  A2_044_SIZE=190000 \
    A2_045=2-sample/A2_045.fa  A2_045_SIZE=190000 \
    A2_046=2-sample/A2_046.fa  A2_046_SIZE=190000 \
    A2_047=2-sample/A2_047.fa  A2_047_SIZE=190000 \
    A2_056=2-sample/A2_056.fa  A2_056_SIZE=190000 \
    A2_057=2-sample/A2_057.fa  A2_057_SIZE=190000 \
    A2_059=2-sample/A2_059.fa  A2_059_SIZE=190000 \
    A2_060=2-sample/A2_060.fa  A2_060_SIZE=190000 \
    A2_061=2-sample/A2_061.fa  A2_061_SIZE=190000 \
    A2_062=2-sample/A2_062.fa  A2_062_SIZE=190000 \
    A2_064=2-sample/A2_064.fa  A2_064_SIZE=190000 \
    A2_065=2-sample/A2_065.fa  A2_065_SIZE=190000 \
    A2_066=2-sample/A2_066.fa  A2_066_SIZE=190000 \
    A2_067=2-sample/A2_067.fa  A2_067_SIZE=190000 \
    A2_069=2-sample/A2_069.fa  A2_069_SIZE=190000 \
    A2_073=2-sample/A2_073.fa  A2_073_SIZE=190000 \
    A2_074=2-sample/A2_074.fa  A2_074_SIZE=190000 \
    A2_076=2-sample/A2_076.fa  A2_076_SIZE=190000 \
    A2_079=2-sample/A2_079.fa  A2_079_SIZE=190000 \
    A2_080=2-sample/A2_080.fa  A2_080_SIZE=190000 \
    A2_084=2-sample/A2_084.fa  A2_084_SIZE=190000 \
    A2_085=2-sample/A2_085.fa  A2_085_SIZE=190000 \
    A2_087=2-sample/A2_087.fa  A2_087_SIZE=190000 \
    A2_091=2-sample/A2_091.fa  A2_091_SIZE=190000 \
    A2_096=2-sample/A2_096.fa  A2_096_SIZE=190000 \
    A2_100=2-sample/A2_100.fa  A2_100_SIZE=190000 \
    A2_101=2-sample/A2_101.fa  A2_101_SIZE=190000 \
    A2_113=2-sample/A2_113.fa  A2_113_SIZE=190000 \
    A2_117=2-sample/A2_117.fa  A2_117_SIZE=190000 \
    A2_118=2-sample/A2_118.fa  A2_118_SIZE=190000 \
    A2_119=2-sample/A2_119.fa  A2_119_SIZE=190000 \
    A2_123=2-sample/A2_123.fa  A2_123_SIZE=190000 \
    A2_124=2-sample/A2_124.fa  A2_124_SIZE=190000 \
    A2_131=2-sample/A2_131.fa  A2_131_SIZE=190000 \
    A2_141=2-sample/A2_141.fa  A2_141_SIZE=190000 \
    A2_142=2-sample/A2_142.fa  A2_142_SIZE=190000 \
    A2_147=2-sample/A2_147.fa  A2_147_SIZE=190000 \
    A2_154=2-sample/A2_154.fa  A2_154_SIZE=190000 \
    A2_156=2-sample/A2_156.fa  A2_156_SIZE=190000 \
    A2_158=2-sample/A2_158.fa  A2_158_SIZE=190000 \
    A2_164=2-sample/A2_164.fa  A2_164_SIZE=190000 \
    A2_255=2-sample/A2_255.fa  A2_255_SIZE=190000 \
    A2_G01=2-sample/A2_G01.fa  A2_G01_SIZE=190000 \
    A2_G02=2-sample/A2_G02.fa  A2_G02_SIZE=190000 \
    A2_G03=2-sample/A2_G03.fa  A2_G03_SIZE=190000 \
    A2_G04=2-sample/A2_G04.fa  A2_G04_SIZE=190000 \
    A2_G30=2-sample/A2_G30.fa  A2_G30_SIZE=190000 \
    A2_G31=2-sample/A2_G31.fa  A2_G31_SIZE=190000 \
    A2_G32=2-sample/A2_G32.fa  A2_G32_SIZE=190000 \
    A2_G33=2-sample/A2_G33.fa  A2_G33_SIZE=190000 \
    A2_G34=2-sample/A2_G34.fa  A2_G34_SIZE=190000 \
    A2_G35=2-sample/A2_G35.fa  A2_G35_SIZE=190000 \
    A2_G36=2-sample/A2_G36.fa  A2_G36_SIZE=190000 \
    A2_G37=2-sample/A2_G37.fa  A2_G37_SIZE=190000 \
    A2_G38=2-sample/A2_G38.fa  A2_G38_SIZE=190000 \
    A2_G39=2-sample/A2_G39.fa  A2_G39_SIZE=190000 \
    A2_G40=2-sample/A2_G40.fa  A2_G40_SIZE=190000 \
    A2_G41=2-sample/A2_G41.fa  A2_G41_SIZE=190000 \
    A2_G42=2-sample/A2_G42.fa  A2_G42_SIZE=190000 \
    A2_G43=2-sample/A2_G43.fa  A2_G43_SIZE=190000 \
    A2_G44=2-sample/A2_G44.fa  A2_G44_SIZE=190000 \
    A2_G45=2-sample/A2_G45.fa  A2_G45_SIZE=190000 \
    A2_G46=2-sample/A2_G46.fa  A2_G46_SIZE=190000 \
    A2_G47=2-sample/A2_G47.fa  A2_G47_SIZE=190000 \
    A2_G48=2-sample/A2_G48.fa  A2_G48_SIZE=190000 \
    A2_G49=2-sample/A2_G49.fa  A2_G49_SIZE=190000 \
    A2_G50=2-sample/A2_G50.fa  A2_G50_SIZE=190000 \
    A2_G51=2-sample/A2_G51.fa  A2_G51_SIZE=190000 \
    A2_G52=2-sample/A2_G52.fa  A2_G52_SIZE=190000 \
    A2_G53=2-sample/A2_G53.fa  A2_G53_SIZE=190000 \
    A2_G54=2-sample/A2_G54.fa  A2_G54_SIZE=190000 \
    A2_G55=2-sample/A2_G55.fa  A2_G55_SIZE=190000 \
    A2_G56=2-sample/A2_G56.fa  A2_G56_SIZE=190000 \
    A2_G57=2-sample/A2_G57.fa  A2_G57_SIZE=190000 \
    A2_G58=2-sample/A2_G58.fa  A2_G58_SIZE=190000 \
    A2_G59=2-sample/A2_G59.fa  A2_G59_SIZE=190000 \
    A2_G60=2-sample/A2_G60.fa  A2_G60_SIZE=190000 \
    A2_G61=2-sample/A2_G61.fa  A2_G61_SIZE=190000 \
    A2_G62=2-sample/A2_G62.fa  A2_G62_SIZE=190000 \
    A2_G63=2-sample/A2_G63.fa  A2_G63_SIZE=190000 \
    A2_G64=2-sample/A2_G64.fa  A2_G64_SIZE=190000 \
    A2_G65=2-sample/A2_G65.fa  A2_G65_SIZE=190000 \
    A2_G66=2-sample/A2_G66.fa  A2_G66_SIZE=190000 \
    A2_N03=2-sample/A2_N03.fa  A2_N03_SIZE=190000 \
    A2_lin=2-sample/A2_lin.fa  A2_lin_SIZE=190000 \
    F1_lon=2-sample/F1_lon.fa  F1_lon_SIZE=145667 \
    $@