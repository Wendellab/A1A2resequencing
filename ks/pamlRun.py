#!/bin/python2.7
import sys
import subprocess

def dS_permutations(w,p,filtLevel):
    fList = ['F1_SRR617704.F1.sort.bam__1','F1_SRR617704.F1.sort.bam__2']
    a1List = ['A1_012.F1.sort.bam__1','A1_012.F1.sort.bam__2','A1_019.F1.sort.bam__1','A1_019.F1.sort.bam__2','A1_028.F1.sort.bam__1','A1_028.F1.sort.bam__2','A1_029.F1.sort.bam__1','A1_029.F1.sort.bam__2','A1_030.F1.sort.bam__1','A1_030.F1.sort.bam__2','A1_037.F1.sort.bam__1','A1_037.F1.sort.bam__2','A1_051.F1.sort.bam__1','A1_051.F1.sort.bam__2','A1_054.F1.sort.bam__1','A1_054.F1.sort.bam__2','A1_073.F1.sort.bam__1','A1_073.F1.sort.bam__2','A1_074.F1.sort.bam__1','A1_074.F1.sort.bam__2','A1_077.F1.sort.bam__1','A1_077.F1.sort.bam__2','A1_079.F1.sort.bam__1','A1_079.F1.sort.bam__2','A1_084.F1.sort.bam__1','A1_084.F1.sort.bam__2','A1_088.F1.sort.bam__1','A1_088.F1.sort.bam__2','A1_097.F1.sort.bam__1','A1_097.F1.sort.bam__2','A1_108.F1.sort.bam__1','A1_108.F1.sort.bam__2','A1_113.F1.sort.bam__1','A1_113.F1.sort.bam__2','A1_123.F1.sort.bam__1','A1_123.F1.sort.bam__2','A1_125.F1.sort.bam__1','A1_125.F1.sort.bam__2','A1_132.F1.sort.bam__1','A1_132.F1.sort.bam__2','A1_133.F1.sort.bam__1','A1_133.F1.sort.bam__2','A1_148.F1.sort.bam__1','A1_148.F1.sort.bam__2','A1_155.F1.sort.bam__1','A1_155.F1.sort.bam__2','A1_158.F1.sort.bam__1','A1_158.F1.sort.bam__2','A1_Af.F1.sort.bam__1','A1_Af.F1.sort.bam__2','A1_Nisa.F1.sort.bam__1','A1_Nisa.F1.sort.bam__2','A1_SRR8283996_Ghe01.F1.sort.bam__1','A1_SRR8283996_Ghe01.F1.sort.bam__2','A1_SRR8283998_Ghe02.F1.sort.bam__1','A1_SRR8283998_Ghe02.F1.sort.bam__2','SRR4456927.F1.sort.bam__1','SRR4456927.F1.sort.bam__2','SRR4456928.F1.sort.bam__1','SRR4456928.F1.sort.bam__2','SRR4456952.F1.sort.bam__1','SRR4456952.F1.sort.bam__2','SRR4456953.F1.sort.bam__1','SRR4456953.F1.sort.bam__2','SRR4456954.F1.sort.bam__1','SRR4456954.F1.sort.bam__2','SRR4456955.F1.sort.bam__1','SRR4456955.F1.sort.bam__2','SRR4456956.F1.sort.bam__1','SRR4456956.F1.sort.bam__2','SRR4456957.F1.sort.bam__1','SRR4456957.F1.sort.bam__2','SRR4456958.F1.sort.bam__1','SRR4456958.F1.sort.bam__2','SRR4456959.F1.sort.bam__1','SRR4456959.F1.sort.bam__2','SRR4456961.F1.sort.bam__1','SRR4456961.F1.sort.bam__2','SRR4456962.F1.sort.bam__1','SRR4456962.F1.sort.bam__2','SRR4456979.F1.sort.bam__1','SRR4456979.F1.sort.bam__2','SRR4456980.F1.sort.bam__1','SRR4456980.F1.sort.bam__2','SRR4456981.F1.sort.bam__1','SRR4456981.F1.sort.bam__2']
    a2List = ['A2_001.F1.sort.bam__1','A2_001.F1.sort.bam__2','A2_004.F1.sort.bam__1','A2_004.F1.sort.bam__2','A2_008.F1.sort.bam__1','A2_008.F1.sort.bam__2','A2_016.F1.sort.bam__1','A2_016.F1.sort.bam__2','A2_019.F1.sort.bam__1','A2_019.F1.sort.bam__2','A2_021.F1.sort.bam__1','A2_021.F1.sort.bam__2','A2_026.F1.sort.bam__1','A2_026.F1.sort.bam__2','A2_029.F1.sort.bam__1','A2_029.F1.sort.bam__2','A2_030.F1.sort.bam__1','A2_030.F1.sort.bam__2','A2_031.F1.sort.bam__1','A2_031.F1.sort.bam__2','A2_034.F1.sort.bam__1','A2_034.F1.sort.bam__2','A2_038.F1.sort.bam__1','A2_038.F1.sort.bam__2','A2_040.F1.sort.bam__1','A2_040.F1.sort.bam__2','A2_041.F1.sort.bam__1','A2_041.F1.sort.bam__2','A2_043.F1.sort.bam__1','A2_043.F1.sort.bam__2','A2_044.F1.sort.bam__1','A2_044.F1.sort.bam__2','A2_045.F1.sort.bam__1','A2_045.F1.sort.bam__2','A2_046.F1.sort.bam__1','A2_046.F1.sort.bam__2','A2_047.F1.sort.bam__1','A2_047.F1.sort.bam__2','A2_056.F1.sort.bam__1','A2_056.F1.sort.bam__2','A2_057.F1.sort.bam__1','A2_057.F1.sort.bam__2','A2_059.F1.sort.bam__1','A2_059.F1.sort.bam__2','A2_060.F1.sort.bam__1','A2_060.F1.sort.bam__2','A2_061.F1.sort.bam__1','A2_061.F1.sort.bam__2','A2_062.F1.sort.bam__1','A2_062.F1.sort.bam__2','A2_064.F1.sort.bam__1','A2_064.F1.sort.bam__2','A2_065.F1.sort.bam__1','A2_065.F1.sort.bam__2','A2_066.F1.sort.bam__1','A2_066.F1.sort.bam__2','A2_067.F1.sort.bam__1','A2_067.F1.sort.bam__2','A2_069.F1.sort.bam__1','A2_069.F1.sort.bam__2','A2_073.F1.sort.bam__1','A2_073.F1.sort.bam__2','A2_074.F1.sort.bam__1','A2_074.F1.sort.bam__2','A2_076.F1.sort.bam__1','A2_076.F1.sort.bam__2','A2_079.F1.sort.bam__1','A2_079.F1.sort.bam__2','A2_080.F1.sort.bam__1','A2_080.F1.sort.bam__2','A2_084.F1.sort.bam__1','A2_084.F1.sort.bam__2','A2_085.F1.sort.bam__1','A2_085.F1.sort.bam__2','A2_087.F1.sort.bam__1','A2_087.F1.sort.bam__2','A2_091.F1.sort.bam__1','A2_091.F1.sort.bam__2','A2_096.F1.sort.bam__1','A2_096.F1.sort.bam__2','A2_100.F1.sort.bam__1','A2_100.F1.sort.bam__2','A2_101.F1.sort.bam__1','A2_101.F1.sort.bam__2','A2_113.F1.sort.bam__1','A2_113.F1.sort.bam__2','A2_117.F1.sort.bam__1','A2_117.F1.sort.bam__2','A2_118.F1.sort.bam__1','A2_118.F1.sort.bam__2','A2_119.F1.sort.bam__1','A2_119.F1.sort.bam__2','A2_123.F1.sort.bam__1','A2_123.F1.sort.bam__2','A2_124.F1.sort.bam__1','A2_124.F1.sort.bam__2','A2_131.F1.sort.bam__1','A2_131.F1.sort.bam__2','A2_141.F1.sort.bam__1','A2_141.F1.sort.bam__2','A2_142.F1.sort.bam__1','A2_142.F1.sort.bam__2','A2_147.F1.sort.bam__1','A2_147.F1.sort.bam__2','A2_154.F1.sort.bam__1','A2_154.F1.sort.bam__2','A2_156.F1.sort.bam__1','A2_156.F1.sort.bam__2','A2_158.F1.sort.bam__1','A2_158.F1.sort.bam__2','A2_164.F1.sort.bam__1','A2_164.F1.sort.bam__2','A2_255.F1.sort.bam__1','A2_255.F1.sort.bam__2','A2_lintless.F1.sort.bam__1','A2_lintless.F1.sort.bam__2','A2_NC503.F1.sort.bam__1','A2_NC503.F1.sort.bam__2','A2_SRR8283969_Gar42.F1.sort.bam__1','A2_SRR8283969_Gar42.F1.sort.bam__2','A2_SRR8283970_Gar43.F1.sort.bam__1','A2_SRR8283970_Gar43.F1.sort.bam__2','A2_SRR8283971_Gar40.F1.sort.bam__1','A2_SRR8283971_Gar40.F1.sort.bam__2','A2_SRR8283972_Gar41.F1.sort.bam__1','A2_SRR8283972_Gar41.F1.sort.bam__2','A2_SRR8283973_Gar38.F1.sort.bam__1','A2_SRR8283973_Gar38.F1.sort.bam__2','A2_SRR8283974_Gar39.F1.sort.bam__1','A2_SRR8283974_Gar39.F1.sort.bam__2','A2_SRR8283975_Gar36.F1.sort.bam__1','A2_SRR8283975_Gar36.F1.sort.bam__2','A2_SRR8283976_Gar37.F1.sort.bam__1','A2_SRR8283976_Gar37.F1.sort.bam__2','A2_SRR8283977_Gar56.F1.sort.bam__1','A2_SRR8283977_Gar56.F1.sort.bam__2','A2_SRR8283978_Gar57.F1.sort.bam__1','A2_SRR8283978_Gar57.F1.sort.bam__2','A2_SRR8283979_Gar58.F1.sort.bam__1','A2_SRR8283979_Gar58.F1.sort.bam__2','A2_SRR8283980_Gar59.F1.sort.bam__1','A2_SRR8283980_Gar59.F1.sort.bam__2','A2_SRR8283981_Gar60.F1.sort.bam__1','A2_SRR8283981_Gar60.F1.sort.bam__2','A2_SRR8283982_Gar61.F1.sort.bam__1','A2_SRR8283982_Gar61.F1.sort.bam__2','A2_SRR8283983_Gar44.F1.sort.bam__1','A2_SRR8283983_Gar44.F1.sort.bam__2','A2_SRR8283984_Gar45.F1.sort.bam__1','A2_SRR8283984_Gar45.F1.sort.bam__2','A2_SRR8283985_Gar34.F1.sort.bam__1','A2_SRR8283985_Gar34.F1.sort.bam__2','A2_SRR8283986_Gar35.F1.sort.bam__1','A2_SRR8283986_Gar35.F1.sort.bam__2','A2_SRR8283987_Gar32.F1.sort.bam__1','A2_SRR8283987_Gar32.F1.sort.bam__2','A2_SRR8283988_Gar33.F1.sort.bam__1','A2_SRR8283988_Gar33.F1.sort.bam__2','A2_SRR8283989_Gar30.F1.sort.bam__1','A2_SRR8283989_Gar30.F1.sort.bam__2','A2_SRR8283990_Gar31.F1.sort.bam__1','A2_SRR8283990_Gar31.F1.sort.bam__2','A2_SRR8283991_Gar03.F1.sort.bam__1','A2_SRR8283991_Gar03.F1.sort.bam__2','A2_SRR8283992_Gar04.F1.sort.bam__1','A2_SRR8283992_Gar04.F1.sort.bam__2','A2_SRR8283993_Gar01.F1.sort.bam__1','A2_SRR8283993_Gar01.F1.sort.bam__2','A2_SRR8283994_Gar02.F1.sort.bam__1','A2_SRR8283994_Gar02.F1.sort.bam__2','A2_SRR8283995_Gar65.F1.sort.bam__1','A2_SRR8283995_Gar65.F1.sort.bam__2','A2_SRR8283997_Gar66.F1.sort.bam__1','A2_SRR8283997_Gar66.F1.sort.bam__2','A2_SRR8283999_Gar55.F1.sort.bam__1','A2_SRR8283999_Gar55.F1.sort.bam__2','A2_SRR8284000_Gar54.F1.sort.bam__1','A2_SRR8284000_Gar54.F1.sort.bam__2','A2_SRR8284001_Gar53.F1.sort.bam__1','A2_SRR8284001_Gar53.F1.sort.bam__2','A2_SRR8284002_Gar52.F1.sort.bam__1','A2_SRR8284002_Gar52.F1.sort.bam__2','A2_SRR8284003_Gar51.F1.sort.bam__1','A2_SRR8284003_Gar51.F1.sort.bam__2','A2_SRR8284004_Gar50.F1.sort.bam__1','A2_SRR8284004_Gar50.F1.sort.bam__2','A2_SRR8284005_Gar49.F1.sort.bam__1','A2_SRR8284005_Gar49.F1.sort.bam__2','A2_SRR8284006_Gar48.F1.sort.bam__1','A2_SRR8284006_Gar48.F1.sort.bam__2','A2_SRR8284007_Gar47.F1.sort.bam__1','A2_SRR8284007_Gar47.F1.sort.bam__2','A2_SRR8284008_Gar46.F1.sort.bam__1','A2_SRR8284008_Gar46.F1.sort.bam__2','A2_SRR8284009_Gar62.F1.sort.bam__1','A2_SRR8284009_Gar62.F1.sort.bam__2','A2_SRR8284010_Gar63.F1.sort.bam__1','A2_SRR8284010_Gar63.F1.sort.bam__2','A2_SRR8284011_Gar64.F1.sort.bam__1','A2_SRR8284011_Gar64.F1.sort.bam__2','SRR4456872.F1.sort.bam__1','SRR4456872.F1.sort.bam__2','SRR4456873.F1.sort.bam__1','SRR4456873.F1.sort.bam__2','SRR4456874.F1.sort.bam__1','SRR4456874.F1.sort.bam__2','SRR4456875.F1.sort.bam__1','SRR4456875.F1.sort.bam__2','SRR4456876.F1.sort.bam__1','SRR4456876.F1.sort.bam__2','SRR4456877.F1.sort.bam__1','SRR4456877.F1.sort.bam__2','SRR4456878.F1.sort.bam__1','SRR4456878.F1.sort.bam__2','SRR4456879.F1.sort.bam__1','SRR4456879.F1.sort.bam__2','SRR4456880.F1.sort.bam__1','SRR4456880.F1.sort.bam__2','SRR4456881.F1.sort.bam__1','SRR4456881.F1.sort.bam__2','SRR4456882.F1.sort.bam__1','SRR4456882.F1.sort.bam__2','SRR4456883.F1.sort.bam__1','SRR4456883.F1.sort.bam__2','SRR4456884.F1.sort.bam__1','SRR4456884.F1.sort.bam__2','SRR4456885.F1.sort.bam__1','SRR4456885.F1.sort.bam__2','SRR4456886.F1.sort.bam__1','SRR4456886.F1.sort.bam__2','SRR4456887.F1.sort.bam__1','SRR4456887.F1.sort.bam__2','SRR4456888.F1.sort.bam__1','SRR4456888.F1.sort.bam__2','SRR4456889.F1.sort.bam__1','SRR4456889.F1.sort.bam__2','SRR4456890.F1.sort.bam__1','SRR4456890.F1.sort.bam__2','SRR4456891.F1.sort.bam__1','SRR4456891.F1.sort.bam__2','SRR4456892.F1.sort.bam__1','SRR4456892.F1.sort.bam__2','SRR4456893.F1.sort.bam__1','SRR4456893.F1.sort.bam__2','SRR4456894.F1.sort.bam__1','SRR4456894.F1.sort.bam__2','SRR4456895.F1.sort.bam__1','SRR4456895.F1.sort.bam__2','SRR4456896.F1.sort.bam__1','SRR4456896.F1.sort.bam__2','SRR4456897.F1.sort.bam__1','SRR4456897.F1.sort.bam__2','SRR4456898.F1.sort.bam__1','SRR4456898.F1.sort.bam__2','SRR4456899.F1.sort.bam__1','SRR4456899.F1.sort.bam__2','SRR4456900.F1.sort.bam__1','SRR4456900.F1.sort.bam__2','SRR4456901.F1.sort.bam__1','SRR4456901.F1.sort.bam__2','SRR4456902.F1.sort.bam__1','SRR4456902.F1.sort.bam__2','SRR4456903.F1.sort.bam__1','SRR4456903.F1.sort.bam__2','SRR4456904.F1.sort.bam__1','SRR4456904.F1.sort.bam__2','SRR4456905.F1.sort.bam__1','SRR4456905.F1.sort.bam__2','SRR4456906.F1.sort.bam__1','SRR4456906.F1.sort.bam__2','SRR4456907.F1.sort.bam__1','SRR4456907.F1.sort.bam__2','SRR4456908.F1.sort.bam__1','SRR4456908.F1.sort.bam__2','SRR4456909.F1.sort.bam__1','SRR4456909.F1.sort.bam__2','SRR4456910.F1.sort.bam__1','SRR4456910.F1.sort.bam__2','SRR4456911.F1.sort.bam__1','SRR4456911.F1.sort.bam__2','SRR4456912.F1.sort.bam__1','SRR4456912.F1.sort.bam__2','SRR4456913.F1.sort.bam__1','SRR4456913.F1.sort.bam__2','SRR4456914.F1.sort.bam__1','SRR4456914.F1.sort.bam__2','SRR4456915.F1.sort.bam__1','SRR4456915.F1.sort.bam__2','SRR4456916.F1.sort.bam__1','SRR4456916.F1.sort.bam__2','SRR4456917.F1.sort.bam__1','SRR4456917.F1.sort.bam__2','SRR4456918.F1.sort.bam__1','SRR4456918.F1.sort.bam__2','SRR4456919.F1.sort.bam__1','SRR4456919.F1.sort.bam__2','SRR4456920.F1.sort.bam__1','SRR4456920.F1.sort.bam__2','SRR4456921.F1.sort.bam__1','SRR4456921.F1.sort.bam__2','SRR4456922.F1.sort.bam__1','SRR4456922.F1.sort.bam__2','SRR4456923.F1.sort.bam__1','SRR4456923.F1.sort.bam__2','SRR4456924.F1.sort.bam__1','SRR4456924.F1.sort.bam__2','SRR4456925.F1.sort.bam__1','SRR4456925.F1.sort.bam__2','SRR4456926.F1.sort.bam__1','SRR4456926.F1.sort.bam__2','SRR4456929.F1.sort.bam__1','SRR4456929.F1.sort.bam__2','SRR4456930.F1.sort.bam__1','SRR4456930.F1.sort.bam__2','SRR4456931.F1.sort.bam__1','SRR4456931.F1.sort.bam__2','SRR4456932.F1.sort.bam__1','SRR4456932.F1.sort.bam__2','SRR4456933.F1.sort.bam__1','SRR4456933.F1.sort.bam__2','SRR4456934.F1.sort.bam__1','SRR4456934.F1.sort.bam__2','SRR4456935.F1.sort.bam__1','SRR4456935.F1.sort.bam__2','SRR4456936.F1.sort.bam__1','SRR4456936.F1.sort.bam__2','SRR4456937.F1.sort.bam__1','SRR4456937.F1.sort.bam__2','SRR4456938.F1.sort.bam__1','SRR4456938.F1.sort.bam__2','SRR4456939.F1.sort.bam__1','SRR4456939.F1.sort.bam__2','SRR4456940.F1.sort.bam__1','SRR4456940.F1.sort.bam__2','SRR4456941.F1.sort.bam__1','SRR4456941.F1.sort.bam__2','SRR4456942.F1.sort.bam__1','SRR4456942.F1.sort.bam__2','SRR4456943.F1.sort.bam__1','SRR4456943.F1.sort.bam__2','SRR4456944.F1.sort.bam__1','SRR4456944.F1.sort.bam__2','SRR4456945.F1.sort.bam__1','SRR4456945.F1.sort.bam__2','SRR4456946.F1.sort.bam__1','SRR4456946.F1.sort.bam__2','SRR4456947.F1.sort.bam__1','SRR4456947.F1.sort.bam__2','SRR4456948.F1.sort.bam__1','SRR4456948.F1.sort.bam__2','SRR4456949.F1.sort.bam__1','SRR4456949.F1.sort.bam__2','SRR4456950.F1.sort.bam__1','SRR4456950.F1.sort.bam__2','SRR4456951.F1.sort.bam__1','SRR4456951.F1.sort.bam__2','SRR4456960.F1.sort.bam__1','SRR4456960.F1.sort.bam__2','SRR4456963.F1.sort.bam__1','SRR4456963.F1.sort.bam__2','SRR4456964.F1.sort.bam__1','SRR4456964.F1.sort.bam__2','SRR4456965.F1.sort.bam__1','SRR4456965.F1.sort.bam__2','SRR4456966.F1.sort.bam__1','SRR4456966.F1.sort.bam__2','SRR4456967.F1.sort.bam__1','SRR4456967.F1.sort.bam__2','SRR4456968.F1.sort.bam__1','SRR4456968.F1.sort.bam__2','SRR4456969.F1.sort.bam__1','SRR4456969.F1.sort.bam__2','SRR4456970.F1.sort.bam__1','SRR4456970.F1.sort.bam__2','SRR4456971.F1.sort.bam__1','SRR4456971.F1.sort.bam__2','SRR4456972.F1.sort.bam__1','SRR4456972.F1.sort.bam__2','SRR4456973.F1.sort.bam__1','SRR4456973.F1.sort.bam__2','SRR4456974.F1.sort.bam__1','SRR4456974.F1.sort.bam__2','SRR4456975.F1.sort.bam__1','SRR4456975.F1.sort.bam__2','SRR4456976.F1.sort.bam__1','SRR4456976.F1.sort.bam__2','SRR4456977.F1.sort.bam__1','SRR4456977.F1.sort.bam__2','SRR4456978.F1.sort.bam__1','SRR4456978.F1.sort.bam__2','SRR4456982.F1.sort.bam__1','SRR4456982.F1.sort.bam__2','SRR4456983.F1.sort.bam__1','SRR4456983.F1.sort.bam__2','SRR4456984.F1.sort.bam__1','SRR4456984.F1.sort.bam__2','SRR4456985.F1.sort.bam__1','SRR4456985.F1.sort.bam__2','SRR4456986.F1.sort.bam__1','SRR4456986.F1.sort.bam__2','SRR4456987.F1.sort.bam__1','SRR4456987.F1.sort.bam__2','SRR4456988.F1.sort.bam__1','SRR4456988.F1.sort.bam__2','SRR4456989.F1.sort.bam__1','SRR4456989.F1.sort.bam__2','SRR4456990.F1.sort.bam__1','SRR4456990.F1.sort.bam__2','SRR4456991.F1.sort.bam__1','SRR4456991.F1.sort.bam__2','SRR4456992.F1.sort.bam__1','SRR4456992.F1.sort.bam__2','SRR4456993.F1.sort.bam__1','SRR4456993.F1.sort.bam__2','SRR4456994.F1.sort.bam__1','SRR4456994.F1.sort.bam__2','SRR4456995.F1.sort.bam__1','SRR4456995.F1.sort.bam__2','SRR4456996.F1.sort.bam__1','SRR4456996.F1.sort.bam__2','SRR4456997.F1.sort.bam__1','SRR4456997.F1.sort.bam__2','SRR4456998.F1.sort.bam__1','SRR4456998.F1.sort.bam__2','SRR4456999.F1.sort.bam__1','SRR4456999.F1.sort.bam__2','SRR4457000.F1.sort.bam__1','SRR4457000.F1.sort.bam__2','SRR4457001.F1.sort.bam__1','SRR4457001.F1.sort.bam__2','SRR4457002.F1.sort.bam__1','SRR4457002.F1.sort.bam__2','SRR4457003.F1.sort.bam__1','SRR4457003.F1.sort.bam__2','SRR4457004.F1.sort.bam__1','SRR4457004.F1.sort.bam__2','SRR4457005.F1.sort.bam__1','SRR4457005.F1.sort.bam__2','SRR4457006.F1.sort.bam__1','SRR4457006.F1.sort.bam__2','SRR4457007.F1.sort.bam__1','SRR4457007.F1.sort.bam__2','SRR4457008.F1.sort.bam__1','SRR4457008.F1.sort.bam__2','SRR4457009.F1.sort.bam__1','SRR4457009.F1.sort.bam__2','SRR4457010.F1.sort.bam__1','SRR4457010.F1.sort.bam__2','SRR4457011.F1.sort.bam__1','SRR4457011.F1.sort.bam__2','SRR4457012.F1.sort.bam__1','SRR4457012.F1.sort.bam__2','SRR4457013.F1.sort.bam__1','SRR4457013.F1.sort.bam__2','SRR4457014.F1.sort.bam__1','SRR4457014.F1.sort.bam__2','SRR4457015.F1.sort.bam__1','SRR4457015.F1.sort.bam__2','SRR4457016.F1.sort.bam__1','SRR4457016.F1.sort.bam__2','SRR4457017.F1.sort.bam__1','SRR4457017.F1.sort.bam__2','SRR4457018.F1.sort.bam__1','SRR4457018.F1.sort.bam__2','SRR4457019.F1.sort.bam__1','SRR4457019.F1.sort.bam__2','SRR4457020.F1.sort.bam__1','SRR4457020.F1.sort.bam__2','SRR4457021.F1.sort.bam__1','SRR4457021.F1.sort.bam__2','SRR4457022.F1.sort.bam__1','SRR4457022.F1.sort.bam__2','SRR4457023.F1.sort.bam__1','SRR4457023.F1.sort.bam__2','SRR4457024.F1.sort.bam__1','SRR4457024.F1.sort.bam__2','SRR4457025.F1.sort.bam__1','SRR4457025.F1.sort.bam__2','SRR4457026.F1.sort.bam__1','SRR4457026.F1.sort.bam__2','SRR4457027.F1.sort.bam__1','SRR4457027.F1.sort.bam__2','SRR4457028.F1.sort.bam__1','SRR4457028.F1.sort.bam__2','SRR4457029.F1.sort.bam__1','SRR4457029.F1.sort.bam__2','SRR4457030.F1.sort.bam__1','SRR4457030.F1.sort.bam__2','SRR4457031.F1.sort.bam__1','SRR4457031.F1.sort.bam__2','SRR4457032.F1.sort.bam__1','SRR4457032.F1.sort.bam__2','SRR4457033.F1.sort.bam__1','SRR4457033.F1.sort.bam__2','SRR4457034.F1.sort.bam__1','SRR4457034.F1.sort.bam__2','SRR4457035.F1.sort.bam__1','SRR4457035.F1.sort.bam__2','SRR4457036.F1.sort.bam__1','SRR4457036.F1.sort.bam__2','SRR4457037.F1.sort.bam__1','SRR4457037.F1.sort.bam__2','SRR4457038.F1.sort.bam__1','SRR4457038.F1.sort.bam__2','SRR4457039.F1.sort.bam__1','SRR4457039.F1.sort.bam__2','SRR4457040.F1.sort.bam__1','SRR4457040.F1.sort.bam__2','SRR4457041.F1.sort.bam__1','SRR4457041.F1.sort.bam__2','SRR4457042.F1.sort.bam__1','SRR4457042.F1.sort.bam__2','SRR4457043.F1.sort.bam__1','SRR4457043.F1.sort.bam__2','SRR4457044.F1.sort.bam__1','SRR4457044.F1.sort.bam__2','SRR4457045.F1.sort.bam__1','SRR4457045.F1.sort.bam__2','SRR4457046.F1.sort.bam__1','SRR4457046.F1.sort.bam__2','SRR4457047.F1.sort.bam__1','SRR4457047.F1.sort.bam__2','SRR4457048.F1.sort.bam__1','SRR4457048.F1.sort.bam__2','SRR4457049.F1.sort.bam__1','SRR4457049.F1.sort.bam__2','SRR4457050.F1.sort.bam__1','SRR4457050.F1.sort.bam__2','SRR4457051.F1.sort.bam__1','SRR4457051.F1.sort.bam__2','SRR4457052.F1.sort.bam__1','SRR4457052.F1.sort.bam__2','SRR4457053.F1.sort.bam__1','SRR4457053.F1.sort.bam__2','SRR4457054.F1.sort.bam__1','SRR4457054.F1.sort.bam__2','SRR4457055.F1.sort.bam__1','SRR4457055.F1.sort.bam__2','SRR4457056.F1.sort.bam__1','SRR4457056.F1.sort.bam__2','SRR4457057.F1.sort.bam__1','SRR4457057.F1.sort.bam__2','SRR4457058.F1.sort.bam__1','SRR4457058.F1.sort.bam__2','SRR4457059.F1.sort.bam__1','SRR4457059.F1.sort.bam__2','SRR4457060.F1.sort.bam__1','SRR4457060.F1.sort.bam__2','SRR4457061.F1.sort.bam__1','SRR4457061.F1.sort.bam__2','SRR4457062.F1.sort.bam__1','SRR4457062.F1.sort.bam__2','SRR4457063.F1.sort.bam__1','SRR4457063.F1.sort.bam__2','SRR4457064.F1.sort.bam__1','SRR4457064.F1.sort.bam__2','SRR4457065.F1.sort.bam__1','SRR4457065.F1.sort.bam__2','SRR4457066.F1.sort.bam__1','SRR4457066.F1.sort.bam__2','SRR4457067.F1.sort.bam__1','SRR4457067.F1.sort.bam__2','SRR4457068.F1.sort.bam__1','SRR4457068.F1.sort.bam__2','SRR4457069.F1.sort.bam__1','SRR4457069.F1.sort.bam__2','SRR4457070.F1.sort.bam__1','SRR4457070.F1.sort.bam__2','SRR4457071.F1.sort.bam__1','SRR4457071.F1.sort.bam__2','SRR4457072.F1.sort.bam__1','SRR4457072.F1.sort.bam__2','SRR4457073.F1.sort.bam__1','SRR4457073.F1.sort.bam__2','SRR4457074.F1.sort.bam__1','SRR4457074.F1.sort.bam__2','SRR4457075.F1.sort.bam__1','SRR4457075.F1.sort.bam__2','SRR4457076.F1.sort.bam__1','SRR4457076.F1.sort.bam__2','SRR4457077.F1.sort.bam__1','SRR4457077.F1.sort.bam__2','SRR4457078.F1.sort.bam__1','SRR4457078.F1.sort.bam__2','SRR4457079.F1.sort.bam__1','SRR4457079.F1.sort.bam__2','SRR4457080.F1.sort.bam__1','SRR4457080.F1.sort.bam__2','SRR4457081.F1.sort.bam__1','SRR4457081.F1.sort.bam__2','SRR4457082.F1.sort.bam__1','SRR4457082.F1.sort.bam__2','SRR4457083.F1.sort.bam__1','SRR4457083.F1.sort.bam__2','SRR4457084.F1.sort.bam__1','SRR4457084.F1.sort.bam__2','SRR4457085.F1.sort.bam__1','SRR4457085.F1.sort.bam__2','SRR4457086.F1.sort.bam__1','SRR4457086.F1.sort.bam__2','SRR4457087.F1.sort.bam__1','SRR4457087.F1.sort.bam__2','SRR4457088.F1.sort.bam__1','SRR4457088.F1.sort.bam__2','SRR4457089.F1.sort.bam__1','SRR4457089.F1.sort.bam__2','SRR4457090.F1.sort.bam__1','SRR4457090.F1.sort.bam__2','SRR4457091.F1.sort.bam__1','SRR4457091.F1.sort.bam__2','SRR4457092.F1.sort.bam__1','SRR4457092.F1.sort.bam__2','SRR4457093.F1.sort.bam__1','SRR4457093.F1.sort.bam__2','SRR4457094.F1.sort.bam__1','SRR4457094.F1.sort.bam__2','SRR4457095.F1.sort.bam__1','SRR4457095.F1.sort.bam__2','SRR4457096.F1.sort.bam__1','SRR4457096.F1.sort.bam__2','SRR4457097.F1.sort.bam__1','SRR4457097.F1.sort.bam__2','SRR4457098.F1.sort.bam__1','SRR4457098.F1.sort.bam__2','SRR4457099.F1.sort.bam__1','SRR4457099.F1.sort.bam__2','SRR4457100.F1.sort.bam__1','SRR4457100.F1.sort.bam__2','SRR4457101.F1.sort.bam__1','SRR4457101.F1.sort.bam__2','SRR4457102.F1.sort.bam__1','SRR4457102.F1.sort.bam__2','SRR4457103.F1.sort.bam__1','SRR4457103.F1.sort.bam__2','SRR4457104.F1.sort.bam__1','SRR4457104.F1.sort.bam__2','SRR4457105.F1.sort.bam__1','SRR4457105.F1.sort.bam__2','SRR4457106.F1.sort.bam__1','SRR4457106.F1.sort.bam__2','SRR4457107.F1.sort.bam__1','SRR4457107.F1.sort.bam__2','SRR4457108.F1.sort.bam__1','SRR4457108.F1.sort.bam__2','SRR4457109.F1.sort.bam__1','SRR4457109.F1.sort.bam__2','SRR4457110.F1.sort.bam__1','SRR4457110.F1.sort.bam__2','SRR4457111.F1.sort.bam__1','SRR4457111.F1.sort.bam__2','SRR4457112.F1.sort.bam__1','SRR4457112.F1.sort.bam__2','SRR4457113.F1.sort.bam__1','SRR4457113.F1.sort.bam__2','SRR4457114.F1.sort.bam__1','SRR4457114.F1.sort.bam__2']
    permutationDict = {}
    i=0
    for f in fList:
        for a1 in a1List:
            for a2 in a2List:
                permutationDict[i] = (f,a1,a2)
                i += 1
    if p >= i:
        return
    ctlFile = open(str(w) + '.' + str(p) + '.ctl','w')
    ctlFile.write('      seqfile = ' + str(w) + '.' + str(p) + '.paml * sequence data filename\n     treefile = Agenome.trees' + '     * tree structure file name\n      outfile = ' + str(w) + '.' + str(p) + '.out           * main result file name\n\n        noisy = 0  * 0,1,2,3,9: how much rubbish on the screen\n      verbose = 0  * 0: concise; 1: detailed, 2: too much\n      runmode = 0  * 0: user tree;  1: semi-automatic;  2: automatic\n                   * 3: StepwiseAddition; (4,5):PerturbationNNI; -2: pairwise\n\n      seqtype = 1  * 1:codons; 2:AAs; 3:codons-->AAs\n    CodonFreq = 2  * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table\n\n*        ndata = 1\n        clock = 0  * 0:no clock, 1:clock; 2:local clock; 3:CombinedAnalysis\n       aaDist = 0  * 0:equal, +:geometric; -:linear, 1-6:G1974,Miyata,c,p,v,a\n   aaRatefile = dat/jones.dat  * only used for aa seqs with model=empirical(_F)\n                   * dayhoff.dat, jones.dat, wag.dat, mtmam.dat, or your own\n\n        model = 0\n                   * models for codons:\n                       * 0:one, 1:b, 2:2 or more dN/dS ratios for branches\n                   * models for AAs or codon-translated AAs:\n                       * 0:poisson, 1:proportional, 2:Empirical, 3:Empirical+F\n                       * 6:FromCodon, 7:AAClasses, 8:REVaa_0, 9:REVaa(nr=189)\n\n      NSsites = 0  * 0:one w;1:neutral;2:selection; 3:discrete;4:freqs;\n                   * 5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma;\n                   * 10:beta&gamma+1; 11:beta&normal>1; 12:0&2normal>1;\n                   * 13:3normal>0\n\n        icode = 0  * 0:universal code; 1:mammalian mt; 2-10:see below\n        Mgene = 0\n                   * codon: 0:rates, 1:separate; 2:diff pi, 3:diff kapa, 4:all diff\n                   * AA: 0:rates, 1:separate\n\n    fix_kappa = 0  * 1: kappa fixed, 0: kappa to be estimated\n        kappa = 2  * initial or fixed kappa\n    fix_omega = 0  * 1: omega or omega_1 fixed, 0: estimate \n        omega = .4 * initial or fixed omega, for codons or codon-based AAs\n        \n        fix_alpha = 1  * 0: estimate gamma shape parameter; 1: fix it at alpha\n        alpha = 0 * initial or fixed alpha, 0:infinity (constant rate)\n       Malpha = 0  * different alphas for genes\n        ncatG = 10  * # of categories in dG of NSsites models\n\n\n        getSE = 0  * 0: dont want them, 1: want S.E.s of estimates\n RateAncestor = 0  * (0,1,2): rates (alpha>0) or ancestral states (1 or 2)\n\n\n   Small_Diff = .5e-6\n    cleandata = 1  * remove sites with ambiguity data (1:yes, 0:no)?\n*  fix_blength = 0  * 0: ignore, -1: random, 1: initial, 2: fixed\n       method = 0  * Optimization method 0: simultaneous; 1: one branch a time\n\n\n* Genetic codes: 0:universal, 1:mammalian mt., 2:yeast mt., 3:mold mt.,\n* 4: invertebrate mt., 5: ciliate nuclear, 6: echinoderm mt., \n* 7: euplotid mt., 8: alternative yeast nu. 9: ascidian mt., \n* 10: blepharisma nu.\n* These codes correspond to transl_table 1 to 11 of GENEBANK.')
    ctlFile.close()
    fofn = open(filtLevel + '.fofn','r')
    lines = fofn.readlines()
    fofn.close()
    fh = lines[w]
    while fh[-1] == '\t' or fh[-1] == '\r' or fh[-1] == '\n':
        fh = fh[0:-1]
    seqDict = buildSeqDict(fh)
    currPermutation = permutationDict[p]
    f = currPermutation[0]
    a1 = currPermutation[1]
    a2 = currPermutation[2]
    pamlFile = open(str(w) + '.' + str(p) + '.paml','w')
    pamlFile.write(' 3 ' +  str(len(seqDict[f])) + '\nF\n' + seqDict[f] + '\nA1\n' + seqDict[a1] + '\nA2\n' + seqDict[a2] + '\n')
    pamlFile.close()
    cmd = subprocess.call(['./codeml', str(w) + '.' + str(p) + '.ctl'],stdout=subprocess.PIPE)
    fSplit = f.split('.')
    a1Split = a1.split('.')
    a2Split = a2.split('.')
    if f[-1] == '1':
        fHap = '1'
    else:
        fHap = '2'
    if a1[-1] == '1':
        a1Hap = '1'
    else:
        a1Hap = '2'
    if a2[-1] == '1':
        a2Hap = '1'
    else:
        a2Hap = '2'
    f = fSplit[0]
    a1 = a1Split[0]
    a2 = a2Split[0]
    tree = '(' + f + '('  + a1 + ', ' + a2 + ');'
    sys.stdout.write(str(w) + '.' + str(p) + '\t' + tree + '\t' + f + '\t' + fHap + '\t' + a1 + '\t' + a1Hap + '\t' + a2 + '\t' + a2Hap + '\t' + str(pairwiseDS(str(w) + '.' + str(p) + '.out')) + '\n')

def pairwiseDS(m0_pamlOut):
    dS = 'NA'
    infile = open(m0_pamlOut,'r')
    dS1=False
    dS2=False
    lines = infile.readlines()
    infile.close()
    lineNum = 0
    values=False
    for line in lines:
        if ' branch          t       N       S   dN/dS      dN      dS  N*dN  S*dS' in line:
            values=True
        elif values == True and '4..2' in line:
            realLine = line
            while realLine[-1] == '\t' or realLine[-1] == '\n' or realLine[-1] == '\r':
                realLine = realLine[0:-1]
            lineSplit = realLine.split(' ')
            while '' in lineSplit:
                lineSplit.remove('')
            dS1 = float(lineSplit[6])
        elif values == True and '4..3' in line:
            realLine = line
            while realLine[-1] == '\t' or realLine[-1] == '\n' or realLine[-1] == '\r':
                realLine = realLine[0:-1]
            lineSplit = realLine.split(' ')
            while '' in lineSplit:
                lineSplit.remove('')
            dS2 = float(lineSplit[6])
    if dS1 != False and dS2 != False:
        dS = dS1 + dS2
    return dS


def buildSeqDict(fasta):
    infile = open(fasta,'r')
    scaffoldDict = {}
    seqName = ''
    currSeq = ''
    for line in infile:
        if line[0] == '>':
            if seqName != '':
                scaffoldDict[seqName] = currSeq
            seqName = line[1:]
            while seqName[-1] == '\n' or seqName[-1] == '\t' or seqName[-1] == '\r':
                seqName = seqName[0:-1]
            currSeq = ''
        else:
            currSeq += line
            while currSeq[-1] == '\n' or currSeq[-1] == '\t' or currSeq[-1] == '\r':
                currSeq = currSeq[0:-1]
    scaffoldDict[seqName] = currSeq 
    return scaffoldDict

dS_permutations(int(sys.argv[3]),int(sys.argv[2]),sys.argv[1]) #w,p