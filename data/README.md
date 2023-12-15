# Synechococcus elongatus Data Processing

The files used within the `build.py` script can be found as follows:
```
data
└── raw_data
   ├── Selon_omics_Pavlo
   │  ├── Description of metabolomics samples.docx
   │  ├── Metabolomics
   │  │  ├── Se_Rt_mixed_CC_EMSL_rem_PO4_HCO3_and_zer.csv
   │  │  ├── Se_Rt_mixed_CC_JHU.csv
   │  │  └── Se_Rt_mixed_CC_metadata.csv
   │  ├── Sample_names.xlsx
   │  └── Transcript
   │     ├── All genes vs model genes in mixed and mBR.pptx
   │     ├── Cyano_Reference.fasta
   │     ├── Cyano_Reference.gff3
   │     └── Se_mixedCC_raw.xlsx
   └── Selon_omics_Pavlo.zip
```

For the BMCA workflow, we want to incorporate files that highlight `Transcriptomics` and `Metabolomics` from the raw data and `Fluxomics` from tools such as EFlux2.

## Transcriptomics

The source transcriptomics file can be found at `data/raw_data/Selon_omics_Pavlo/Transcript/Se_mixedCC_raw.xlsx`. This file contains the structure:

```
                    Se_Rt_CC_90_1  Se_Rt_CC_90_2  Se_Rt_CC_90_3  ...  Se_ax_190_1  Se_ax_190_2  Se_ax_190_3
Label                                                            ...                                       
SYNPCC7942_RS00005          10732           9114          10103  ...         7016         8855         8689
SYNPCC7942_RS00010           5481           3735           4422  ...         3747         6670         6960
SYNPCC7942_RS00015          18664          17223          17564  ...        13440        17253        17053
SYNPCC7942_RS00020          19609          21833          26411  ...        24835        39769        30280
SYNPCC7942_RS00025           5652           4231           4710  ...         3326         4146         3180
...                           ...            ...            ...  ...          ...          ...          ...
HTX97_RS00025                 354            278            342  ...          302          350          354
HTX97_RS00030                1016            782            992  ...          606         1046          996
HTX97_RS00035                  10             10              8  ...            6            8           10
HTX97_RS00040                 226            136            180  ...          154          182          216
HTX97_RS00005                 138            120            114  ...           72          144          122
```

## Metabolomics

The source metabolomics file can be found at `data/raw_data/Selon_omics_Pavlo/Metabolomics/Se_Rt_mixed_CC_EMSL_rem_PO4_HCO3_and_zer.csv` and `data/raw_data/Selon_omics_Pavlo/Metabolomics/Se_Rt_mixed_CC_JHU.csv`. These files contain the structure:

```
          Se_Rt_cc_d1_1  Se_Rt_cc_d1_2  Se_Rt_cc_d1_3  Se_Rt_cc_d2_1  ...  Se_axen_d8_3  Se_axen_d9_1  Se_axen_d9_2  Se_axen_d9_3
Sample                                                                ...                                                        
114)PYR          123575         142656         110177          64895  ...        207550        245882        253436        174657
115)LACA           5617           6200           6400           6995  ...         20469         15861         65168         16261
116)ALA          362804         346109         281133         151103  ...        545044        629706        695662        477540
117)GLY          252244         234559         182766          82099  ...        442981        470467        564856        362018
118)VAL          207291         210857         168531          85338  ...        350909        391540        452276        319635
119)LEU          298167         294771         237533         123703  ...        431127        561591        627749        459052
120)ILE          558263         559182         476994         255027  ...        685940        788712        926085        739087
```
